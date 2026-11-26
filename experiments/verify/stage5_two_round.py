"""Section 6 of arXiv:2507.02181v4: exact two-round key recovery.

  E2(P) = L(S(L(S(P ^ K1)) ^ K2)) ^ K3,   K1||K2 = 256-bit master key.

Part A re-derives the paper's numbers (sigma_perp, det D* mod 5, Table 5, N1,
rho_2, N2). Part B runs the actual chosen-plaintext attack on a vectorised
two-round model that is first checked against the reference implementation.

Run from the repo root:  python3 experiments/verify/stage5_two_round.py [n_random_keys]
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf import POLY_K, mul_table

sys.path.insert(0, "src/gost")
import kuznechik as KZ

OUT = "results/verify"
Q = 256
C = 0x04
S = np.frombuffer(open("data/original/pi.bin", "rb").read(), dtype=np.uint8).astype(np.int64)
SINV = np.argsort(S)
MK = np.array(mul_table(POLY_K), dtype=np.int64)
X = np.arange(256)


def ext_table(c):
    N = np.zeros((256, 256), dtype=np.int64)
    for a in range(256):
        N[a] = np.bincount(S[MK[c, X] ^ a] ^ S[X], minlength=256)
    return N


def matrix_of(fn):
    """F256-matrix of an F256-linear map given on hex states: M[i, t] = fn(e_t)[i]."""
    M = np.zeros((16, 16), dtype=np.int64)
    for t in range(16):
        e = bytearray(16)
        e[t] = 1
        M[:, t] = list(bytes.fromhex(fn(e.hex())))
    return M


LM = matrix_of(KZ.l_layer_r16)
LINV = matrix_of(KZ.l_inv_layer)


def apply(M, V):
    """V: (n, 16) array of bytes -> M V over F256."""
    out = np.zeros_like(V)
    for i in range(16):
        for t in range(16):
            out[:, i] ^= MK[M[i, t], V[:, t]]
    return out


def det_mod_p(A, p):
    A = A.copy() % p
    n = len(A)
    det = 1
    for col in range(n):
        piv = next((r for r in range(col, n) if A[r, col]), None)
        if piv is None:
            return 0
        if piv != col:
            A[[col, piv]] = A[[piv, col]]
            det = -det
        det = det * A[col, col] % p
        inv = pow(int(A[col, col]), p - 2, p)
        A[col + 1:] = (A[col + 1:] - np.outer(A[col + 1:, col] * inv % p, A[col])) % p
    return det % p


def bhattacharyya_rho(P):
    R = np.sqrt(P)
    B = R @ R.T
    np.fill_diagonal(B, -1)
    return B.max()


def main():
    os.makedirs(OUT, exist_ok=True)
    out = {}
    # ---------- sanity: matrices and model ----------
    assert (apply(LINV, apply(LM, np.eye(16, dtype=np.int64))) == np.eye(16)).all()
    print("L is F256-linear with all 256 coefficients nonzero:", bool((LM != 0).all()))

    D = np.zeros((256, 256), dtype=np.int64)
    for d in range(256):
        D[d] = np.bincount(S[X ^ d] ^ S[X], minlength=256)
    N = ext_table(C)
    W1, Wc = D / Q, N / Q

    print("\n=== Part A: theory numbers ===")
    Ws = W1[1:, 1:]
    Pz = np.eye(255) - 1 / 255
    sigma = np.linalg.svd(Pz @ Ws, compute_uv=False)[0]
    print(f"sigma_perp(S) = {sigma:.6f} (paper 0.210276), squared {sigma**2:.6f} (paper 0.044216)")
    print(f"det D*_S mod 5 = {det_mod_p(D[1:, 1:], 5)} (paper 1)")
    out["sigma_perp"] = sigma

    # Channel Q_{c,l} = W_c P_l W_1 for every coefficient value l present in L.
    dmin = {}
    for l in sorted(set(LM.flatten().tolist())):
        Qcl = Wc @ W1[MK[l, X]]
        dmin[l] = -math.log(bhattacharyya_rho(Qcl))
    # Paper indexes bytes as in RFC a_15..a_0 or reversed; test both orders.
    paper_t5 = {0: (11, 0xbd, 0.001020865), 1: (12, 0xe9, 0.001003868), 2: (4, 0x6b, 0.000995577),
                3: (9, 0x2f, 0.000976755), 4: (11, 0xe9, 0.001003868), 5: (5, 0xbc, 0.000964278),
                6: (13, 0x5a, 0.000981690), 7: (13, 0xe0, 0.000986993), 8: (0, 0x0a, 0.000996427),
                9: (9, 0x54, 0.000989626), 10: (8, 0x2a, 0.001016197), 11: (7, 0xb1, 0.001003688),
                12: (2, 0x9f, 0.000965553), 13: (7, 0x2a, 0.001016197), 14: (3, 0xe8, 0.001012424),
                15: (10, 0xbd, 0.001020865)}
    for name, idx in [("state order", lambda k: k), ("reversed (a_0 first)", lambda k: 15 - k)]:
        hits = sum(LM[idx(i), idx(j)] == l for j, (i, l, _) in paper_t5.items())
        print(f"Table 5 coefficients l_ij matching under {name}: {hits}/16")
    order = (lambda k: k) if sum(LM[i, j] == l for j, (i, l, _) in paper_t5.items()) == 16 else (lambda k: 15 - k)
    design = {}
    print(" j  i  l_ij  Dmin(ours)    Dmin(paper)")
    for j in range(16):
        jj = order(j)
        i_best = max(range(16), key=lambda i: dmin[int(LM[order(i), jj])])
        l = int(LM[order(i_best), jj])
        pi, pl, pd = paper_t5[j]
        design[j] = (i_best, l)
        mark = "" if (i_best, l) == (pi, pl) and abs(dmin[l] - pd) < 5e-9 else "  <-- differs"
        print(f"{j:2d} {i_best:2d}  {l:02x}  {dmin[l]:.9f}  {pd:.9f}{mark}")
    d1 = min(dmin[l] for _, l in design.values())
    N1 = math.ceil(math.log(16 * 255 / 0.005) / d1)
    rho2 = bhattacharyya_rho(Wc)
    N2 = math.ceil(math.log(16 * 255 / 0.005) / -math.log(rho2))
    print(f"worst D(1)min = {d1:.9f} (paper 0.000964278) -> N1 = {N1} (paper 14117)")
    print(f"rho_2 = {rho2:.16f} (paper 0.7074682424211026) -> N2 = {N2} (paper 40)")
    print(f"pairs = {16 * (N1 + N2)} (paper 226512), queries <= {17 * (N1 + N2)} (paper 240669)")
    out.update({"N1": N1, "N2": N2, "rho2": rho2, "D1min": d1})

    # ---------- Part B: run the attack ----------
    print("\n=== Part B: attack on the real two-round cipher ===")
    rng = np.random.default_rng(2026)

    def make_oracle(master_hex):
        rk = KZ.key_schedule(master_hex)
        K = [np.frombuffer(bytes.fromhex(k), dtype=np.uint8).astype(np.int64) for k in rk[:3]]
        def enc(P):
            return apply(LM, S[apply(LM, S[P ^ K[0]]) ^ K[1]]) ^ K[2]
        # cross-check the vectorised model with the reference code on a few blocks
        P = rng.integers(0, 256, (4, 16))
        for p, c in zip(P, enc(P)):
            ref = KZ.x_layer(rk[2], KZ.lsx_layer(rk[1], KZ.lsx_layer(rk[0], bytes(p.tolist()).hex())))
            assert ref == bytes(c.tolist()).hex(), "model != reference"
        return enc, K

    logQ = {}
    for j, (i, l) in design.items():
        Qcl = Wc @ W1[MK[l, X]]
        logQ[l] = np.log(np.maximum(Qcl, 1e-300))
    logW = np.log(np.maximum(Wc, 1e-300))
    one_c = 1 ^ C
    keys_alpha = MK[one_c, X]  # (1+c)k for every hypothesis k

    def attack(enc):
        queries = 0
        # Stage 1
        P = rng.integers(0, 256, (N1, 16))
        CP = enc(P)
        queries += N1
        k1 = np.zeros(16, dtype=np.int64)
        for j, (i, l) in design.items():
            jj, ii = order(j), order(i)
            a = int(rng.integers(0, 256))
            P2 = P.copy()
            P2[:, jj] = MK[C, P[:, jj]] ^ a
            G = apply(LINV, CP ^ enc(P2))[:, ii]
            queries += N1
            h = np.bincount(G, minlength=256)
            ll = logQ[l][a ^ keys_alpha] @ h
            k1[jj] = int(np.argmax(ll))
        # Stage 2
        V = rng.integers(0, 256, (N2, 16))
        PV = k1 ^ SINV[apply(LINV, V)]
        CV = enc(PV)
        queries += N2
        k2 = np.zeros(16, dtype=np.int64)
        for j in range(16):
            a = int(rng.integers(0, 256))
            V2 = V.copy()
            V2[:, j] = MK[C, V[:, j]] ^ a
            B = apply(LINV, CV ^ enc(k1 ^ SINV[apply(LINV, V2)]))[:, j]
            queries += N2
            ll = logW[a ^ keys_alpha] @ np.bincount(B, minlength=256)
            k2[j] = int(np.argmax(ll))
        return k1, k2, queries

    rfc = "8899aabbccddeeff0011223344556677fedcba98765432100123456789abcdef"
    n_rand = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    masters = [rfc] + [bytes(rng.integers(0, 256, 32).tolist()).hex() for _ in range(n_rand)]
    ok1 = ok2 = 0
    for m in masters:
        enc, K = make_oracle(m)
        k1, k2, qn = attack(enc)
        s1, s2 = (k1 == K[0]).all(), (k2 == K[1]).all()
        ok1 += s1
        ok2 += s1 and s2
        rec = (bytes(k1.tolist()) + bytes(k2.tolist())).hex()
        print(f"  key {m[:16]}..: K1 {'OK' if s1 else 'FAIL'}, K2 {'OK' if s2 else 'FAIL'}, "
              f"master recovered: {rec == m}, queries {qn}")
    print(f"full 256-bit key recovered for {ok2}/{len(masters)} keys")
    out.update({"keys_tested": len(masters), "K1_ok": int(ok1), "full_ok": int(ok2)})
    json.dump(out, open(f"{OUT}/stage5_two_round.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
