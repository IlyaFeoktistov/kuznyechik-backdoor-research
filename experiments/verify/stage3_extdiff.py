"""Stage 3+4: complete extended differential tables of the Kuznyechik S-box and
check against every number of Stanica-Dutta-Mandal (arXiv:2507.02181v4),
Tables 1-4 and Section 4.4.

N_c(a, b) = #{x : S(c*x ^ a) ^ S(x) = b}, multiplication in the Kuznyechik
field K (0x1C3). For each (c, a) one pass over x fills a 256-counter row.

Run from the repo root:  python3 experiments/verify/stage3_extdiff.py
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf import POLY_K, POLY_T, mul_table, f2_rank

sys.path.insert(0, "src/gost")
import kuznechik as KZ  # used ONLY as an independent cross-check of the field

OUT = "results/verify"
S = np.frombuffer(open("data/original/pi.bin", "rb").read(), dtype=np.uint8).astype(np.int64)
MK = np.array(mul_table(POLY_K), dtype=np.int64)
MT = np.array(mul_table(POLY_T), dtype=np.int64)
X = np.arange(256)
RESULTS = {}
FAIL = []


def check(label, got, want, tol=0.0):
    ok = abs(got - want) <= tol if isinstance(want, float) else got == want
    if not ok:
        FAIL.append(label)
    print(f"  {'OK ' if ok else 'BAD'} {label}: got {got}, paper {want}")
    RESULTS[label] = {"got": got, "paper": want, "ok": bool(ok)}


def ext_table(c, M=MK):
    """Full 256x256 table N_c[a, b]."""
    Sx = S[X]
    cx = M[c, X]
    N = np.zeros((256, 256), dtype=np.int64)
    for a in range(256):
        N[a] = np.bincount(S[cx ^ a] ^ Sx, minlength=256)
    return N


def naive(c, a, b):
    """Direct definition with the cipher's own gf_mul (independent code path)."""
    s = list(S)
    return sum(1 for x in range(256) if s[KZ.gf_mul(c, x) ^ a] ^ s[x] == b)


def ddt():
    D = np.zeros((256, 256), dtype=np.int64)
    for d in range(256):
        D[d] = np.bincount(S[X ^ d] ^ S[X], minlength=256)
    return D


def main():
    os.makedirs(OUT, exist_ok=True)
    print("Field cross-check: clmul(0x1C3) vs kuznechik.gf_mul on all 65536 pairs:",
          all(MK[a, b] == KZ.gf_mul(a, b) for a in range(256) for b in range(256)))

    tables = {c: ext_table(c) for c in range(1, 256)}
    D = tables[1]
    assert (D == ddt()).all()

    print("\n--- Sanity (Stage 4) ---")
    rng = np.random.default_rng(1)
    bad = 0
    for _ in range(300):
        c, a, b = (int(v) for v in rng.integers(1, 256, 3))
        bad += naive(c, a, b) != tables[c][a, b]
    print(f"  naive vs table on 300 random (c,a,b): {300 - bad}/300 agree")
    ds = all((N.sum(0) == 256).all() and (N.sum(1) == 256).all() for N in tables.values())
    print(f"  every N_c has all row and column sums = 256 (doubly stochastic): {ds}")
    z = all((tables[c][:, 0] == 1).all() for c in range(2, 256))
    print(f"  Lemma 25, N_c(a,0)=1 for all a, c!=1: {z}")
    distinct = all(len({r.tobytes() for r in tables[c]}) == 256 for c in range(2, 256))
    print(f"  rows of N_c pairwise distinct for every c!=1: {distinct}")

    print("\n--- Table 1: extended differential uniformity ---")
    delta = {1: int(D[1:].max())}
    delta.update({c: int(tables[c].max()) for c in range(2, 256)})
    check("delta_01 (ordinary DDT, a!=0)", delta[1], 8)
    for c, w in [(0x02, 64), (0xE1, 64), (0x04, 33), (0x91, 33), (0x03, 21), (0xBE, 21)]:
        check(f"delta_{c:02x}", delta[c], w)
    spectrum = sorted(((d, c) for c, d in delta.items() if c != 1), reverse=True)
    print("  top-12 of whole spectrum:", ", ".join(f"{c:02x}:{d}" for d, c in spectrum[:12]))
    print(f"  median delta_c over c!=1: {int(np.median([d for d, _ in spectrum]))}")

    print("\n--- Table 2: hidden-field certificates ---")
    for c, b, h, r, lb, dc in [(0x02, 0x12, 120, 1, 60, 64), (0x04, 0x26, 120, 2, 30, 33),
                                (0x10, 0x24, 120, 4, 8, 11), (0x1D, 0x30, 120, 4, 8, 12),
                                (0x03, 0xA2, 32, 1, 16, 21)]:
        hc = int(np.sum((S[MT[c, X]] ^ S[X]) == b))
        rank = f2_rank([int(MK[c, 1 << k] ^ MT[c, 1 << k]) for k in range(8)])
        check(f"c={c:02x} hidden h(b={b:02x})", hc, h)
        check(f"c={c:02x} rank(M_K ^ M_T)", rank, r)
        check(f"c={c:02x} bound ceil(h/2^r)", math.ceil(hc / 2 ** rank), lb)
        check(f"c={c:02x} delta_c", delta[c], dc)

    print("\n--- Section 4.4: how hidden solutions split over offsets ---")
    for c, b in [(0x02, 0x12), (0x04, 0x26), (0x03, 0xA2)]:
        hidden_x = [x for x in range(256) if S[MT[c, x]] ^ S[x] == b]
        parts = {}
        for x in hidden_x:
            a = int(MK[c, x] ^ MT[c, x])  # offset making cx^a = c_T x in K
            parts[a] = parts.get(a, 0) + 1
        col = {a: int(tables[c][a, b]) for a in parts}
        print(f"  c={c:02x} b={b:02x}: hidden parts {{offset: count}} = "
              + ", ".join(f"{a:02x}:{n}" for a, n in sorted(parts.items(), key=lambda t: -t[1]))
              + " | full table N_c(offset,b): " + ", ".join(f"{a:02x}:{v}" for a, v in col.items()))
    check("N_02(de,12)", int(tables[2][0xDE, 0x12]), 64)
    check("N_02(00,12)", int(tables[2][0x00, 0x12]), 56)
    n04 = sorted(int(tables[4][a, 0x26]) for a in (0x00, 0x7F, 0xDE, 0xA1))
    check("N_04(*,26) at 00,7f,de,a1 (sorted)", n04, [24, 32, 33, 33])

    print("\n--- Table 4 / Table 3: energy and mutual information ---")
    def energy_info(c):
        N = tables[c]
        E = int((N * N).sum())
        nz = N[N > 0].astype(float)
        I = float((nz * np.log2(nz)).sum() / 256 ** 2)
        return E, E / 2 ** 16 - 1, I
    for c, E0, eta0, I0 in [(0x02, 145324, 1.217468, 0.887664), (0x04, 138108, 1.107361, 0.858480),
                            (0x03, 130792, 0.995728, 0.816620)]:
        for cc in (c, int(KZ.gf_inv(c))):
            E, eta, I = energy_info(cc)
            check(f"E_{cc:02x}", E, E0)
            check(f"eta_{cc:02x}", round(eta, 6), eta0, 1e-6)
            check(f"I_{cc:02x}(K;B)", round(I, 6), I0, 1e-6)
    for c, I0 in [(0x10, 0.832516), (0x1D, 0.818254)]:
        check(f"I_{c:02x}(K;B) (Table 3)", round(energy_info(c)[2], 6), I0, 1e-6)
    # Theorem 11: E_c = sum_{alpha,beta} D(alpha,beta) D(c*alpha,beta)
    th11 = all(int((tables[c] ** 2).sum()) == int((D * D[MK[c, X]]).sum()) for c in range(1, 256))
    print(f"  Theorem 11 identity holds for all 255 multipliers: {th11}")

    print("\n--- Section 4.4: Pearson corr(B(c), delta_c) over 254 multipliers ---")
    Hmax = {lam: int(np.bincount(S[MT[lam, X]] ^ S[X], minlength=256)[1:].max() if lam == 1
                     else np.bincount(S[MT[lam, X]] ^ S[X], minlength=256).max())
            for lam in range(1, 256)}
    Bc, dc = [], []
    for c in range(2, 256):
        best = 0.0
        for lam in range(1, 256):
            r = f2_rank([int(MK[c, 1 << k] ^ MT[lam, 1 << k]) for k in range(8)])
            best = max(best, math.ceil(Hmax[lam] / 2 ** r))  # paper uses ceiling
        Bc.append(best)
        dc.append(delta[c])
    corr = float(np.corrcoef(Bc, dc)[0, 1])
    check("corrPearson(B(c), delta_c)", round(corr, 10), 0.9788451787, 1e-9)
    viol = sum(1 for b, d in zip(Bc, dc) if d < b)
    print(f"  Theorem 6 lower bound violated for {viol} multipliers (want 0)")

    np.save(f"{OUT}/ext_tables.npy", np.stack([tables[c] for c in range(1, 256)]).astype(np.uint16))
    json.dump({"delta": {f"{c:02x}": d for c, d in delta.items()}, "checks": RESULTS},
              open(f"{OUT}/stage3_checks.json", "w"), indent=1)
    print(f"\nsaved {OUT}/ext_tables.npy (index c-1) and stage3_checks.json")
    print(f"FAILED checks: {FAIL if FAIL else 'none'}")


if __name__ == "__main__":
    main()
