"""Stages B and D: linear layer, constants, fields, key schedule.

B1  Round constants: C_i == L(Vec128(i)) for i = 1..32 and the RFC test values.
    (No free parameters: the constants are a fixed function of L.)
B2  Position of both field polynomials among the 16 primitive degree-8 polynomials.
B3  Fixed/periodic structure of L and R over F2: dim ker(M^k + I), k = 1..32,
    compared with random invertible 16x16 matrices over GF(2^8)_K.
D1  Master key from any consecutive round-key pair (K_{2i+1}, K_{2i+2}) by
    inverting the key-schedule Feistel, on 20 random keys (seed 31).

Usage (repo root): python3 experiments/backdoor/layer_keys_audit.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verify"))
sys.path.insert(0, "src/gost")
import kuznechik as KZ
from gf import POLY_K, mul_table
from stage5_two_round import LM

MK = np.array(mul_table(POLY_K), dtype=np.int64)
OUT = {}


def b1():
    rfc = {1: "6ea276726c487ab85d27bd10dd849401", 2: "dc87ece4d890f4b3ba4eb92079cbeb02",
           8: "f6593616e6055689adfba18027aa2a08"}
    ok = all(KZ.l_layer_r16(i.to_bytes(16, "big").hex()) == KZ.c_i(i) for i in range(1, 33))
    rfc_ok = all(KZ.c_i(i) == v for i, v in rfc.items())
    print(f"B1 C_i = L(Vec128(i)) for all 32: {ok}; matches RFC test values C_1, C_2, C_8: {rfc_ok}")
    OUT["B1"] = {"constants_from_L": ok, "rfc_match": rfc_ok}


def is_primitive(p):
    # order of X modulo p must be 255
    x, k = 2, 1
    while x != 1:
        x <<= 1
        if x & 0x100:
            x ^= p
        k += 1
        if k > 255:
            return False
    return k == 255


def b2():
    prims = [p for p in range(0x101, 0x200, 2) if is_primitive(p)]
    iT, iK = prims.index(0x11D), prims.index(0x1C3)
    print(f"B2 primitive degree-8 polynomials: {len(prims)}; "
          f"0x11D is #{iT + 1} (lexicographically), 0x1C3 is #{iK + 1}")
    print("   list:", " ".join(f"{p:#x}" for p in prims))
    OUT["B2"] = {"primitive": prims, "rank_11D": iT + 1, "rank_1C3": iK + 1}


def f2_matrix(M16):
    """128x128 F2 matrix (rows as Python ints) of a 16x16 matrix over GF(2^8)_K."""
    cols = []
    for t in range(16):
        for b in range(8):
            out = 0
            for i in range(16):
                out |= int(MK[M16[i, t], 1 << b]) << (8 * i)
            cols.append(out)
    return cols  # column images


def apply_cols(cols, v):
    r = 0
    k = 0
    while v:
        if v & 1:
            r ^= cols[k]
        v >>= 1
        k += 1
    return r


def kernel_dim(cols):
    # rank of the 128 column vectors
    piv = {}
    r = 0
    for v in cols:
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                piv[h] = v
                r += 1
                break
    return 128 - r


def periodic_dims(M16, kmax=32):
    cols = f2_matrix(M16)
    dims = []
    cur = [1 << j for j in range(128)]  # images of basis under M^k
    for k in range(1, kmax + 1):
        cur = [apply_cols(cols, v) for v in cur]
        dims.append(kernel_dim([c ^ (1 << j) for j, c in enumerate(cur)]))
    return dims


def b3():
    R = np.zeros((16, 16), dtype=np.int64)
    R[0, :] = KZ.L_COEFF
    for i in range(1, 16):
        R[i, i - 1] = 1
    assert (np.linalg.matrix_power(np.eye(1), 1) == 1).all()
    dL, dR = periodic_dims(LM), periodic_dims(R)
    print("B3 dim ker(L^k + I), k=1..32:", dL)
    print("   dim ker(R^k + I), k=1..32:", dR)
    rng = np.random.default_rng(41)
    ctrl = []
    for _ in range(20):
        M = rng.integers(0, 256, (16, 16))
        ctrl.append(periodic_dims(M, 16))
    print("   random 16x16 over GF(2^8) (n=20), max over controls per k=1..16:",
          [max(c[k] for c in ctrl) for k in range(16)])
    OUT["B3"] = {"L": dL, "R": dR, "controls_k1_16": ctrl}


def d1():
    rng = np.random.default_rng(31)
    ok = 0
    for _ in range(20):
        master = bytes(rng.integers(0, 256, 32).tolist()).hex()
        rk = KZ.key_schedule(master)
        good = True
        for i in range(1, 5):  # pair (K_{2i+1}, K_{2i+2}) at list index (2i, 2i+1)
            a1, a0 = rk[2 * i], rk[2 * i + 1]
            for step in range(i * 8, 0, -1):  # undo F[C_step]
                c = KZ.c_i(step)
                a1, a0 = a0, KZ.x_layer(a1, KZ.lsx_layer(c, a0))
            good &= (a1 + a0 == master)
        ok += good
    print(f"D1 master key recovered from every consecutive pair (K3,K4)...(K9,K10): {ok}/20 keys")
    OUT["D1"] = {"keys": 20, "recovered": ok}


def main():
    b1()
    b2()
    b3()
    d1()
    os.makedirs("results/backdoor", exist_ok=True)
    json.dump(OUT, open("results/backdoor/layer_keys_audit.json", "w"), indent=1)


if __name__ == "__main__":
    main()
