"""Sanity checks for the affine controls and an independent (slow, pure
Python, no numpy) recomputation of the 2-round score from char_search.py.

1. G = A o pi o B (A, B random affine) keeps the TKlog structure, only in a
   conjugated field: with B = M x ^ k, the relation
       G(c' x ^ a) ^ G(x) = A_lin(b)   for c' = M^-1 o (c_T *) o M
   has exactly as many solutions as pi has for the hidden multiplier c_T.
   So controls are "same structure, not aligned with the Kuznyechik field".
2. Pure-Python 2-round score for pi and two controls vs char_search.scores.

Run from the repo root:  python3 experiments/verify/check_controls.py
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf import POLY_K, POLY_T, mul_clmul
from control_large import KUZ, rand_affine
from char_search import scores
sys.path.insert(0, "src/gost")
import kuznechik as KZ

pi = [int(v) for v in KUZ]


def hidden_count(S, mult):
    """max_b #{x : S(mult(x)) ^ S(x) = b}, mult an arbitrary permutation of bytes."""
    cnt = {}
    for x in range(256):
        d = S[mult[x]] ^ S[x]
        cnt[d] = cnt.get(d, 0) + 1
    return max(cnt.values())


def main():
    rng = np.random.default_rng(11)
    cT = [mul_clmul(2, x, POLY_T) for x in range(256)]
    print("pi, hidden multiplier 02 in field T: max count", hidden_count(pi, cT))
    for t in range(3):
        A, B = rand_affine(rng), rand_affine(rng)
        G = [int(A[pi[int(B[x])]]) for x in range(256)]
        Binv = [0] * 256
        for x in range(256):
            Binv[int(B[x])] = x
        conj = [Binv[cT[int(B[x])]] for x in range(256)]  # B^-1 o c_T o B
        k_cnt = max(hidden_count(G, [mul_clmul(c, x, POLY_K) for x in range(256)]) for c in range(2, 256))
        print(f"control {t}: structure in conjugated field: {hidden_count(G, conj)} "
              f"(same as pi); best over Kuznyechik-field multipliers c*x: {k_cnt}")

    print("\nIndependent pure-Python 2-round score (classical c=1 and c=02):")
    LI = [[0] * 16 for _ in range(16)]
    for k in range(16):
        e = bytearray(16)
        e[k] = 1
        col = bytes.fromhex(KZ.l_inv_layer(e.hex()))
        for i in range(16):
            LI[i][k] = col[i]

    def score2(S, c):
        D = [[0] * 256 for _ in range(256)]
        for d in range(256):
            for x in range(256):
                D[d][S[x ^ d] ^ S[x]] += 1
        mstar = [0] * 256  # max over a != 0 of N_c(a, b)
        for a in range(1, 256):
            row = [0] * 256
            for x in range(256):
                row[S[mul_clmul(c, x, POLY_K) ^ a] ^ S[x]] += 1
            mstar = [max(m, r) for m, r in zip(mstar, row)]
        best = -math.inf
        for k in range(16):
            for beta in range(1, 256):
                s = sum(math.log2(mstar[mul_clmul(LI[i][k], beta, POLY_K)]) for i in range(16)) - 128
                s += math.log2(max(D[beta][1:])) - 8
                best = max(best, s)
        return best

    rng = np.random.default_rng(12)
    boxes = [("pi", pi), ("random", [int(v) for v in rng.permutation(256)]),
             ("affine", [int(v) for v in rand_affine(rng)[KUZ[rand_affine(rng)]]])]
    for name, S in boxes:
        fast = scores(S)
        for c in (1, 2):
            print(f"  {name:7s} c={c:02x}: slow {score2(S, c):.4f}  fast {fast[c][0]:.4f}")


if __name__ == "__main__":
    main()
