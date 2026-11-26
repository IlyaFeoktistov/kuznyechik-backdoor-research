"""Random MDS linear layers over GF(2^8)/0x1C3 for Stage F.

Construction: Cauchy matrix M[i,j] = 1/(x_i + y_j) with x, y two disjoint sets
of 16 distinct field elements (every square submatrix of a Cauchy matrix is
nonsingular, so M is MDS), followed by random row/column permutations and
random nonzero row/column scalings (both preserve MDS).

Limitation: this samples the Cauchy-equivalent subclass of MDS matrices, NOT
the uniform distribution over all 16x16 MDS matrices. It is generated without
reference to the Kuznyechik matrix.

Independent checks (check_mds):
  * inverse exists (Gaussian elimination over GF(2^8));
  * every 1x1, 2x2 and 3x3 minor nonzero (exhaustive) and 2000 random larger
    square minors nonzero;
  * branch number: all weight-1 and weight-2 inputs give output weight >= 16
    and >= 15, i.e. input + output weight >= 17 for those (exhaustive).
"""
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verify"))
from gf import POLY_K, mul_table, log_tables

MK = np.array(mul_table(POLY_K), dtype=np.int64)
EXP, LOG = log_tables(POLY_K)
INV = np.zeros(256, dtype=np.int64)
for _x in range(1, 256):
    INV[_x] = EXP[(255 - LOG[_x]) % 255]


def cauchy_mds(rng):
    pts = rng.permutation(256)[:32]
    x, y = pts[:16], pts[16:]
    M = INV[np.bitwise_xor.outer(x, y)]
    M = M[rng.permutation(16)][:, rng.permutation(16)]
    r, c = rng.integers(1, 256, 16), rng.integers(1, 256, 16)
    M = MK[MK[M, r[:, None]], c[None, :]]
    return M


def gf_inverse(M):
    n = len(M)
    A = np.concatenate([M.copy(), np.eye(n, dtype=np.int64)], axis=1)
    for col in range(n):
        piv = next((r for r in range(col, n) if A[r, col]), None)
        if piv is None:
            return None
        A[[col, piv]] = A[[piv, col]]
        A[col] = MK[INV[A[col, col]], A[col]]
        for r in range(n):
            if r != col and A[r, col]:
                A[r] ^= MK[A[r, col], A[col]]
    return A[:, n:]


def gf_det_nonzero(M):
    A = M.copy()
    n = len(A)
    for col in range(n):
        piv = next((r for r in range(col, n) if A[r, col]), None)
        if piv is None:
            return False
        A[[col, piv]] = A[[piv, col]]
        inv = INV[A[col, col]]
        for r in range(col + 1, n):
            if A[r, col]:
                A[r] ^= MK[MK[A[r, col], inv], A[col]]
    return True


def apply(M, V):
    out = np.zeros_like(V)
    for i in range(16):
        for t in range(16):
            out[..., i] ^= MK[M[i, t], V[..., t]]
    return out


def check_mds(M, rng, n_random=2000):
    Minv = gf_inverse(M)
    ok = {"invertible": Minv is not None and (apply(M, apply(Minv, np.eye(16, dtype=np.int64))) == np.eye(16)).all()}
    ok["minors_1"] = bool((M != 0).all())
    ok["minors_2_3"] = all(gf_det_nonzero(M[np.ix_(r, c)])
                           for k in (2, 3)
                           for r in itertools.combinations(range(16), k)
                           for c in itertools.combinations(range(16), k))
    good = True
    for _ in range(n_random):
        k = int(rng.integers(4, 17))
        r, c = np.sort(rng.permutation(16)[:k]), np.sort(rng.permutation(16)[:k])
        good &= gf_det_nonzero(M[np.ix_(r, c)])
    ok["minors_random_4_16"] = bool(good)
    # branch number on weight-1 and weight-2 inputs, exhaustive
    w1 = min(int(np.count_nonzero(MK[M[:, t], b])) for t in range(16) for b in range(1, 256))
    vals = np.arange(1, 256)
    w2 = 16
    for t1, t2 in itertools.combinations(range(16), 2):
        out = MK[M[:, t1][:, None, None], vals[None, :, None]] ^ MK[M[:, t2][:, None, None], vals[None, None, :]]
        w2 = min(w2, int(np.count_nonzero(out, axis=0).min()))
    ok["branch_w1"] = 1 + w1
    ok["branch_w2"] = 2 + w2
    ok["mds"] = ok["invertible"] and ok["minors_1"] and ok["minors_2_3"] and ok["minors_random_4_16"] \
        and ok["branch_w1"] >= 17 and ok["branch_w2"] >= 17
    return ok, Minv
