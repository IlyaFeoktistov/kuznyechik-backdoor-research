"""Stage A: S-box audit for design-level trapdoors (pi and pi^-1) against controls.

Tests (all exhaustive, no heuristics):
  part   Partition trapdoors / subspace trails. A partition invariant under XOR
         with every key must be the coset partition of a subspace V. S maps
         V-cosets into W-cosets iff D_v S(x) ^ D_v S(0) in W for all v in V, x.
         We compute, for every v != 0, r(v) = dim span{D_v S(x) ^ D_v S(0)}.
         min_v r(v) = 8  =>  no nontrivial (V, W) with dim W < 8 exists.
  lstr   Linear structures: pairs (u != 0, v != 0) with u . D_v S constant.
  q2,q3  Number of linearly independent F2 relations of degree <= 2 / <= 3
         in (x, S(x)) (algebraic-attack surface; AES S-box has 39 of degree 2).
  inv2   Nonlinear invariants g of degree <= 2 with g(S(x)) = g(x) + c.
  polyT, polyK  Number of nonzero coefficients of S as a univariate polynomial
         over the hidden field T (0x11D) and over the Kuznyechik field K (0x1C3).

Positive controls: a planted partition-trapdoor S-box and the AES S-box.
Usage (repo root): python3 experiments/backdoor/sbox_audit.py --n 1000 --seed 5
"""
import argparse
import json
import os
import sys
import time
from itertools import combinations
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verify"))
from gf import POLY_K, POLY_T, mul_clmul, log_tables, f2_rank
from control_large import rand_affine

PI = [int(b) for b in open("data/original/pi.bin", "rb").read()]


# ---------------------------------------------------------------- helpers
def rank_rows(rows):
    """F2 rank of a list of Python-int bit rows."""
    pivots = {}
    r = 0
    for v in rows:
        while v:
            h = v.bit_length() - 1
            if h in pivots:
                v ^= pivots[h]
            else:
                pivots[h] = v
                r += 1
                break
    return r


def monomials(nvars, deg):
    out = []
    for d in range(deg + 1):
        out.extend(combinations(range(nvars), d))
    return out


MON16 = {d: monomials(16, d) for d in (2, 3)}
MON8_2 = [m for m in monomials(8, 2) if m]  # deg 1..2, constant handled separately


def mon_bits(z, mons):
    bits = [(z >> i) & 1 for i in range(16)]
    v = 0
    for k, m in enumerate(mons):
        if all(bits[i] for i in m):
            v |= 1 << k
    return v


# ---------------------------------------------------------------- tests
def partition_test(S):
    r = []
    for v in range(1, 256):
        base = S[0] ^ S[v]
        r.append(f2_rank([S[x] ^ S[x ^ v] ^ base for x in range(256)]))
    return min(r)


def linear_structures(S):
    cnt = 0
    for v in range(1, 256):
        d = [S[x] ^ S[x ^ v] for x in range(256)]
        for u in range(1, 256):
            first = bin(u & d[0]).count("1") & 1
            if all((bin(u & y).count("1") & 1) == first for y in d):
                cnt += 1
    return cnt


def relations(S, deg):
    mons = MON16[deg]
    rows = []
    # columns = monomials, rows = points; rank over points equals column rank
    cols = [0] * len(mons)
    for x in range(256):
        b = mon_bits(x | (S[x] << 8), mons)
        rows.append(b)
    return len(mons) - rank_rows(rows)


def invariants2(S):
    """dim of {(g, c) : g deg<=2 without constant, g(S(x)) ^ g(x) = c}, minus trivial g=0."""
    nm = len(MON8_2)
    # unknown vector: coefficients of g (nm bits) + c (1 bit). One equation per x.
    eqs = []
    for x in range(256):
        a = mon_bits(x, MON8_2) ^ mon_bits(S[x], MON8_2)
        eqs.append(a | (1 << nm))
    # kernel dim = unknowns - rank
    return (nm + 1) - rank_rows(eqs)


def poly_weight(S, poly):
    exp, log = log_tables(poly)
    f = np.array(S)
    nz = np.array([x for x in range(1, 256) if S[x] != 0])
    lf, lx = np.array([log[S[x]] for x in nz]), np.array([log[x] for x in nz])
    cnt = int(S[0] != 0)  # a_0 = f(0)
    expa = np.array(exp)
    for k in range(1, 255):
        acc = np.bitwise_xor.reduce(expa[(lf - k * lx) % 255])
        cnt += int(acc != 0)
    a_last = np.bitwise_xor.reduce(f)  # a_255 = sum_x f(x)  (char 2)
    return cnt + int(a_last != 0)


def audit(S):
    Sinv = [0] * 256
    for x in range(256):
        Sinv[S[x]] = x
    return {
        "part": partition_test(S), "part_inv": partition_test(Sinv),
        "lstr": linear_structures(S), "q2": relations(S, 2), "q3": relations(S, 3),
        "inv2": invariants2(S), "inv2_inv": invariants2(Sinv),
        "polyT": poly_weight(S, POLY_T), "polyK": poly_weight(S, POLY_K),
    }


def aes_sbox():
    s = []
    for x in range(256):
        inv = 0 if x == 0 else next(y for y in range(256) if mul_clmul(x, y, 0x11B) == 1)
        r = inv
        for i in range(1, 5):
            r ^= ((inv << i) | (inv >> (8 - i))) & 0xFF
        s.append(r ^ 0x63)
    return s


def planted_trapdoor(seed):
    """S(h, l) = (P(h), Q_h(l)): maps cosets of the low-nibble space onto cosets of itself."""
    rng = np.random.default_rng(seed)
    P = rng.permutation(16)
    Q = [rng.permutation(16) for _ in range(16)]
    return [int(P[x >> 4] << 4 | Q[x >> 4][x & 15]) for x in range(256)]


def job(args):
    fam, seed = args
    rng = np.random.default_rng(seed)
    if fam == "R":
        S = [int(v) for v in rng.permutation(256)]
    else:
        A, B = rand_affine(rng), rand_affine(rng)
        S = [int(A[PI[int(B[x])]]) for x in range(256)]
    r = audit(S)
    r.update(family=fam, seed=seed)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--procs", type=int, default=16)
    a = ap.parse_args()
    t0 = time.time()
    # sanity check of the univariate interpolation on a random point set
    pi_res = audit(PI)
    print("pi       :", pi_res)
    print("AES      :", audit(aes_sbox()), " (positive control: q2 should be 39)")
    print("planted  :", audit(planted_trapdoor(1)), " (positive control: part < 8)")
    jobs = [(f, a.seed * 10 ** 6 + (0 if f == "R" else 5 * 10 ** 5) + i) for f in ("R", "A") for i in range(a.n)]
    with Pool(a.procs) as pool:
        rows = pool.map(job, jobs, chunksize=4)
    summary = {}
    for fam in ("R", "A"):
        fr = [r for r in rows if r["family"] == fam]
        print(f"\n== family {fam} (n={len(fr)}) ==")
        for k in pi_res:
            v = np.array([r[k] for r in fr])
            lo, hi = (1 + int(np.sum(v <= pi_res[k]))) / (1 + len(v)), (1 + int(np.sum(v >= pi_res[k]))) / (1 + len(v))
            summary[f"{fam}.{k}"] = {"pi": pi_res[k], "min": int(v.min()), "median": float(np.median(v)),
                                     "max": int(v.max()), "p_le": lo, "p_ge": hi}
            print(f"  {k:9s}: pi={pi_res[k]:>4}  controls min/median/max = {v.min()}/{np.median(v)}/{v.max()}"
                  f"  p(<=pi)={lo:.4f} p(>=pi)={hi:.4f}")
    print(f"time {time.time() - t0:.0f}s")
    os.makedirs("results/backdoor", exist_ok=True)
    json.dump({"params": vars(a), "pi": pi_res, "summary": summary, "rows": rows},
              open(f"results/backdoor/sbox_audit_n{a.n}_seed{a.seed}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
