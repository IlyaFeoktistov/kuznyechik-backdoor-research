"""Priority 2 (bounded): product-model characteristic scores of Sec. 7 of
arXiv:2507.02181v4 for pi and for controls, maximised over ALL multipliers.

Family (no pre-whitening core, as in the paper):
  2R: 16 -> 1      log2 P = sum_i log2 m*_c(D'_i) - 128 + log2 D(beta, gamma) - 8
  3R: 16 -> 1 -> 16  + sum_i log2 d*(L(gamma e_k)_i) - 128
where D' = L^{-1}(beta e_k), m*_c(b) = max_{a != 0} N_c(a, b), d*(x) = max_y D(x, y).
c = 1 is the classical case. The attacker's best is max over all c.

Usage (repo root):
  OMP_NUM_THREADS=1 python3 experiments/verify/char_search.py --n 200 --seed 3
"""
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf import POLY_K, mul_table
from control_large import KUZ, rand_affine
from stage5_two_round import LM, LINV

MK = np.array(mul_table(POLY_K), dtype=np.int64)
X = np.arange(256)
NZ = np.arange(1, 256)
# DL[k, beta, i] = (L^{-1}(beta e_k))_i ;  DF[k, gamma, i] = (L(gamma e_k))_i
DL = np.stack([MK[LINV[:, k][None, :], NZ[:, None]] for k in range(16)])
DF = np.stack([MK[LM[:, k][None, :], NZ[:, None]] for k in range(16)])


def scores(S):
    S = np.asarray(S, dtype=np.int64)
    D = np.stack([np.bincount(S[X ^ d] ^ S[X], minlength=256) for d in range(256)])
    with np.errstate(divide="ignore"):
        logD = np.log2(D[1:, 1:].astype(float))            # [beta-1, gamma-1]
    dstar = np.log2(np.maximum(D.max(1), 1).astype(float))  # d*(x); x=0 unused
    third = dstar[DF].sum(2) - 128                          # [k, gamma-1]
    best2_round2 = logD.max(1) - 8                          # best gamma for each beta
    out = {}
    for c in range(1, 256):
        N = np.stack([np.bincount(S[MK[c, X] ^ a] ^ S[X], minlength=256) for a in range(256)])
        mstar = np.log2(np.maximum(N[1:].max(0), 1).astype(float))  # a != 0
        mstar[N[1:].max(0) == 0] = -np.inf
        first = mstar[DL].sum(2) - 128                       # [k, beta-1]
        s2 = (first + best2_round2[None, :]).max()
        s3 = max((first[k][:, None] + logD - 8 + third[k][None, :]).max() for k in range(16))
        out[c] = (float(s2), float(s3))
    return out


def summarize(sc):
    c2 = max(sc, key=lambda c: sc[c][0])
    c3 = max(sc, key=lambda c: sc[c][1])
    return {"classical2": sc[1][0], "classical3": sc[1][1],
            "best2": sc[c2][0], "best2_c": c2, "best3": sc[c3][1], "best3_c": c3}


def job(args):
    fam, seed = args
    rng = np.random.default_rng(seed)
    S = rng.permutation(256) if fam == "R" else rand_affine(rng)[KUZ[rand_affine(rng)]]
    r = summarize(scores(S))
    r.update(family=fam, seed=seed)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--procs", type=int, default=16)
    a = ap.parse_args()
    t0 = time.time()
    sc = scores(KUZ)
    print("Paper check (Tables 6, 7):")
    for c, w2, w3 in [(0x01, -89.1504, -174.3), (0x02, -83.9801, -169.7), (0xE1, -83.9801, -169.7),
                      (0x04, -89.9814, None), (0x03, -91.1955, -177.4)]:
        print(f"  c={c:02x}: 2R {sc[c][0]:.4f} (paper {w2})   3R {sc[c][1]:.2f} (paper {w3})")
    kuz = summarize(sc)
    print("pi:", kuz)
    jobs = [(f, a.seed * 10 ** 6 + (0 if f == "R" else 5 * 10 ** 5) + i) for f in ("R", "A") for i in range(a.n)]
    with Pool(a.procs) as pool:
        rows = pool.map(job, jobs, chunksize=1)
    summary = {}
    for fam in ("R", "A"):
        fr = [r for r in rows if r["family"] == fam]
        print(f"\n== family {fam} (n={len(fr)}) ==")
        for k in ("classical2", "best2", "classical3", "best3"):
            v = np.array([r[k] for r in fr])
            p = (1 + int(np.sum(v >= kuz[k]))) / (1 + len(v))
            summary[f"{fam}.{k}"] = {"kuz": kuz[k], "min": v.min(), "median": float(np.median(v)),
                                     "max": v.max(), "p_control_at_least_as_good_for_attacker": p}
            print(f"  {k:10s}: pi={kuz[k]:8.2f}  controls min/median/max = "
                  f"{v.min():.2f}/{np.median(v):.2f}/{v.max():.2f}  p = {p:.4f}")
    print(f"time {time.time() - t0:.0f}s")
    json.dump({"params": vars(a), "kuznyechik": kuz, "kuz_all_c": {f"{c:02x}": v for c, v in sc.items()},
               "summary": summary, "rows": rows},
              open(f"results/verify/char_search_n{a.n}_seed{a.seed}.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
