"""Priority 1: large-scale controls for every attack-relevant metric of the
extended-differential approach, maximised over ALL 254 non-classical
multipliers (max-statistic handles the multiple comparisons over c).

Per S-box metrics (direction "more dangerous"):
  dmax  = max_c delta_c                        (higher)
  emax  = max_c E_c = sum N_c^2                (higher; = max Pearson info about whitening byte)
  n2min = min_c N2(c), stage-2 samples (Sec. 6.4) (lower)
  n1min = min_c N1(c), stage-1 samples (Sec. 6.3, full 2-round attack) (lower; optional, expensive)

Control families:
  R: uniformly random permutations
  A: A o pi o B with random invertible affine A, B over F2^8. Same DDT/LAT
     spectra, degree, NL as pi, but the hidden field is no longer aligned
     with the Kuznyechik field.

Usage (repo root):
  OMP_NUM_THREADS=1 python3 experiments/verify/control_large.py --n 1000 --seed 1
  OMP_NUM_THREADS=1 python3 experiments/verify/control_large.py --n 64 --seed 2 --n1
"""
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf import POLY_K, mul_table, f2_rank

MK = np.array(mul_table(POLY_K), dtype=np.int64)
X = np.arange(256)
KUZ = np.frombuffer(open("data/original/pi.bin", "rb").read(), dtype=np.uint8).astype(np.int64)
LOGB = math.log(16 * 255 / 0.005)
_LM = None


def lm():
    global _LM
    if _LM is None:
        from stage5_two_round import LM
        _LM = LM
    return _LM


def rho(P):
    R = np.sqrt(P)
    B = R @ R.T
    np.fill_diagonal(B, -1)
    return B.max()


def metrics(S, with_n1=False):
    S = np.asarray(S, dtype=np.int64)
    D = np.stack([np.bincount(S[X ^ d] ^ S[X], minlength=256) for d in range(256)])
    res = {"delta1": int(D[1:].max())}
    dmax = emax = -1
    n2min = n1min = 10 ** 9
    argc = {}
    if with_n1:
        LM = lm()
        coeffs = sorted(set(LM.flatten().tolist()))
        W1 = D / 256
    for c in range(2, 256):
        N = np.stack([np.bincount(S[MK[c, X] ^ a] ^ S[X], minlength=256) for a in range(256)])
        d, e = int(N.max()), int((N * N).sum())
        n2 = math.ceil(LOGB / -math.log(rho(N / 256)))
        if d > dmax: dmax, argc["dmax"] = d, c
        if e > emax: emax, argc["emax"] = e, c
        if n2 < n2min: n2min, argc["n2min"] = n2, c
        if with_n1:
            Wc = N / 256
            dm = {l: -math.log(rho(Wc @ W1[MK[l, X]])) for l in coeffs}
            worst = min(max(dm[int(LM[i, j])] for i in range(16)) for j in range(16))
            n1 = math.ceil(LOGB / worst)
            if n1 < n1min: n1min, argc["n1min"] = n1, c
    res.update(dmax=dmax, emax=emax, n2min=n2min, argc=argc)
    if with_n1:
        res["n1min"] = n1min
    return res


def rand_affine(rng):
    while True:
        cols = [int(v) for v in rng.integers(1, 256, 8)]
        if f2_rank(cols) == 8:
            break
    lin = np.zeros(256, dtype=np.int64)
    for x in range(256):
        for k in range(8):
            if (x >> k) & 1:
                lin[x] ^= cols[k]
    return lin ^ int(rng.integers(0, 256))


def make_sbox(args):
    family, seed, with_n1 = args
    rng = np.random.default_rng(seed)
    if family == "R":
        S = rng.permutation(256)
    else:
        A, B = rand_affine(rng), rand_affine(rng)
        S = A[KUZ[B]]
    r = metrics(S, with_n1)
    r.update(family=family, seed=seed)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n1", action="store_true")
    ap.add_argument("--procs", type=int, default=16)
    a = ap.parse_args()
    tag = f"n{a.n}_seed{a.seed}{'_n1' if a.n1 else ''}"
    t0 = time.time()
    kuz = metrics(KUZ, a.n1)
    print("Kuznyechik:", kuz, flush=True)
    jobs = [(f, a.seed * 10 ** 6 + (0 if f == "R" else 5 * 10 ** 5) + i, a.n1)
            for f in ("R", "A") for i in range(a.n)]
    with Pool(a.procs) as pool:
        rows = pool.map(make_sbox, jobs, chunksize=1)
    keys = [("dmax", 1), ("emax", 1), ("n2min", -1)] + ([("n1min", -1)] if a.n1 else [])
    summary = {}
    for fam in ("R", "A"):
        fr = [r for r in rows if r["family"] == fam]
        print(f"\n== family {fam} (n={len(fr)}) ==")
        for k, sign in keys:
            v = np.array([r[k] for r in fr])
            # one-sided empirical p: fraction of controls at least as "dangerous" as pi
            worse = int(np.sum(sign * v >= sign * kuz[k]))
            p = (1 + worse) / (1 + len(v))
            summary[f"{fam}.{k}"] = {"kuz": kuz[k], "min": int(v.min()), "median": float(np.median(v)),
                                     "max": int(v.max()), "p_one_sided": p}
            print(f"  {k:6s}: pi={kuz[k]:>8}  controls min/median/max = {v.min()}/{np.median(v)}/{v.max()}"
                  f"  p(control at least as dangerous) = {p:.4f}")
    m = len(keys) * 2
    print(f"\nBonferroni threshold over {m} tests at 0.05: p < {0.05 / m:.4f}")
    print(f"time {time.time() - t0:.0f}s")
    os.makedirs("results/verify", exist_ok=True)
    json.dump({"params": vars(a), "kuznyechik": kuz, "summary": summary, "rows": rows},
              open(f"results/verify/control_large_{tag}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
