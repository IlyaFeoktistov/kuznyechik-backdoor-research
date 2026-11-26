"""Stage F: interaction of the S-box pi with the linear layer L.

Groups (same metrics, same attacker for all):
  A   pi            + Kuznyechik L
  B   pi            + random Cauchy MDS layer (field 0x1C3)
  Bb  pi            + the other 63 palindromic BCH layers of field 0x1C3
  C   A o pi o B    + Kuznyechik L       (random affine A, B)
  D   A o pi o B    + random Cauchy MDS layer (independent)
  R   random permutation + Kuznyechik L (baseline without hidden structure; the
      T multipliers are then just 254 extra linear maps for the attacker)

Attackers (multiplier sets, best taken over the set):
  K   multiplication by c in the Kuznyechik field, c = 2..255 (paper)
  KT  K plus multiplication in the hidden field T (for affine copies: T
      conjugated by B's linear part, i.e. the same structure seen from G's input)
The attacker's best is always max(classical c = 1, set), reported as K* and KT*.
Offset variants: nz (a != 0, paper family) and all (a = 0 allowed).
Rounds: 2, 3, 4 (16->1, 16->1->16, 16->1->16->1).

R3/R5 metrics (delta_c, E_c, I(K;B)) do not involve L at all, so they are
identical across A/B and across C/D by definition and are not recomputed.

Usage (repo root):
  OMP_NUM_THREADS=1 python3 experiments/interaction/run_stage_f.py --n 100 --seed 1
"""
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "verify"))
sys.path.insert(0, os.path.join(HERE, "..", "backdoor"))
from char_general import MK, MT, X, tails, best_over, k_multipliers, t_multipliers
from mds import cauchy_mds, check_mds, gf_inverse
from stage5_two_round import LM as L_ORIG, LINV as LINV_ORIG
from control_large import rand_affine

PI = np.frombuffer(open("data/original/pi.bin", "rb").read(), dtype=np.uint8).astype(np.int64)
KEYS = [f"{v}{r}" for v in ("nz", "all") for r in (2, 3, 4)]


def bch_layers():
    from bch_params import family
    from gf import POLY_K
    import kuznechik as KZ
    out = []
    for coeffs in sorted({tuple(c) for _, l, c in family(POLY_K) if list(c[:15]) == list(c[:15])[::-1]}):
        if list(coeffs) == KZ.L_COEFF:
            continue
        R = np.zeros((16, 16), dtype=np.int64)
        R[0, :] = coeffs
        for i in range(1, 16):
            R[i, i - 1] = 1
        Lm = np.eye(16, dtype=np.int64)
        for _ in range(16):  # Lm = R * Lm over GF(2^8)
            Lm = np.array([[np.bitwise_xor.reduce(MK[R[i], Lm[:, j]]) for j in range(16)] for i in range(16)])
        out.append((list(coeffs), Lm))
    return out


def evaluate(S, Lm, Linv, Blin):
    tl = tails(S, Lm, Linv)
    cls = best_over(S, {"id": X}, tl)
    k = best_over(S, k_multipliers(), tl)
    t = best_over(S, t_multipliers(Blin), tl)
    kt = {key: max(k[key], t[key]) for key in KEYS}
    return {"classical": {key: cls[key][0] for key in KEYS},
            "K*": {key: max(cls[key][0], k[key][0]) for key in KEYS},
            "KT*": {key: max(cls[key][0], kt[key][0]) for key in KEYS},
            "K": {key: k[key][0] for key in KEYS}, "K_arg": {key: k[key][1] for key in KEYS},
            "KT": {key: kt[key][0] for key in KEYS}, "KT_arg": {key: kt[key][1] for key in KEYS}}


def job(args):
    group, seed, extra = args
    rng = np.random.default_rng(seed)
    S, Blin = PI, None
    if group == "R":
        S = rng.permutation(256)
    if group in ("C", "D"):
        A, B = rand_affine(rng), rand_affine(rng)
        S = A[PI[B]]
        Blin = B ^ B[0]
        # the conjugated hidden-field structure must be fully present
        Mc = t_multipliers(Blin)["T02"]
        assert max(np.bincount(S[Mc ^ a] ^ S[X], minlength=256).max() for a in range(256)) == 120
    if group in ("B", "D"):
        Lm = cauchy_mds(rng)
        Linv = gf_inverse(Lm)
    elif group == "Bb":
        Lm = np.array(extra)
        Linv = gf_inverse(Lm)
    else:
        Lm, Linv = L_ORIG, LINV_ORIG
    r = evaluate(S, Lm, Linv, Blin)
    r.update(group=group, seed=seed)
    if group in ("B", "D"):
        r["L"] = Lm.tolist()
    return r


def describe(vals, ref):
    v = np.array(vals)
    q = np.quantile(v, [0.05, 0.25, 0.5, 0.75, 0.95])
    return {"n": len(v), "min": v.min(), "q05": q[0], "q25": q[1], "median": q[2], "q75": q[3], "q95": q[4],
            "max": v.max(), "ref": ref, "p_ge_ref": (1 + int(np.sum(v >= ref - 1e-9))) / (1 + len(v)),
            "ref_minus_median": ref - q[2]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--procs", type=int, default=16)
    ap.add_argument("--check", type=int, default=5, help="Cauchy matrices to run the full MDS check on")
    a = ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(a.seed + 999)
    checks = []
    for i in range(a.check):
        ok, _ = check_mds(cauchy_mds(np.random.default_rng(a.seed * 10 ** 6 + 100000 + i)), rng, 500)
        checks.append(ok)
    print("MDS checks on generated Cauchy layers:", all(c["mds"] for c in checks), flush=True)

    base = a.seed * 10 ** 6
    jobs = [("A", base, None)]
    jobs += [("B", base + 100000 + i, None) for i in range(a.n)]
    jobs += [("Bb", base + 200000 + i, Lm.tolist()) for i, (_, Lm) in enumerate(bch_layers())]
    jobs += [("C", base + 300000 + i, None) for i in range(a.n)]
    jobs += [("D", base + 400000 + i, None) for i in range(a.n)]
    jobs += [("R", base + 500000 + i, None) for i in range(a.n)]
    with Pool(a.procs) as pool:
        rows = pool.map(job, jobs, chunksize=1)
    A = rows[0]
    print("A (pi + Kuznyechik L):")
    for att in ("classical", "K*", "KT*"):
        print(f"  {att:9s}", {k: round(v, 2) for k, v in A[att].items()})
    summary = {}
    for g in ("B", "Bb", "C", "D", "R"):
        gr = [r for r in rows if r["group"] == g]
        print(f"\n== group {g} (n={len(gr)}) : median [q05, q95] max ; A's value ; p(control >= A) ==")
        for att in ("classical", "K*", "KT*"):
            for key in KEYS:
                d = describe([r[att][key] for r in gr], A[att][key])
                summary[f"{g}.{att}.{key}"] = d
                print(f"  {att:9s} {key}: {d['median']:8.2f} [{d['q05']:8.2f}, {d['q95']:8.2f}] {d['max']:8.2f} ;"
                      f" A {A[att][key]:8.2f} ; p {d['p_ge_ref']:.3f}")
        # gain of the attacker over the classical attack, per instance
        for att in ("K*", "KT*"):
            for key in KEYS:
                gains = [r[att][key] - r["classical"][key] for r in gr]
                summary[f"{g}.gain_{att}.{key}"] = describe(gains, A[att][key] - A["classical"][key])
    print(f"\ntime {time.time() - t0:.0f}s")
    os.makedirs("results/interaction", exist_ok=True)
    json.dump({"params": vars(a), "mds_checks": checks, "A": A, "summary": summary, "rows": rows},
              open(f"results/interaction/stage_f_n{a.n}_seed{a.seed}.json", "w"), default=float)


if __name__ == "__main__":
    main()
