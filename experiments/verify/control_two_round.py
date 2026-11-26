"""Control for Section 6: Stage-1 sample count N1 of the same two-round attack
when the Kuznyechik S-box is replaced by random permutations (same L layer).

If random S-boxes need about as many samples, the attack is generic for
two rounds and is not powered by the hidden TKlog structure.

Run from the repo root:  python3 experiments/verify/control_two_round.py [n_random]
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stage5_two_round import LM, MK, X, bhattacharyya_rho

KUZ = np.frombuffer(open("data/original/pi.bin", "rb").read(), dtype=np.uint8).astype(np.int64)
COEFFS = sorted(set(LM.flatten().tolist()))
MULTIPLIERS = [0x02, 0x03, 0x04, 0x4F, 0x91, 0xAA, 0xE1, 0x55]


def n1(S, c):
    D = np.stack([np.bincount(S[X ^ d] ^ S[X], minlength=256) for d in range(256)]) / 256
    N = np.stack([np.bincount(S[MK[c, X] ^ a] ^ S[X], minlength=256) for a in range(256)]) / 256
    dmin = {l: -math.log(bhattacharyya_rho(N @ D[MK[l, X]])) for l in COEFFS}
    worst = min(max(dmin[int(LM[i, j])] for i in range(16)) for j in range(16))
    return math.ceil(math.log(16 * 255 / 0.005) / worst)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    rng = np.random.default_rng(7)
    res = {"kuznyechik": {f"{c:02x}": n1(KUZ, c) for c in MULTIPLIERS}, "random": []}
    print("Kuznyechik N1 per c:", res["kuznyechik"])
    for k in range(n):
        S = rng.permutation(256)
        row = {f"{c:02x}": n1(S, c) for c in MULTIPLIERS}
        res["random"].append(row)
        print(f"random #{k}: min N1 {min(row.values())}, c=04: {row['04']}")
    best = [min(r.values()) for r in res["random"]]
    print(f"\nKuznyechik best N1 over these c: {min(res['kuznyechik'].values())}")
    print(f"random S-boxes best N1: min {min(best)}, median {int(np.median(best))}, max {max(best)}")
    json.dump(res, open("results/verify/control_two_round.json", "w"), indent=1)


if __name__ == "__main__":
    main()
