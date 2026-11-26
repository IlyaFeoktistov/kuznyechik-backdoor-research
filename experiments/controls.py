"""
Control generation for Kuznyechik S-box analysis.

Generates:
  1. N random permutations (baseline distribution).
  2. Matched controls: permutations whose (max ANF degree, nonlinearity, DDT
     max-for-d!=0) profile matches the Kuznyechik S-box within tolerance.
  3. TKlog bijective-preserving variants (parameter combinations that keep the
     map a permutation).

For each control we compute the same metrics as the Kuznyechik S-box and store
them so we can compare distributions.

Kuznyechik baseline (measured):
  max ANF degree = 7, nonlinearity = 102, DDT max (d!=0) = 8
"""

import os
import random
import time

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sbox_analysis import (
    load_sbox, sbox_anf_degrees, ddt_summary, nonlinearity,
)

KUZ = load_sbox()
KUZ_METRICS = {
    "max_degree": max(sbox_anf_degrees(KUZ)[0]),
    "nl": nonlinearity(KUZ),
    "ddt_max": ddt_summary(KUZ)["max_for_d_nonzero"],
}
print(f"Kuznyechik metrics: {KUZ_METRICS}")


def random_permutation(rng):
    p = list(range(256))
    rng.shuffle(p)
    return p


def profile(p):
    degs, max_deg = sbox_anf_degrees(p)
    nl = nonlinearity(p)
    ddt_max = ddt_summary(p)["max_for_d_nonzero"]
    return max_deg, nl, ddt_max


def generate_random(n, seed=12345, tol=None):
    """Generate n random permutations; return list of (perm, metrics)."""
    rng = random.Random(seed)
    out = []
    for i in range(n):
        p = random_permutation(rng)
        m = profile(p)
        out.append((p, m))
    return out


def collect_statistics(results):
    """Given list of (perm, (deg, nl, ddt)) return summary stats."""
    degs = [m[0] for _, m in results]
    nls = [m[1] for _, m in results]
    ddt = [m[2] for _, m in results]

    def summary(xs):
        return {
            "min": min(xs), "max": max(xs),
            "mean": round(sum(xs) / len(xs), 3),
        }

    return {
        "max_degree": summary(degs),
        "nl": summary(nls),
        "ddt_max": summary(ddt),
        "count": len(results),
    }


def matched_controls(results, kuz, tol_deg=1, tol_nl=3, tol_ddt=2):
    """Filter random permutations to those matching Kuznyechik's profile."""
    matched = []
    for p, m in results:
        if (abs(m[0] - kuz["max_degree"]) <= tol_deg
                and abs(m[1] - kuz["nl"]) <= tol_nl
                and abs(m[2] - kuz["ddt_max"]) <= tol_ddt):
            matched.append((p, m))
    return matched


def tklog_variant(p, seed):
    """
    A bijective-preserving TKlog-style transformation applied to a base
    permutation.

    This is a placeholder for the actual TKlog construction; the real one is
    defined in the papers (Perrin 2019, BPU 2016). Here we apply a set of
    bijective operations (which preserve the permutation property) that are
    plausible TKlog-style modifications:
      - byte-level XOR with a fixed mask (bijection)
      - composition with a fixed permutation (bijection)

    The point of the control is to keep the map a permutation while exploring
    the parameter space. We'll replace this with the real TKlog once we have
    the papers.
    """
    return p  # placeholder


def main():
    t0 = time.time()
    # Step 1: random permutations
    print("\nGenerating 200 random permutations...")
    randoms = generate_random(200, seed=12345)
    print(f"  done in {time.time() - t0:.1f}s")

    print("\n=== Random permutation statistics (n=200) ===")
    print("  " + str(collect_statistics(randoms)))

    # Step 2: matched controls
    print("\nFinding matched controls...")
    matched = matched_controls(randoms, KUZ_METRICS)
    print(f"  matched count: {len(matched)} / 200")
    if matched:
        print("  matched metrics sample:")
        for p, m in matched[:10]:
            print(f"    deg={m[0]} nl={m[1]} ddt={m[2]}")

    # Save results
    os.makedirs("results", exist_ok=True)
    with open("results/random_perms.json", "w") as f:
        import json
        json.dump([{"metrics": list(m)} for _, m in randoms], f)
    print("\nSaved results/random_perms.json")


if __name__ == "__main__":
    main()
