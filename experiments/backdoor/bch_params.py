"""B5: how much freedom did the designers have in choosing L?

Augot-Finiasz recursive MDS: R = companion matrix of g(X) = prod_{i=0}^{15} (X + beta^(l+i)),
beta of order 255 in a field GF(2^8). L = R^16 is then MDS (BCH bound, length 255 >= 32).

1. Rebuild g from (field 0x1C3, beta = alpha^64, l = 120) and check it reproduces
   the RFC coefficients exactly (no extra parameter).
2. Enumerate the whole family in field 0x1C3: all beta of order 255, all l.
   Palindromic g (roots closed under inversion) forces l = 120.
3. Among the palindromic candidates, see what singles out Kuznyechik's choice
   (count of coefficients equal to 1, i.e. free multiplications in hardware),
   and check that the best ones are one Frobenius orbit (c -> c^2).
4. Same count over all 16 primitive field polynomials.

Usage (repo root): python3 experiments/backdoor/bch_params.py
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verify"))
sys.path.insert(0, "src/gost")
import kuznechik as KZ
from gf import POLY_K, mul_clmul, log_tables
from layer_keys_audit import is_primitive


def poly_from_roots(roots, poly):
    g = [1]  # coefficients, g[k] = coeff of X^k
    for r in roots:
        ng = [0] * (len(g) + 1)
        for k, c in enumerate(g):
            ng[k + 1] ^= c                       # X * g
            ng[k] ^= mul_clmul(c, r, poly)        # r * g  (char 2: minus = plus)
        g = ng
    return g


def lfsr_coeffs(g):
    """Kuznyechik R: new byte = sum_i L_COEFF[i] * state[i], state shifts right.
    For a monic g of degree 16 with g[0] = 1 this corresponds to L_COEFF[i] = g[15 - i]
    for i < 15 and L_COEFF[15] = g[0] (checked against the RFC below)."""
    return [g[15 - i] for i in range(15)] + [g[0]]


def family(poly):
    exp, log = log_tables(poly)
    out = []
    for s in range(1, 255):
        if math.gcd(s, 255) != 1:
            continue
        for l in range(255):
            roots = [exp[(s * (l + i)) % 255] for i in range(16)]
            out.append((s, l, lfsr_coeffs(poly_from_roots(roots, poly))))
    return out


def main():
    res = {}
    exp, _ = log_tables(POLY_K)
    g = poly_from_roots([exp[(64 * (120 + i)) % 255] for i in range(16)], POLY_K)
    ours = lfsr_coeffs(g)
    match = ours == KZ.L_COEFF
    print(f"1. (0x1C3, beta=alpha^64, l=120) reproduces RFC L coefficients: {match}")
    print(f"   {ours}")
    res["reproduces_rfc"] = match

    fam = family(POLY_K)
    distinct = {tuple(c) for _, _, c in fam}
    pal = [(s, l, c) for s, l, c in fam if c[:15] == c[:15][::-1]]
    pal_distinct = {tuple(c) for _, _, c in pal}
    print(f"2. field 0x1C3: {len(fam)} (beta, l) pairs -> {len(distinct)} distinct R; "
          f"palindromic: {len(pal_distinct)} distinct, l values used: {sorted({l for _, l, _ in pal})}")
    ones = sorted(((sum(1 for v in c[:15] if v == 1), c) for c in pal_distinct), reverse=True)
    hist = {}
    for k, _ in ones:
        hist[k] = hist.get(k, 0) + 1
    ours_ones = sum(1 for v in KZ.L_COEFF[:15] if v == 1)
    print(f"3. palindromic candidates by number of coefficients equal to 1 (excluding the forced last one): {hist}")
    print(f"   Kuznyechik has {ours_ones}; candidates with at least that many: "
          f"{sum(1 for k, _ in ones if k >= ours_ones)}")
    for k, c in ones[:4]:
        print(f"     {k} ones: {c}{'   <- Kuznyechik' if list(c) == KZ.L_COEFF else ''}")
    res.update(family=len(distinct), palindromic=len(pal_distinct), ones_hist=hist, ours_ones=ours_ones)
    best = {c for k, c in ones if k == ours_ones}
    sq = lambda c: tuple(mul_clmul(v, v, POLY_K) for v in c)
    orbit = [tuple(KZ.L_COEFF)]
    while sq(orbit[-1]) != orbit[0]:
        orbit.append(sq(orbit[-1]))
    same = set(orbit) == best
    print(f"   Frobenius orbit of Kuznyechik's L has size {len(orbit)}; equals the set of best candidates: {same}")
    res["best_is_frobenius_orbit"] = same

    prims = [p for p in range(0x101, 0x200, 2) if is_primitive(p)]
    per_field = {}
    for p in prims:
        f = family(p)
        pd = {tuple(c) for _, l, c in f if c[:15] == c[:15][::-1]}
        per_field[hex(p)] = {"palindromic": len(pd),
                             "max_ones": max(sum(1 for v in c[:15] if v == 1) for c in pd),
                             "n_with_max": sum(1 for c in pd if sum(1 for v in c[:15] if v == 1) ==
                                               max(sum(1 for v in d[:15] if v == 1) for d in pd))}
    print("4. per field: palindromic count, max #ones, #candidates attaining it")
    for p, v in per_field.items():
        print(f"   {p}: {v}")
    res["per_field"] = per_field
    os.makedirs("results/backdoor", exist_ok=True)
    json.dump(res, open("results/backdoor/bch_params.json", "w"), indent=1)


if __name__ == "__main__":
    main()
