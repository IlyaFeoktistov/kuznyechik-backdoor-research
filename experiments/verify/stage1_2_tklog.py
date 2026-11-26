"""Stage 1: byte-exact TKlog reconstruction. Stage 2: its linear geometry.

Sage's GF(2^8) is backed by Givaro, whose `_log_repr()` stores elements as
logarithms with the convention zero -> 0 and one -> q-1 = 255 (NOT 0).
Every other nonzero element alpha^k (1 <= k <= 254) is stored as k.
Perrin's script feeds that value straight into i = l % 17, j = l // 17,
so x = 1 gives l = 255, i = 0, j = 15, kappa(1) = 0xEE.

Run from the repo root:  python3 experiments/verify/stage1_2_tklog.py
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf import POLY_T, mul_clmul, log_tables, f2_rank, f2_span

OUT = "results/verify"
S_ORIG = list(open("data/original/pi.bin", "rb").read())

S_PERM = [0, 12, 9, 8, 7, 4, 14, 6, 5, 10, 2, 11, 1, 3, 13]
LAMBDA = [0x12, 0x26, 0x24, 0x30]
CSTT = 0xFC

exp, log = log_tables(POLY_T)


def givaro_log_repr(x):
    """Emulates Sage/Givaro F.fetch_int(x)._log_repr() for x != 0."""
    return 255 if x == 1 else log[x]


def kappa(x):
    r = 0
    for j in range(4):
        if (x >> j) & 1:
            r ^= LAMBDA[j]
    return r ^ CSTT


def tklog(log_fn):
    pi = [kappa(0)]
    for x in range(1, 256):
        l = log_fn(x)
        i, j = l % 17, l // 17
        if i == 0:
            pi.append(kappa(16 - j))
        else:
            pi.append(kappa(16 - i) ^ exp[(17 * S_PERM[j]) % 255])
    return pi


def report(name, pi):
    mism = [x for x in range(256) if pi[x] != S_ORIG[x]]
    print(f"[{name}] mismatches: {len(mism)}; first: "
          + ", ".join(f"x={x:#04x}: got {pi[x]:#04x} want {S_ORIG[x]:#04x}" for x in mism[:5]))
    print(f"    pi(0)=0x{pi[0]:02x} (want 0xfc)  pi(1)=0x{pi[1]:02x} (want 0xee)  "
          f"permutation: {len(set(pi)) == 256}")
    return mism


def main():
    os.makedirs(OUT, exist_ok=True)
    print("=== Stage 1: TKlog reconstruction ===")
    report("ordinary dlog, log(1)=0", tklog(lambda x: log[x]))
    pi = tklog(givaro_log_repr)
    mism = report("Givaro _log_repr, log(1)=255", pi)
    if mism:
        sys.exit("STOP: TKlog does not reproduce the S-box")
    print("    -> byte-exact match with data/original/pi.bin")

    print("\n=== Stage 2: linear geometry ===")
    F16 = sorted(x for x in range(256) if exp_pow(x, 16) == x)
    print(f"GF(16) = {{x : x^16 = x}}: size {len(F16)}")
    f16set = set(F16)
    closed_add = all((a ^ b) in f16set for a in F16 for b in F16)
    closed_mul = all(mul_clmul(a, b, POLY_T) in f16set for a in F16 for b in F16)
    print(f"  closed under + : {closed_add}; closed under * : {closed_mul} -> subfield: {closed_add and closed_mul}")
    print(f"  GF(16)* = <alpha^17>: {set(exp[17 * k % 255] for k in range(15)) == f16set - {0}}")

    U = f2_span(LAMBDA)
    print(f"U = span(lambda) : dim {f2_rank(LAMBDA)}, |U| = {len(U)}")
    print(f"U ∩ GF(16) = {sorted(U & f16set)}")
    print(f"U + GF(16) = GF(256): {len({u ^ v for u in U for v in F16}) == 256}")

    # Every nonzero x = alpha^l, coset index i = l mod 17 (Givaro convention for l).
    # Perrin: coset i != 0 maps into the additive coset kappa(16-i) + GF(16);
    # coset 0 (= GF(16)*) together with 0 maps into kappa(U) i.e. CSTT + U.
    rows, ok = [], True
    for x in range(256):
        if x == 0:
            cls, layer_ok = "zero", pi[0] in {CSTT ^ u for u in U}
        else:
            l = givaro_log_repr(x)
            i = l % 17
            if i == 0:
                cls, layer_ok = 0, pi[x] in {CSTT ^ u for u in U}
            else:
                cls, layer_ok = i, (pi[x] ^ kappa(16 - i)) in f16set
        ok &= layer_ok
        rows.append({"x": x, "log_repr": None if x == 0 else givaro_log_repr(x),
                     "class": cls, "pi": pi[x], "ok": layer_ok})
    print(f"all 256 elements land in the predicted additive layer: {ok}")
    # The 16 nontrivial classes go to 16 distinct cosets of GF(16)
    cosets = {min(pi[r['x']] ^ f for f in F16) for r in rows if r["class"] not in ("zero", 0)}
    print(f"distinct additive cosets for classes 1..16: {len(cosets)} (want 16)")

    json.dump({"match": True, "pi": pi, "F16": F16, "U": sorted(U), "rows": rows},
              open(f"{OUT}/tklog_classes.json", "w"), indent=1)
    print(f"saved {OUT}/tklog_classes.json")


def exp_pow(x, e):
    r = 1
    for _ in range(e):
        r = mul_clmul(r, x, POLY_T)
    return r


if __name__ == "__main__":
    main()
