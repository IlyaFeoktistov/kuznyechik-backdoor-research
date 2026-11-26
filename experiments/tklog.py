"""Independent verification of the TKlog structure of the Kuznyechik/Streebog S-box.

Port of the SAGE script published by L. Perrin (FAQ https://who.paris.inria.fr/Leo.Perrin/pi.html,
section 2.1.1) and the paper "Partitions in the S-Box of Streebog and Kuznyechik" (IACR ePrint 2019/092).

TKlog parameters (from Perrin 2019):
  hidden field F = GF(2^8)/(X^8 + X^4 + X^3 + X^2 + 1)   -> reduction 0x11D
  alpha = generator of F^*
  s = [0, 12, 9, 8, 7, 4, 14, 6, 5, 10, 2, 11, 1, 3, 13]   (perm of 0..14)
  lambda_vectors = [0x12, 0x26, 0x24, 0x30]
  cstt = 0xFC

Core property (Perrin 2019, sec 2.2):
  pi maps the multiplicative cosets of F_{16}^* in F_{28}^* to additive cosets of F_{16}.
  For a != 1: pi(a * b) = f(a) XOR g(b), with g a permutation of F_{16}.
"""
import sys
sys.path.insert(0, "src/gost")
sys.path.insert(0, "experiments")
import kuznechik as K
from sbox_analysis import load_sbox

# ---------------- Hidden field GF(2^8)/(X^8+X^4+X^3+X^2+1) = 0x11D ----------------
POLY = 0x11D
RED = 0x1D  # reduction constant = lower 8 bits of X^8+X^4+X^3+X^2+1 = 0x11D & 0xFF

def gf_mul(a, b):
    """GF(2^8) multiplication with polynomial X^8+X^4+X^3+X^2+1 (0x11D)."""
    result = 0
    for _ in range(8):
        if b & 1:
            result ^= a
        b >>= 1
        carry = a & 0x80
        a = (a << 1) & 0xFF
        if carry:
            a ^= 0x1D
    return result

def gf_pow(a, e):
    result = 1
    while e:
        if e & 1:
            result = gf_mul(result, a)
        a = gf_mul(a, a)
        e >>= 1
    return result

def gf_inv(a):
    return gf_pow(a, 254)

def main():
    # verify alpha = X is a generator of F^* (order 255)
    alpha = 0x02  # X = 0b00000010
    # build log table
    log = [-1] * 256
    el = 1
    for l in range(255):
        log[el] = l
        el = gf_mul(el, alpha)
    assert el == 1, "alpha is not a generator (order != 255)"
    assert all(log[x] >= 0 for x in range(1, 256)), "log table incomplete"
    print("alpha=0x02 is a generator of F^* (order 255): OK")

    # ----- TKlog generation (port of SAGE script) -----
    s = [0, 12, 9, 8, 7, 4, 14, 6, 5, 10, 2, 11, 1, 3, 13]
    lambda_vectors = [0x12, 0x26, 0x24, 0x30]
    cstt = 0xFC

    def kappa(x):
        result = 0
        for j in range(0, 4):
            if (x >> j) & 1 == 1:
                result ^= lambda_vectors[j]
        return result ^ cstt

    pi = [kappa(0)]
    for x in range(1, 256):
        l = log[x]
        i, j = l % 17, l // 17
        if i == 0:
            y = kappa(16 - j)
        else:
            gf_elmt = gf_pow(gf_pow(alpha, 17), s[j])
            y = kappa(16 - i) ^ gf_elmt
        pi.append(y)

    pi = pi[:256]
    assert len(pi) == 256
    # is pi a permutation?
    assert len(set(pi)) == 256, "TKlog generation is not a permutation!"
    print("TKlog generation produces a permutation: OK")

    # compare with the actual S-box
    actual = load_sbox()
    mismatches = [x for x in range(256) if pi[x] != actual[x]]
    print(f"TKlog vs actual pi.bin: {len(mismatches)} mismatches")
    if mismatches:
        print(f"  first mismatches: {mismatches[:10]}")
        for x in mismatches[:10]:
            print(f"    x={x:#04x}: tklog={pi[x]:#04x} actual={actual[x]:#04x}")
    else:
        print("  *** TKlog EXACTLY reproduces the Kuznyechik S-box! ***")

    # ----- verify the coset-to-coset property -----
    # F_16 = {x : x^16 == x}
    F16 = sorted(x for x in range(256) if gf_pow(x, 16) == x)
    print(f"\nF_16 size = {len(F16)} (expect 16)")
    F16set = set(F16)
    F16_star = [x for x in F16 if x != 0]
    # multiplicative cosets of F_16^* in F_28^*: coset_i = {alpha^i * F_16^*} for i in 0..16
    # compute coset index of each nonzero element
    def coset_index(x):
        # x = alpha^i * (alpha^17)^j ; i = log(x) % 17
        return log[x] % 17
    # group outputs by input coset
    from collections import defaultdict
    in_cosets = defaultdict(list)
    for x in range(1, 256):
        in_cosets[coset_index(x)].append(x)
    print(f"Number of nontrivial multiplicative cosets (i=1..16): {len(in_cosets)-1} (expect 16), incl i=0: {len(in_cosets)} (expect 17)")

    # For each input coset i, check that all outputs lie in ONE additive coset of F_16
    print("\nCoset-to-coset verification (each input coset maps into a single additive coset of F_16):")
    ok_all = True
    for i in sorted(in_cosets):
        outs = [pi[x] for x in in_cosets[i]]
        # check all outs are in the same additive coset: outs[0] + F_16 contains all
        base = outs[0]
        same = all(((o ^ base) in F16set) for o in outs)
        if not same:
            ok_all = False
            print(f"  coset i={i}: NOT in single additive coset")
    print(f"  All {len(in_cosets)} input cosets map into a single additive coset of F_16: {ok_all}")

    # distinctness of output cosets
    coset_rep = {}
    for i in sorted(in_cosets):
        outs = [pi[x] for x in in_cosets[i]]
        coset_rep[i] = outs[0] ^ F16[0]  # representative (coset = base + F16)
    distinct = len(set(coset_rep.values())) == len(in_cosets)
    print(f"  All {len(in_cosets)} output cosets are distinct: {distinct}")

    # ----- verify pi(a*b) = f(a) XOR g(b) for a != 1 -----
    # For fixed i (coset of a), pi(alpha^i * (alpha^17)^j) = kappa(16-i) XOR (alpha^17)^s[j]
    # So f(i) = kappa(16-i), and g(j) = (alpha^17)^s[j] (a permutation of F_16^*)
    print("\nDecomposition pi(a*b) = f(a) XOR g(b) for a != 1:")
    # check g is a permutation of F_16^*
    g_vals = [gf_pow(gf_pow(alpha, 17), s[j]) for j in range(15)]
    print(f"  g(j) values are in F_16^*: {all(v in F16_star for v in g_vals)}")
    print(f"  g is a permutation (distinct): {len(set(g_vals)) == 15}")
    # verify the formula for a few samples
    ok = True
    for i in range(1, 17):
        for j in range(15):
            a = gf_pow(alpha, i)
            b = gf_pow(gf_pow(alpha, 17), j)
            x = gf_mul(a, b)
            f_a = kappa(16 - i)
            g_b = gf_pow(gf_pow(alpha, 17), s[j])
            if pi[x] != (f_a ^ g_b):
                ok = False
    print(f"  pi(a*b) = f(a) XOR g(b) holds for all a in coset 1..16, b in F_16^*: {ok}")

    # ----- field mismatch with Kuznyechik linear layer -----
    # Kuznyechik uses poly 0x1C3, hidden field uses 0x11D.
    # Show that the multiplication maps differ (disagreement rank) - see Stanica 2025.
    print("\nField mismatch (hidden 0x11D vs Kuznyechik 0x1C3):")
    # disagreement map D = M_K XOR M_T for a fixed multiplier z
    # M_K(z): x -> z*x in K (0x1C3); M_T(z): x -> z*x in T (0x11D)
    def mul_K(a, b):
        return K.gf_mul(a, b)  # 0x1C3
    z = 0x04  # one of the exceptional multipliers
    D = [mul_K(z, x) ^ gf_mul(z, x) for x in range(256)]
    # rank of D as an F2-linear map (D(0)=0): dimension of span of {D(x)}
    basis = []
    for d in D:
        v = d
        for b in basis:
            if v & (1 << (v.bit_length() - 1)) and (b & (1 << (b.bit_length() - 1))) == (v & (1 << (v.bit_length() - 1))):
                v ^= b
        if v != 0:
            basis.append(v)
    rank = len(basis)
    print(f"  disagreement rank D = M_K(z) XOR M_T(z) for z={z:#04x}: {rank}")

    print("\n=== TKlog verification complete ===")

if __name__ == "__main__":
    main()
