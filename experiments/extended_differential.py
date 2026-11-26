"""Verify the extended differential results of Stanica, Dutta, Mandal (arXiv:2507.02181v4).

Reproduces Table 2 of the paper:
    c   b   hidden count h   rank r   lower bound   actual uniformity
    02  12  120              1        60            64
    04  26  120              2        30            33
    10  24  120              4        8             11
    1d  30  120              4        8             12
    03  a2  32               1        16            21

Definitions (Stanica 2025, Theorem 6):
  - hidden transition (field T = GF(2^8)/(X^8+X^4+X^3+X^2+1), poly 0x11D):
        h_c(b) = #{x : S(c_T * x) ^ S(x) = b}
  - Kuznyechik extended differential (field K = GF(2^8)/(X^8+X^7+X^6+X+1), poly 0x1C3):
        delta_c(b) = max_a #{x : S(c_K * x ^ a) ^ S(x) = b}
  - disagreement rank r = rank_F2(M_K^c ^ M_T^c)  (Theorem 7: r = deg(c))
  - lower bound: delta_c(b) >= h_c(b) / 2^r
"""
import sys
sys.path.insert(0, "src/gost")
sys.path.insert(0, "experiments")
import kuznechik as K
from sbox_analysis import load_sbox

# hidden field T: poly 0x11D, reduction 0x1D
def gf_mul_T(a, b):
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

def gf_pow_T(a, e):
    result = 1
    while e:
        if e & 1:
            result = gf_mul_T(result, a)
        a = gf_mul_T(a, a)
        e >>= 1
    return result

def gf_inv_K(a):
    return K.gf_pow(a, 254) if hasattr(K, "gf_pow") else None

def rank_of_linear_map(values):
    """rank of the F2-linear span of the given 8-bit vectors"""
    basis = []
    for v in values:
        x = v
        for b in basis:
            if x and x.bit_length() - 1 == b.bit_length() - 1:
                x ^= b
        if x:
            basis.append(x)
    return len(basis)

def main():
    S = load_sbox()
    # K multiplication = kuznechik.gf_mul (0x1C3)
    mul_K = K.gf_mul

    print("=== Hidden transitions h_c(b) in field T (0x11D) ===")
    hidden_cases = [(0x02, 0x12, 120), (0x04, 0x26, 120), (0x10, 0x24, 120), (0x1d, 0x30, 120), (0x03, 0xa2, 32)]
    for c, b, expected_h in hidden_cases:
        counts = {}
        for x in range(256):
            d = S[gf_mul_T(c, x)] ^ S[x]
            counts[d] = counts.get(d, 0) + 1
        maxb = max(counts.values())
        hb = counts.get(b, 0)
        print(f"  c={c:#04x} b={b:#04x}: h(b)={hb} (expected {expected_h}), max_b h={maxb}")

    print("\n=== Disagreement ranks r = rank(M_K^c ^ M_T^c) ===")
    for c, expected_r in [(0x02, 1), (0x04, 2), (0x10, 4), (0x1d, 4), (0x03, 1)]:
        D = [mul_K(c, x) ^ gf_mul_T(c, x) for x in range(256)]
        r = rank_of_linear_map(D)
        # Theorem 7: rank = deg(c). compute deg of c.
        deg = c.bit_length() - 1
        print(f"  c={c:#04x}: rank={r} (expected {expected_r}), deg(c)={deg}")

    print("\n=== Kuznyechik extended differential uniformities delta_c(b) (field K, 0x1C3) ===")
    ed_cases = [(0x02, 0x12, 64), (0x04, 0x26, 33), (0x10, 0x24, 11), (0x1d, 0x30, 12), (0x03, 0xa2, 21)]
    for c, b, expected_d in ed_cases:
        # delta_c(b) = max_a #{x : S(c_K*x ^ a) ^ S(x) = b}
        best = 0
        best_a = -1
        for a in range(256):
            cnt = 0
            for x in range(256):
                if S[mul_K(c, x) ^ a] ^ S[x] == b:
                    cnt += 1
            if cnt > best:
                best = cnt
                best_a = a
        print(f"  c={c:#04x} b={b:#04x}: delta={best} (expected {expected_d}), at a={best_a:#04x}")

    print("\n=== Inverse multiplier classes (Proposition 3) ===")
    for c in [0x02, 0x04, 0x03]:
        cinv = K.gf_inv(c)
        b = {0x02: 0x12, 0x04: 0x26, 0x03: 0xa2}[c]
        for cc in (c, cinv):
            best = 0
            for a in range(256):
                cnt = 0
                for x in range(256):
                    if S[mul_K(cc, x) ^ a] ^ S[x] == b:
                        cnt += 1
                if cnt > best:
                    best = cnt
            print(f"  multiplier {cc:#04x} (b={b:#04x}): delta={best}")
        print(f"  inverse of {c:#04x} in K = {cinv:#04x}")

    print("\n=== Extended differential uniformity (max over a AND b) for all c ===")
    # For fixed c: for each a, build count distribution over b in 256*x, take max, then max over a.
    # 255 * 256 * 256 = ~16.7M ops.
    spec = []
    for c in range(1, 256):
        best = 0
        for a in range(256):
            counts = [0] * 256
            for x in range(256):
                counts[S[mul_K(c, x) ^ a] ^ S[x]] += 1
            mb = max(counts)
            if mb > best:
                best = mb
        spec.append((best, c))
    spec.sort(reverse=True)
    print("  top 15 (uniformity, c):")
    for d, c in spec[:15]:
        print(f"    c={c:#04x}: {d}")
    # baseline: random permutation would have uniformity ~256/256*... actually max of 256*256 binomial.
    print("  (for reference, a random permutation has extended diff uniformity ~ 3-4)")

    print("\n=== Extended differential verification complete ===")

if __name__ == "__main__":
    main()
