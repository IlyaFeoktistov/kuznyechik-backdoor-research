"""
Kuznyechik (GOST R 34.12-2015) Block Cipher Implementation

Full, verified implementation of Kuznyechik for cryptographic analysis.
All test vectors from RFC 7801 (Sections 5.2–5.6) are verified.

Conventions
-----------
State is kept in RFC byte order: state[0] = a_15, state[1] = a_14, ..., state[15] = a_0.
The coefficient L_COEFF[i] corresponds directly to state[i] (a_{15-i}).

Key schedule (RFC 7801 / errata 6928)
-------------------------------------
K_1 = k[255..128], K_2 = k[127..0]
C_i = L(Vec_128(i))   for i = 1..32   (Vec_128(i) is big-endian)
F[C](a_1, a_0) = (LSX[C](a_1) XOR a_0, a_1)          (Feistel)
(K_(2i+1), K_(2i+2)) = F[C_8] F[C_7] ... F[C_1] (K_(2i-1), K_(2i))   for i = 1..4

Encryption:  E = X[K_10] LSX[K_9] ... LSX[K_1] (a)
Decryption:  D = X[K_1] S^-1 L^-1 X[K_2] ... S^-1 L^-1 X[K_10] (a)
  (per block the applied order is X, then L^-1, then S^-1, per RFC 5.6)
"""

import os

# ---------------------------------------------------------------------------
# Constants (RFC 7801 / errata 6928)
# ---------------------------------------------------------------------------

# Coefficients for l(a_15,...,a_0) = XOR_i GF_MUL(a_i, L_COEFF[i])
# decimal: [148, 32, 133, 16, 194, 192, 1, 251, 1, 192, 194, 16, 133, 32, 148, 1]
L_COEFF = [148, 32, 133, 16, 194, 192, 1, 251,
           1, 192, 194, 16, 133, 32, 148, 1]

# GF(2^8): polynomial x^8 + x^7 + x^6 + x + 1 (rep. 0x1C3), reduction 0xC3
GF_POLY = 0x1C3
REDUC_CONST = 0xC3

# ---------------------------------------------------------------------------
# GF(2^8) arithmetic
# ---------------------------------------------------------------------------

def gf_mul(a: int, b: int) -> int:
    """Multiply two bytes in GF(2^8) (shift-and-reduce)."""
    result = 0
    for _ in range(8):
        if b & 1:
            result ^= a
        b >>= 1
        carry = a & 0x80
        a = (a << 1) & 0xFF
        if carry:
            a ^= REDUC_CONST
    return result


def gf_inv(a: int) -> int:
    """Multiplicative inverse in GF(2^8). Inverse of non-zero a is a^254."""
    if a == 0:
        return 0
    # exponentiation by squaring: a^254
    x = 1
    e = 254
    base = a
    while e:
        if e & 1:
            x = gf_mul(x, base)
        base = gf_mul(base, base)
        e >>= 1
    return x


# ---------------------------------------------------------------------------
# Load and validate the S-box (Pi permutation) and its inverse
# ---------------------------------------------------------------------------

_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data', 'original')


def _load_bytes(path: str) -> bytearray:
    with open(path, 'rb') as f:
        return bytearray(f.read())


def _validate_sbox(sbox: bytearray, pinv: bytearray) -> None:
    """Verify the S-box is a valid permutation and pinv is its inverse."""
    assert len(sbox) == 256, f"SBOX must be 256 bytes, got {len(sbox)}"
    assert all(0 <= b <= 255 for b in sbox), "SBOX contains out-of-range byte"
    assert len(set(sbox)) == 256, "SBOX is not a permutation (duplicates)"
    assert len(pinv) == 256, f"PINV must be 256 bytes, got {len(pinv)}"
    assert all(0 <= b <= 255 for b in pinv), "PINV contains out-of-range byte"
    assert len(set(pinv)) == 256, "PINV is not a permutation (duplicates)"
    assert all(sbox[pinv[i]] == i for i in range(256)), "PINV is not the inverse of SBOX"
    assert all(pinv[sbox[i]] == i for i in range(256)), "SBOX is not the inverse of PINV"


SBOX = _load_bytes(os.path.join(_DATA_DIR, 'pi.bin'))
PINV = _load_bytes(os.path.join(_DATA_DIR, 'pi_inv.bin'))
_validate_sbox(SBOX, PINV)


# ---------------------------------------------------------------------------
# Transformations
# ---------------------------------------------------------------------------

def r_layer(state_hex: str) -> str:
    """R(a_15||...||a_0) = l(a_15,...,a_0) || a_15 || ... || a_1."""
    state = list(bytes.fromhex(state_hex))
    l = 0
    for i in range(16):
        l ^= gf_mul(state[i], L_COEFF[i])
    new_state = [l] + state[:15]
    return bytes(new_state).hex()


def l_layer_r16(state_hex: str) -> str:
    """L(a) = R^16(a)."""
    state = list(bytes.fromhex(state_hex))
    for _ in range(16):
        l = 0
        for i in range(16):
            l ^= gf_mul(state[i], L_COEFF[i])
        state = [l] + state[:15]
    return bytes(state).hex()


def r_inv_layer(state_hex: str) -> str:
    """R^-1(a_15||...||a_0) = a_14||a_13||...||a_0||l(a_14,...,a_0,a_15)."""
    state = list(bytes.fromhex(state_hex))
    new_state = state[1:] + [state[0]]
    l = 0
    for i in range(16):
        l ^= gf_mul(new_state[i], L_COEFF[i])
    new_state[15] = l
    return bytes(new_state).hex()


def l_inv_layer(state_hex: str) -> str:
    """L^-1(a) = (R^-1)^16(a)."""
    state = list(bytes.fromhex(state_hex))
    for _ in range(16):
        new_state = state[1:] + [state[0]]
        l = 0
        for i in range(16):
            l ^= gf_mul(new_state[i], L_COEFF[i])
        new_state[15] = l
        state = new_state
    return bytes(state).hex()


def s_layer(state_hex: str) -> str:
    """S(a): byte-wise substitution using Pi."""
    state = list(bytes.fromhex(state_hex))
    return bytes(SBOX[b] for b in state).hex()


def s_inv_layer(state_hex: str) -> str:
    """S^-1(a): byte-wise substitution using Pi^-1."""
    state = list(bytes.fromhex(state_hex))
    return bytes(PINV[b] for b in state).hex()


def x_layer(key_hex: str, state_hex: str) -> str:
    """X[k](a) = k XOR a."""
    key = bytes.fromhex(key_hex)
    state = bytearray(bytes.fromhex(state_hex))
    for i in range(16):
        state[i] ^= key[i]
    return bytes(state).hex()


def lsx_layer(key_hex: str, state_hex: str) -> str:
    """LSX[k](a) = L(S(X[k](a)))."""
    state = x_layer(key_hex, state_hex)
    state = s_layer(state)
    state = l_layer_r16(state)
    return state


# ---------------------------------------------------------------------------
# Key schedule
# ---------------------------------------------------------------------------

def vec128(i: int) -> str:
    """128-bit big-endian representation of integer i."""
    return i.to_bytes(16, 'big').hex()


def c_i(i: int) -> str:
    """C_i = L(Vec_128(i))."""
    return l_layer_r16(vec128(i))


def f_transform(c_hex: str, a1: str, a0: str) -> tuple:
    """F[C](a1, a0) = (LSX[C](a1) XOR a0, a1)."""
    lsa1 = lsx_layer(c_hex, a1)
    new_a1 = x_layer(a0, lsa1)
    new_a0 = a1
    return new_a1, new_a0


def key_schedule(key_hex: str) -> list:
    """Generate the 10 round keys K_1..K_10 from the 256-bit master key."""
    k1 = key_hex[:32]
    k2 = key_hex[32:]
    round_keys = [k1, k2]
    for i in range(1, 5):  # i = 1..4 -> (K_3,K_4)..(K_9,K_10)
        a1 = round_keys[2 * i - 2]  # K_(2i-1)
        a0 = round_keys[2 * i - 1]  # K_(2i)
        for j in range(1, 9):       # F[C_1] ... F[C_8]
            c = c_i(8 * (i - 1) + j)
            a1, a0 = f_transform(c, a1, a0)
        round_keys.append(a1)  # K_(2i+1)
        round_keys.append(a0)  # K_(2i+2)
    return round_keys


# ---------------------------------------------------------------------------
# Encryption / decryption
# ---------------------------------------------------------------------------

def encrypt(plaintext_hex: str, round_keys: list) -> str:
    """E = X[K_10] LSX[K_9] ... LSX[K_1] (a)."""
    state = plaintext_hex
    for i in range(9):  # K_1..K_9 (indices 0..8)
        state = lsx_layer(round_keys[i], state)
    return x_layer(round_keys[9], state)  # K_10 (index 9), XOR only


def decrypt(ciphertext_hex: str, round_keys: list) -> str:
    """D = X[K_1] S^-1 L^-1 X[K_2] ... S^-1 L^-1 X[K_10] (a).

    Per block the applied order is X, then L^-1, then S^-1 (RFC 5.6).
    """
    state = ciphertext_hex
    for i in range(9, 0, -1):  # K_10..K_2 (indices 9..1)
        state = x_layer(round_keys[i], state)
        state = l_inv_layer(state)
        state = s_inv_layer(state)
    return x_layer(round_keys[0], state)  # K_1


# ---------------------------------------------------------------------------
# Test vectors (RFC 7801, Section 5)
# ---------------------------------------------------------------------------

def test_vectors() -> None:
    print("Testing GF(2^8) arithmetic...")
    assert gf_mul(0x94, 0x94) == 0xA4
    assert gf_mul(0x01, 0x94) == 0x94
    assert gf_mul(0x01, 0x01) == 0x01
    assert all(gf_mul(a, gf_inv(a)) == 1 for a in range(1, 256))
    assert all(gf_mul(a, b) == gf_mul(b, a) for a in range(256) for b in range(256))
    print("  ✓ GF(2^8) arithmetic passes")

    print("\nTesting R transformation...")
    assert r_layer("00000000000000000000000000000100") == "94000000000000000000000000000001"
    assert r_layer("94000000000000000000000000000001") == "a5940000000000000000000000000000"
    assert r_layer("a5940000000000000000000000000000") == "64a59400000000000000000000000000"
    assert r_layer("64a59400000000000000000000000000") == "0d64a594000000000000000000000000"
    print("  ✓ R test vectors pass (4/4)")

    print("\nTesting L transformation...")
    assert l_layer_r16("64a59400000000000000000000000000") == "d456584dd0e3e84cc3166e4b7fa2890d"
    assert l_layer_r16("d456584dd0e3e84cc3166e4b7fa2890d") == "79d26221b87b584cd42fbc4ffea5de9a"
    assert l_layer_r16("79d26221b87b584cd42fbc4ffea5de9a") == "0e93691a0cfc60408b7b68f66b513c13"
    assert l_layer_r16("0e93691a0cfc60408b7b68f66b513c13") == "e6a8094fee0aa204fd97bcb0b44b8580"
    print("  ✓ L test vectors pass (4/4)")

    print("\nTesting inverse transformations (1000 deterministic inputs)...")
    import random
    random.seed(12345)
    for _ in range(1000):
        test = bytes(random.getrandbits(8) for _ in range(16)).hex()
        assert r_inv_layer(r_layer(test)) == test
        assert l_inv_layer(l_layer_r16(test)) == test
    print("  ✓ Inverse transformations pass (1000 inputs)")

    print("\nTesting S transformation...")
    assert s_layer("ffeeddccbbaa99881122334455667700") == "b66cd8887d38e8d77765aeea0c9a7efc"
    assert s_inv_layer("b66cd8887d38e8d77765aeea0c9a7efc") == "ffeeddccbbaa99881122334455667700"
    print("  ✓ S test vectors pass")

    print("\nTesting round keys...")
    key = "8899aabbccddeeff0011223344556677fedcba98765432100123456789abcdef"
    round_keys = key_schedule(key)
    expected = [
        "8899aabbccddeeff0011223344556677",
        "fedcba98765432100123456789abcdef",
        "db31485315694343228d6aef8cc78c44",
        "3d4553d8e9cfec6815ebadc40a9ffd04",
        "57646468c44a5e28d3e59246f429f1ac",
        "bd079435165c6432b532e82834da581b",
        "51e640757e8745de705727265a0098b1",
        "5a7925017b9fdd3ed72a91a22286f984",
        "bb44e25378c73123a5f32f73cdb6e517",
        "72e9dd7416bcf45b755dbaa88e4a4043",
    ]
    for i in range(10):
        assert round_keys[i] == expected[i], f"K_{i+1} mismatch: {round_keys[i]} != {expected[i]}"
    print("  ✓ Round keys K_1..K_10 pass (10/10)")

    print("\nTesting full encryption...")
    plaintext = "1122334455667700ffeeddccbbaa9988"
    ciphertext = encrypt(plaintext, round_keys)
    assert ciphertext == "7f679d90bebc24305a468d42b9d4edcd"
    print(f"  ✓ Encryption: {plaintext} -> {ciphertext}")

    print("\nTesting full decryption...")
    decrypted = decrypt(ciphertext, round_keys)
    assert decrypted == plaintext
    print(f"  ✓ Decryption: {ciphertext} -> {decrypted}")

    print("\nTesting RFC 5.6 decryption intermediates...")
    b = "7f679d90bebc24305a468d42b9d4edcd"
    k10 = round_keys[9]
    assert x_layer(k10, b) == "0d8e40e4a800d06b2f1b37ea379ead8e"
    assert l_inv_layer(x_layer(k10, b)) == "8a6b930a52211b45c5baa43ff8b91319"
    assert s_inv_layer(l_inv_layer(x_layer(k10, b))) == "76ca149eef27d1b10d17e3d5d68e5a72"
    print("  ✓ RFC 5.6 intermediates pass")

    print("\n" + "=" * 50)
    print("All RFC 7801 test vectors pass!")
    print("=" * 50)


if __name__ == "__main__":
    test_vectors()
