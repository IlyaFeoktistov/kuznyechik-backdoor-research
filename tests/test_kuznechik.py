"""
Independent test suite for the Kuznyechik implementation.

Run with:
    cd /home/ifeoktistov/projects/pets/kuznechik && python3 -m pytest tests/test_kuznechik.py -v
    or: python3 tests/test_kuznechik.py
"""

import sys
import os
import random

# Add src to path so we can import the implementation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src', 'gost'))

import kuznechik as K


def _random_block(rng):
    return bytes(rng.getrandbits(8) for _ in range(16)).hex()


class TestGFArithmetic:
    def test_known_products(self):
        assert K.gf_mul(0x94, 0x94) == 0xA4
        assert K.gf_mul(0x01, 0x94) == 0x94
        assert K.gf_mul(0x01, 0x01) == 0x01

    def test_inverse_all_nonzero(self):
        for a in range(1, 256):
            assert K.gf_mul(a, K.gf_inv(a)) == 1

    def test_commutativity(self):
        for a in range(256):
            for b in range(256):
                assert K.gf_mul(a, b) == K.gf_mul(b, a)

    def test_zero(self):
        for a in range(256):
            assert K.gf_mul(0, a) == 0
            assert K.gf_mul(a, 0) == 0

    def test_independent_table_based(self):
        """Cross-check against a second, independently-written GF mul."""
        def gf_mul_table(a, b):
            r = 0
            for i in range(8):
                if (b >> i) & 1:
                    r ^= a << i
            for i in range(15, 7, -1):
                if (r >> i) & 1:
                    r ^= (0xC3 << (i - 8))
            return r & 0xFF

        for a in range(256):
            for b in range(256):
                assert K.gf_mul(a, b) == gf_mul_table(a, b)


class TestRAndL:
    def test_r_vectors(self):
        assert K.r_layer("00000000000000000000000000000100") == "94000000000000000000000000000001"
        assert K.r_layer("94000000000000000000000000000001") == "a5940000000000000000000000000000"
        assert K.r_layer("a5940000000000000000000000000000") == "64a59400000000000000000000000000"
        assert K.r_layer("64a59400000000000000000000000000") == "0d64a594000000000000000000000000"

    def test_l_vectors(self):
        assert K.l_layer_r16("64a59400000000000000000000000000") == "d456584dd0e3e84cc3166e4b7fa2890d"
        assert K.l_layer_r16("d456584dd0e3e84cc3166e4b7fa2890d") == "79d26221b87b584cd42fbc4ffea5de9a"
        assert K.l_layer_r16("79d26221b87b584cd42fbc4ffea5de9a") == "0e93691a0cfc60408b7b68f66b513c13"
        assert K.l_layer_r16("0e93691a0cfc60408b7b68f66b513c13") == "e6a8094fee0aa204fd97bcb0b44b8580"

    def test_r_inverse(self):
        rng = random.Random(12345)
        for _ in range(1000):
            x = _random_block(rng)
            assert K.r_inv_layer(K.r_layer(x)) == x

    def test_l_inverse(self):
        rng = random.Random(12345)
        for _ in range(1000):
            x = _random_block(rng)
            assert K.l_inv_layer(K.l_layer_r16(x)) == x

    def test_l_equals_r16(self):
        rng = random.Random(99)
        for _ in range(20):
            x = _random_block(rng)
            # manually compute R^16 and compare with l_layer_r16
            cur = x
            for _ in range(16):
                cur = K.r_layer(cur)
            assert cur == K.l_layer_r16(x)


class TestSBox:
    def test_sbox_is_permutation(self):
        assert len(K.SBOX) == 256
        assert all(0 <= b <= 255 for b in K.SBOX)
        assert len(set(K.SBOX)) == 256

    def test_sbox_inverse(self):
        assert all(K.SBOX[K.PINV[i]] == i for i in range(256))
        assert all(K.PINV[K.SBOX[i]] == i for i in range(256))

    def test_s_vectors(self):
        assert K.s_layer("ffeeddccbbaa99881122334455667700") == "b66cd8887d38e8d77765aeea0c9a7efc"
        assert K.s_inv_layer("b66cd8887d38e8d77765aeea0c9a7efc") == "ffeeddccbbaa99881122334455667700"
        assert K.s_inv_layer(K.s_layer("ffeeddccbbaa99881122334455667700")) == "ffeeddccbbaa99881122334455667700"


class TestKeySchedule:
    def test_round_keys(self):
        key = "8899aabbccddeeff0011223344556677fedcba98765432100123456789abcdef"
        round_keys = K.key_schedule(key)
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
        assert len(round_keys) == 10
        assert round_keys == expected

    def test_c_i_big_endian(self):
        assert K.c_i(1) == "6ea276726c487ab85d27bd10dd849401"
        assert K.c_i(2) == "dc87ece4d890f4b3ba4eb92079cbeb02"
        assert K.c_i(8) == "f6593616e6055689adfba18027aa2a08"


class TestEncryptDecrypt:
    def test_full_encryption_vector(self):
        key = "8899aabbccddeeff0011223344556677fedcba98765432100123456789abcdef"
        plaintext = "1122334455667700ffeeddccbbaa9988"
        round_keys = K.key_schedule(key)
        assert K.encrypt(plaintext, round_keys) == "7f679d90bebc24305a468d42b9d4edcd"

    def test_decrypt_recovers_plaintext(self):
        key = "8899aabbccddeeff0011223344556677fedcba98765432100123456789abcdef"
        plaintext = "1122334455667700ffeeddccbbaa9988"
        round_keys = K.key_schedule(key)
        ciphertext = K.encrypt(plaintext, round_keys)
        assert K.decrypt(ciphertext, round_keys) == plaintext

    def test_round_trip_random(self):
        rng = random.Random(42)
        key = "8899aabbccddeeff0011223344556677fedcba98765432100123456789abcdef"
        round_keys = K.key_schedule(key)
        for _ in range(50):
            plaintext = _random_block(rng)
            assert K.decrypt(K.encrypt(plaintext, round_keys), round_keys) == plaintext

    def test_rfc56_intermediates(self):
        key = "8899aabbccddeeff0011223344556677fedcba98765432100123456789abcdef"
        round_keys = K.key_schedule(key)
        b = "7f679d90bebc24305a468d42b9d4edcd"
        k10 = round_keys[9]
        assert K.x_layer(k10, b) == "0d8e40e4a800d06b2f1b37ea379ead8e"
        assert K.l_inv_layer(K.x_layer(k10, b)) == "8a6b930a52211b45c5baa43ff8b91319"
        assert K.s_inv_layer(K.l_inv_layer(K.x_layer(k10, b))) == "76ca149eef27d1b10d17e3d5d68e5a72"


if __name__ == "__main__":
    # Minimal runner without pytest
    import pytest
    sys.exit(pytest.main([__file__, "-v"]))
