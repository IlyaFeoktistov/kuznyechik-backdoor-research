"""
S-box analysis toolkit for Kuznyechik (GOST R 34.12-2015).

Computes the standard security metrics for an 8-bit S-box (permutation of
GF(2^8) represented as a 256-byte table):

  * algebraic degree (max ANF degree over the 8 output bits)
  * differential distribution table (DDT) + max differential probability
  * linear approximation table (LAT) / Walsh spectrum + nonlinearity
  * partition invariants (sizes of the partition of the input cube induced by
    each output bit)
  * nonlinearity (NL)

All functions take the S-box as a sequence of 256 ints (0..255).
"""

import math
from collections import defaultdict


# ---------------------------------------------------------------------------
# ANF / algebraic degree
# ---------------------------------------------------------------------------

def anf_coefficients(bit_fn):
    """Return the ANF coefficients c_S for every subset S of {0..7}.

    bit_fn: a function mapping an 8-bit integer x (0..255) to a bit (0/1).
    c_S = XOR_{T subset S} bit_fn(T)   (Möbius transform over the cube).
    The coefficient of the monomial prod_{i in S} x_i is c_S.
    """
    f = [bit_fn(x) & 1 for x in range(256)]
    c = f[:]  # copy
    for i in range(8):
        for x in range(256):
            if (x >> i) & 1:
                c[x] ^= c[x ^ (1 << i)]
    return c


def anf_degree(bit_fn):
    """Maximum number of input variables in any nonzero ANF monomial."""
    c = anf_coefficients(bit_fn)
    for size in range(8, -1, -1):
        # any subset of exactly `size` bits with c_S == 1
        for mask in range(256):
            if c[mask] and mask.bit_count() == size:
                return size
    return 0


def sbox_anf_degrees(sbox):
    """Return the ANF degree of each of the 8 output bits, plus the max."""
    degs = []
    for out_bit in range(8):
        def f(x, ob=out_bit):
            return (sbox[x] >> ob) & 1
        degs.append(anf_degree(f))
    return degs, max(degs)


# ---------------------------------------------------------------------------
# Differential distribution table (DDT)
# ---------------------------------------------------------------------------

def ddt(sbox):
    """Return the DDT as a flat list of 256*256 ints and the (d, e) pairs.

    ddt[d][e] = #{ x : sbox[x] XOR sbox[x XOR d] == e }.
    Returns (table, max_val, d_list, e_list) where table is a dict keyed by
    (d, e), plus the max entry value and the (d, e) achieving it.
    """
    table = [[0] * 256 for _ in range(256)]
    for d in range(256):
        for x in range(256):
            e = sbox[x] ^ sbox[(x ^ d) & 0xFF]
            table[d][e] += 1
    max_val = 0
    max_pair = None
    for d in range(256):
        for e in range(256):
            if table[d][e] > max_val:
                max_val = table[d][e]
                max_pair = (d, e)
    return table, max_val, max_pair


def ddt_summary(sbox):
    """Return a dict of DDT statistics."""
    table, max_val, max_pair = ddt(sbox)
    # distribution of the DDT entries (how often each count value occurs)
    entry_dist = defaultdict(int)
    for d in range(256):
        for e in range(256):
            entry_dist[table[d][e]] += 1
    # max differential probability for d != 0
    max_d = 0
    for d in range(1, 256):
        for e in range(256):
            if table[d][e] > max_d:
                max_d = table[d][e]
    return {
        "max_value": max_val,
        "max_pair": max_pair,
        "max_for_d_nonzero": max_d,
        "entry_distribution": dict(sorted(entry_dist.items())),
    }


# ---------------------------------------------------------------------------
# Linear approximation table / Walsh spectrum / nonlinearity
# ---------------------------------------------------------------------------

def walsh_spectrum(sbox):
    """For each of the 8 output bits, compute the Walsh spectrum over all
    256 input linear masks. Returns a list of 8 lists, each of length 256.

    Walsh_l,u = sum_x (-1)^{ u·x + l·sbox(x) }
    """
    spectrum = []
    for l in range(8):
        spec = [0] * 256
        for u in range(256):
            s = 0
            for x in range(256):
                val = ((u & x).bit_count() + ((sbox[x] >> l) & 1)) & 1
                s += -1 if val else 1
            spec[u] = s
        spectrum.append(spec)
    return spectrum


def nonlinearity(sbox):
    """Nonlinearity = min over output bits l, input masks u of
    | 2^8 - |Walsh_l,u| | / 2  ... equivalently (256 - max|Walsh|) / 2."""
    spec = walsh_spectrum(sbox)
    min_nl = 10**9
    for l in range(8):
        for u in range(256):
            nl = (256 - abs(spec[l][u])) // 2
            if nl < min_nl:
                min_nl = nl
    return min_nl


def lat_summary(sbox):
    """Return summary statistics of the LAT (per output bit)."""
    spec = walsh_spectrum(sbox)
    per_bit = []
    for l in range(8):
        vals = [abs(v) for v in spec[l]]
        per_bit.append({
            "max_abs": max(vals),
            "nonlinearity": (256 - max(vals)) // 2,
            "num_max": sum(1 for v in vals if v == max(vals)),
            "num_zero": sum(1 for v in spec[l] if v == 0),
        })
    return per_bit


# ---------------------------------------------------------------------------
# Partition invariants
# ---------------------------------------------------------------------------

def partition_sizes(sbox, use_output=True):
    """For each output bit l, compute the multiset of partition cell sizes of
    the partition of the input cube induced by the l-th output bit.

    The l-th output bit splits the 256 inputs into two cells (0 / 1). But to get
    a richer invariant, we consider the partition by the full output value:
    each output byte defines a partition into up to 256 cells (the preimage of
    each output value). Since the S-box is a permutation, every cell has size 1.

    To get a non-trivial partition, we use the *linear* partition: for each
    output bit l, consider the two sets {x : bit_l(sbox(x))=0} and =1.

    We return, for each l, the sizes of the two cells.
    """
    sizes = []
    for l in range(8):
        c0 = sum(1 for x in range(256) if ((sbox[x] >> l) & 1) == 0)
        c1 = 256 - c0
        sizes.append((c0, c1))
    return sizes


def linear_partition(sbox):
    """For every non-zero input mask u and every output bit l, count how many
    inputs map to a given output bit. This gives the balance of each output
    bit across affine subcubes. Returns a list of (u, l, count0, count1).

    A simpler, more useful invariant: for each output bit l, the number of
    inputs x with sbox(x) having bit l set. Already covered by partition_sizes.

    We also compute, for each pair (u, l), the Walsh coefficient which measures
    how balanced the l-th output bit is on the coset structure of u.
    """
    return walsh_spectrum(sbox)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_sbox(path="data/original/pi.bin"):
    with open(path, "rb") as f:
        return list(f.read())


def main():
    sbox = load_sbox()
    print("=== Kuznyechik S-box baseline ===")

    degs, max_deg = sbox_anf_degrees(sbox)
    print(f"ANF degrees per output bit: {degs}")
    print(f"Max ANF degree: {max_deg}")

    ds = ddt_summary(sbox)
    print(f"DDT max value: {ds['max_value']} (pair {ds['max_pair']})")
    print(f"DDT max for d!=0: {ds['max_for_d_nonzero']}")
    print(f"DDT entry distribution: {ds['entry_distribution']}")

    nl = nonlinearity(sbox)
    print(f"Nonlinearity: {nl}")

    lat = lat_summary(sbox)
    print("LAT per output bit:")
    for l, v in enumerate(lat):
        print(f"  bit {l}: max|Walsh|={v['max_abs']}, NL={v['nonlinearity']}, "
              f"num_max={v['num_max']}, num_zero={v['num_zero']}")

    print("Partition sizes (per output bit):")
    for l, (c0, c1) in enumerate(partition_sizes(sbox)):
        print(f"  bit {l}: c0={c0}, c1={c1}")


if __name__ == "__main__":
    main()
