"""Independent GF(2^8) arithmetic (does not import the cipher implementation).

Two different algorithms are provided so results can be cross-checked:
  * mul_clmul: carry-less multiply + polynomial long division;
  * log/antilog tables built from a primitive element.
"""

POLY_T = 0x11D  # hidden TKlog field, X^8+X^4+X^3+X^2+1 (Perrin 2019)
POLY_K = 0x1C3  # Kuznyechik linear-layer field, X^8+X^7+X^6+X+1 (RFC 7801)


def mul_clmul(a, b, poly):
    """Carry-less product, then reduction modulo `poly` by long division."""
    p = 0
    for i in range(8):
        if (b >> i) & 1:
            p ^= a << i
    for deg in range(14, 7, -1):
        if (p >> deg) & 1:
            p ^= poly << (deg - 8)
    return p


def mul_table(poly):
    """Full 256x256 multiplication table as a list of lists."""
    return [[mul_clmul(a, b, poly) for b in range(256)] for a in range(256)]


def log_tables(poly, gen=0x02):
    """exp[i] = gen^i, log[x] = i (ordinary discrete log, log[1] = 0)."""
    exp = [0] * 255
    log = [-1] * 256
    x = 1
    for i in range(255):
        exp[i] = x
        log[x] = i
        x = mul_clmul(x, gen, poly)
    assert x == 1 and -1 not in log[1:], "generator is not primitive"
    return exp, log


def f2_rank(vectors):
    """Rank over F2 of a list of 8-bit integers (Gaussian elimination)."""
    basis = []  # kept reduced: each element has a unique leading bit
    for v in vectors:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    return len(basis)


def f2_span(vectors):
    span = {0}
    for v in vectors:
        span |= {s ^ v for s in span}
    return span
