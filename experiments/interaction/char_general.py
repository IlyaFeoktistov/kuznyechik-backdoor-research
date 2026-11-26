"""Product-model characteristic scores (paper Sec. 7) for an arbitrary S-box,
an arbitrary MDS layer L and an arbitrary set of input multipliers M
(any F2-linear byte permutation chosen by the attacker).

  first(k,beta) = sum_i log2 m*_M((L^-1 beta e_k)_i) - 128
  2R: first + max_gamma log2 D(beta,gamma) - 8
  3R: first + log2 D(beta,gamma) - 8 + sum_i log2 d*((L gamma e_k)_i) - 128
  4R (16->1->16->1): first + log2 D(beta,gamma) - 8
        + max_{k',beta'} [ sum_i log2 D((L gamma e_k)_i, (L^-1 beta' e_k')_i) - 128 + log2 d*(beta') - 8 ]
m*_M(b) = max over offsets a of N_M(a,b); variant "nz" restricts to a != 0 (paper),
variant "all" allows a = 0.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verify"))
from gf import POLY_K, POLY_T, mul_table

MK = np.array(mul_table(POLY_K), dtype=np.int64)
MT = np.array(mul_table(POLY_T), dtype=np.int64)
X = np.arange(256)
NZ = np.arange(1, 256)


def layer_tables(LM, LINV):
    DL = np.stack([MK[LINV[:, k][None, :], NZ[:, None]] for k in range(16)])  # [k, beta-1, i]
    DF = np.stack([MK[LM[:, k][None, :], NZ[:, None]] for k in range(16)])    # [k, gamma-1, i]
    return DL, DF


def sbox_tables(S):
    D = np.stack([np.bincount(S[X ^ d] ^ S[X], minlength=256) for d in range(256)])
    with np.errstate(divide="ignore"):
        logD = np.log2(D.astype(float))
    dstar = np.log2(np.maximum(D.max(1), 1).astype(float))
    return D, logD, dstar


def tails(S, LM, LINV):
    """S-box/L dependent parts that do not depend on the multiplier."""
    DL, DF = layer_tables(LM, LINV)
    D, logD, dstar = sbox_tables(S)
    third3 = dstar[DF].sum(2) - 128                                    # [k, gamma-1]
    flatDL = DL.reshape(-1, 16)                                        # [(k',beta'), i]
    tail4 = dstar[np.tile(NZ, 16)] - 8                                 # d*(beta') - 8
    third4 = np.empty((16, 255))
    for k in range(16):
        for g0 in range(0, 255, 51):
            blk = DF[k, g0:g0 + 51]                                    # [g, i]
            s = logD[blk[:, None, :], flatDL[None, :, :]].sum(2) - 128 + tail4[None, :]
            third4[k, g0:g0 + 51] = s.max(1)
    return {"DL": DL, "logD12": logD[1:, 1:], "best2": logD[1:, 1:].max(1) - 8,
            "third3": third3, "third4": third4}


def score_multiplier(S, Mmap, tl):
    N = np.stack([np.bincount(S[Mmap ^ a] ^ S[X], minlength=256) for a in range(256)])
    out = {}
    for var, rows in (("nz", N[1:]), ("all", N)):
        m = rows.max(0)
        with np.errstate(divide="ignore"):
            mstar = np.log2(m.astype(float))
        first = mstar[tl["DL"]].sum(2) - 128                           # [k, beta-1]
        s2 = float((first + tl["best2"][None, :]).max())
        mid = tl["logD12"] - 8
        s3 = max(float((first[k][:, None] + mid + tl["third3"][k][None, :]).max()) for k in range(16))
        s4 = max(float((first[k][:, None] + mid + tl["third4"][k][None, :]).max()) for k in range(16))
        out[var] = (s2, s3, s4)
    return out


def best_over(S, maps, tl):
    """maps: dict label -> 256-array. Returns per variant/round the best score and its label."""
    res = {}
    for lab, Mmap in maps.items():
        sc = score_multiplier(S, Mmap, tl)
        for var, triple in sc.items():
            for r, v in zip((2, 3, 4), triple):
                key = f"{var}{r}"
                if key not in res or v > res[key][0]:
                    res[key] = (v, lab)
    return res


def k_multipliers():
    return {f"K{c:02x}": MK[c] for c in range(2, 256)}


def t_multipliers(Blin=None):
    """Hidden-field multipliers; for an affine copy A o pi o B conjugate by B's linear part."""
    if Blin is None:
        return {f"T{c:02x}": MT[c] for c in range(2, 256)}
    Binv = np.argsort(Blin)
    return {f"T{c:02x}": Binv[MT[c][Blin]] for c in range(2, 256)}
