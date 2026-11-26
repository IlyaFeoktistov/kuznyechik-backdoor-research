"""B4: eigenvalues of the Kuznyechik LFSR matrix R (L = R^16) over GF(2^8)_K.

Finding: the characteristic polynomial of R splits into 16 distinct linear
factors whose roots are beta^120..beta^135, beta = alpha^64 = 0xa0, i.e. the
BCH-style recursive MDS construction (Augot-Finiasz, FSE 2014). This explains
the palindromic coefficients and the order-17 eigenvalues (dim ker(L^17+I)=16).
Controls: 2000 palindromic-shape and 2000 fully random coefficient vectors.
Usage (repo root): python3 experiments/backdoor/linear_eigen.py   (~6 min)
"""
import sys, numpy as np
sys.path.insert(0,"experiments/verify"); sys.path.insert(0,"src/gost"); sys.path.insert(0,"experiments/backdoor")
import kuznechik as KZ
from gf import POLY_K, POLY_T, mul_table, log_tables
from layer_keys_audit import f2_matrix, kernel_dim
MK=np.array(mul_table(POLY_K)); expK,logK=log_tables(POLY_K)
def eig(coeffs):
    R=np.zeros((16,16),dtype=np.int64); R[0,:]=coeffs
    for i in range(1,16): R[i,i-1]=1
    out=[]
    for lam in range(1,256):
        M=R.copy()
        for i in range(16): M[i,i]^=lam
        d=kernel_dim(f2_matrix(M))//8
        if d: out.append((lam,d,255//np.gcd(255,logK[lam])))
    return out

print("R eigenvalues (lambda, geom mult, order):", eig(KZ.L_COEFF))
# controls: random coefficient vectors with the same palindromic shape a,b,c,d,e,f,1,g,1,f,e,d,c,b,a,1
rng=np.random.default_rng(51); hits17=0; anyroot=0; N=2000
for _ in range(N):
    a=[int(v) for v in rng.integers(1,256,7)]
    co=[a[0],a[1],a[2],a[3],a[4],a[5],1,a[6],1,a[5],a[4],a[3],a[2],a[1],a[0],1]
    e=eig(co); anyroot+=bool(e); hits17+=sum(1 for _,_,o in e if o==17)>=2
print(f"palindromic-shape controls n={N}: any eigenvalue {anyroot/N:.3f}, >=2 eigenvalues of order 17: {hits17/N:.4f}")
rng=np.random.default_rng(52); h=0; ar=0
for _ in range(N):
    co=[int(v) for v in rng.integers(1,256,16)]
    e=eig(co); ar+=bool(e); h+=sum(1 for _,_,o in e if o==17)>=2
print(f"fully random coefficients n={N}: any eigenvalue {ar/N:.3f}, >=2 of order 17: {h/N:.4f}")
import json, math
expK_, logK_ = expK, logK
Lg = sorted(logK[l] for l, _, _ in eig(KZ.L_COEFF))
cons = [(s, l) for s in range(1, 255) if math.gcd(s, 255) == 1 for l in range(255)
        if sorted((s * (l + i)) % 255 for i in range(16)) == Lg]
print("consecutive-power representations (s, l): beta = alpha^s, roots beta^l..beta^(l+15):", cons)
json.dump({"eigen_logs": Lg, "consecutive": cons,
           "palindromic_controls": {"n": N, "seed": 51, "any_eig": anyroot / N, "ge2_order17": hits17 / N},
           "random_controls": {"n": N, "seed": 52, "any_eig": ar / N, "ge2_order17": h / N}},
          open("results/backdoor/linear_eigen.json", "w"), indent=1)
