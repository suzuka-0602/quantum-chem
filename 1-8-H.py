#代数計算ライブラリ
import sympy as sp
#数値計算ライブラリ
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
import numpy as np

z, r=sp.symbols('z r', real=True, positive=True)
#試行関数
def trial_function(type):
    if type=="Slater":
        return sp.sqrt(z**3 / sp.pi)*sp.exp(-z*r)
    elif type=="Gauss":
        return (2*z/sp.pi)**sp.Rational(3, 4)*sp.exp(-z*r**2)
    else:
        print("error : not deffined type")
        exit()
psi=trial_function(type="Slater")

#ハミルトニアンの作用(水素原子を考える)
h_psi = -1/sp.S(2) * (1/r**2 * sp.diff(r**2 * sp.diff(psi, r), r)) - (1/r) * psi


#被積分関数
integrand_sym=4 * sp.pi * r**2 * psi * h_psi #ヤコビアン*<psi|H|psi>

# Symbolic (記号) integration over r
E_sym = sp.integrate(integrand_sym, (r, 0, sp.oo))
E_sym = sp.simplify(E_sym)
print("E(z) =", E_sym)

# sympyで書いてみる．
dE = sp.diff(E_sym, z)
z_star = sp.solve(sp.Eq(dE, 0), z)
print("optimal z candidates =", z_star)

# Pick the physically relevant positive solution and evaluate energy
z_opt = [sol for sol in z_star if sol.is_real][0]
E_opt = sp.simplify(E_sym.subs(z, z_opt))
print("z_opt =", z_opt)
print("E_opt =", E_opt)
