#V12の計算の実装は難しいので、教科書にあった解析積分の結果を使う

import sympy as sym
from scipy.integrate import nquad
from scipy.optimize import minimize_scalar
import numpy as np

z, r1, r2,=sym.symbols('z r1 r2')
#試行関数
psi = (z**3/sym.pi) * (sym.exp(-z*r1)) *(sym.exp(-z*r2))

#ハミルトニアン
def laplacian_psi(r):
     return (1/r**2 * sym.diff(r**2 * sym.diff(psi, r), r))

h_psi=(-1/sym.S(2)*laplacian_psi(r1)) + (-1/sym.S(2)*laplacian_psi(r2)) + (-2/r1*psi) + (-2/r2*psi)# + (1/r12*psi)←解析積分の結果をあとで加える

integrand=psi*h_psi * (r1**2) * (r2**2) * (2*sym.pi)**2 * 4

integrand=sym.lambdify((r1, r2, z), integrand, 'numpy')#代数計算sympyで求めたものを数値計算scipyで使えるように変換

def calculate_energy(z_val):
    options = [[0, np.inf], [0, np.inf]]
    val, err=nquad(integrand, options, args=[z_val])
    return val+5/8*z_val #V12の期待値を加えた

result = minimize_scalar(calculate_energy, bounds=(0.1, 3.0), method='bounded')
optimal_zeta = result.x
min_energy = result.fun

print("変分パラメータの最適値 : {:.5f}".format(optimal_zeta))
print("基底エネルギーの近似値 : {:.5f} hartree".format(min_energy))