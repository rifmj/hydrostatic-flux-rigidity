#!/usr/bin/env python3
"""Independent symbolic checks for the paper (fresh; no imports from the
research repo):

  (1) Lemma 3.4 (flux formula, polar sector): for W = b(Z,t) + A(Z) sin(theta),
      theta = xi + phi0(Z) + psi(Z,t):
        <W_Z u^2>_xi = (1/2) b_Z A^2      (u = W_xi)
      -- the xi-average kills every oscillatory term.
  (1b) The same identity in Fourier form: W = b + P cos xi + Q sin xi gives
        <W_Z u^2>_xi = (1/2) b_Z (P^2+Q^2).
  (2) Triple-channel parity: <cos^3> = <sin cos^2> = 0 over a period.
  (3) Proposition 7.1 integral:
        I = ∮ (cos Z + cos 2Z)(1 + cos Z)^2 dZ = 5*pi/2,
      hence Phi_bar = (eps/(4*pi)) * I = 5*eps/8.
Exit 0 iff all three are exact zeros/equalities."""
import sympy as sp

xi, Z, t, eps = sp.symbols("xi Z t epsilon", real=True)
b = sp.Function("b")(Z, t)
A = sp.Function("A")(Z)
phi0 = sp.Function("phi0")(Z)
psi = sp.Function("psi")(Z, t)
theta = xi + phi0 + psi

W = b + A * sp.sin(theta)
u = sp.diff(W, xi)
WZ = sp.diff(W, Z)

avg = lambda f: sp.integrate(f, (xi, 0, 2 * sp.pi)) / (2 * sp.pi)

lhs = sp.simplify(avg(sp.expand_trig(WZ * u ** 2)))
rhs = sp.Rational(1, 2) * sp.diff(b, Z) * A ** 2
c1 = sp.simplify(lhs - rhs) == 0

# (1b) Cartesian/Fourier form of the same identity (Lemma 3.1/3.4):
#      W = b + P cos xi + Q sin xi  =>  <W_Z u^2>_xi = (1/2) b_Z (P^2+Q^2)
Pc = sp.Function("P")(Z, t)
Qc = sp.Function("Q")(Z, t)
Wc = b + Pc * sp.cos(xi) + Qc * sp.sin(xi)
uc = sp.diff(Wc, xi)
lhs_c = sp.simplify(avg(sp.expand_trig(sp.diff(Wc, Z) * uc ** 2)))
rhs_c = sp.Rational(1, 2) * sp.diff(b, Z) * (Pc ** 2 + Qc ** 2)
c1b = sp.simplify(lhs_c - rhs_c) == 0

c2 = (sp.integrate(sp.cos(xi) ** 3, (xi, 0, 2 * sp.pi)) == 0
      and sp.integrate(sp.sin(xi) * sp.cos(xi) ** 2,
                       (xi, 0, 2 * sp.pi)) == 0)

I = sp.integrate((sp.cos(Z) + sp.cos(2 * Z)) * (1 + sp.cos(Z)) ** 2,
                 (Z, 0, 2 * sp.pi))
c3 = sp.simplify(I - sp.Rational(5, 2) * sp.pi) == 0
phibar = sp.simplify(eps / (4 * sp.pi) * I)
c4 = sp.simplify(phibar - sp.Rational(5, 8) * eps) == 0

print("(1)  <W_Z u^2>_xi == (1/2) b_Z A^2       :", c1)
print("(1b) Cartesian: == (1/2) b_Z (P^2+Q^2)   :", c1b)
print("(2)  triple-channel parity zeros         :", c2)
print("(3)  counterexample integral == 5*pi/2   :", c3, "| I =", I)
print("     Phi_bar == 5*eps/8                  :", c4)
ok = c1 and c1b and c2 and c3 and c4
print("ALL-CHECKS-", "PASS" if ok else "FAIL", sep="")
raise SystemExit(0 if ok else 1)
