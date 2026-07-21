#!/usr/bin/env python3
"""Symbolic (sympy) verification of the single-mode classification
(Lemma 3.1, Fourier form, general harmonic k) and the flux formula
(Lemma 3.4).

Checks (all exact symbolic identities; nonzero exit on any failure):
  (C1) rotation law, SYMBOLIC k: for
       W = b(Z,t) + P(Z,t) cos(k xi) + Q(Z,t) sin(k xi),
       the PDE residual w_t + {W,w} (w = W_xixi) equals
       -k^2 (P_t - k b_Z Q) cos(k xi) - k^2 (Q_t + k b_Z P) sin(k xi)
       identically.  This single identity carries both directions of the
       equivalence: at every point, residual = 0 iff the rotation law
       (P,Q)_t = k b_Z (Q,-P) holds there.
  (C2) solution form, SYMBOLIC k: P = P0 cos(k beta) + Q0 sin(k beta),
       Q = -P0 sin(k beta) + Q0 cos(k beta), where beta stands for dZ B,
       DEFINED by the substitution beta_t := b_Z (the only property of
       B used); C2 checks the rotated datum solves the law under that
       definition.
  (C3) modal-energy conservation: d/dt (P^2 + Q^2) == 0 under the law.
  (C4) v-formula: v = -W_Z = -b_Z - P_Z cos(k xi) - Q_Z sin(k xi).
  (C5) flux formula at k = 1, 2, 3:
       <W_Z u^2>_xi = (k^2/2) b_Z (P^2 + Q^2)  (u = W_xi);
       the T^2-normalized Phi = (k^2/4pi) oint b_Z S dZ follows by
       averaging over Z, and the 1/(4pi) factor is gated exactly in
       verify_flux_formula.py (checks (1b) and Phi_bar = 5 eps/8).
  (C6) polar compatibility, SYMBOLIC k: with
       W0' = A sin(k xi + k phi0), the composed/rotated form equals
       A sin(k xi + k phi0 + k beta).
  (C7) parity: <u^3>_xi == 0 at k = 1, 2, 3 (the conserved-energy flux
       of Remark 2.2 vanishes identically on the class).
"""
import sympy as sp

xi, Z, t = sp.symbols('xi Z t', real=True)
k = sp.Symbol('k', integer=True, positive=True)
b = sp.Function('b')(Z, t)
P = sp.Function('P')(Z, t)
Q = sp.Function('Q')(Z, t)

failures = []


def check(name, expr_zero):
    val = sp.simplify(sp.expand_trig(sp.simplify(expr_zero)))
    ok = sp.simplify(val) == 0
    print(f"{name}: residual = {val}  {'PASS' if ok else 'FAIL'}")
    if not ok:
        failures.append(name)


# ---- (C1) rotation law <=> PDE, ARBITRARY P, Q, b, SYMBOLIC k ---------------
W = b + P*sp.cos(k*xi) + Q*sp.sin(k*xi)
w = sp.diff(W, xi, 2)
bracket = sp.diff(W, xi)*sp.diff(w, Z) - sp.diff(W, Z)*sp.diff(w, xi)
resid = sp.diff(w, t) + bracket
law_form = -k**2*(sp.diff(P, t) - k*sp.diff(b, Z)*Q)*sp.cos(k*xi) \
           - k**2*(sp.diff(Q, t) + k*sp.diff(b, Z)*P)*sp.sin(k*xi)
check("(C1) residual == -k^2[(P_t - k b_Z Q)cos + (Q_t + k b_Z P)sin]",
      resid - law_form)

# ---- (C2) the rotated datum solves the law, SYMBOLIC k ----------------------
# beta(Z,t) stands for dZ B; its defining property is beta_t = b_Z.
P0 = sp.Function('P0')(Z)
Q0 = sp.Function('Q0')(Z)
beta = sp.Function('beta')(Z, t)


def sub_beta_t(e):
    return e.replace(sp.Derivative(beta, t), sp.diff(b, Z))


Psol = P0*sp.cos(k*beta) + Q0*sp.sin(k*beta)
Qsol = -P0*sp.sin(k*beta) + Q0*sp.cos(k*beta)
check("(C2a) P_t - k b_Z Q for the rotated datum",
      sub_beta_t(sp.diff(Psol, t)) - k*sp.diff(b, Z)*Qsol)
check("(C2b) Q_t + k b_Z P for the rotated datum",
      sub_beta_t(sp.diff(Qsol, t)) + k*sp.diff(b, Z)*Psol)

# ---- (C3) S conservation under the rotation law -----------------------------
check("(C3) d/dt (P^2+Q^2) under the law",
      2*P*(k*sp.diff(b, Z)*Q) + 2*Q*(-k*sp.diff(b, Z)*P))

# ---- (C4) v-formula ---------------------------------------------------------
v = -sp.diff(W, Z)
v_claim = -sp.diff(b, Z) - sp.diff(P, Z)*sp.cos(k*xi) \
          - sp.diff(Q, Z)*sp.sin(k*xi)
check("(C4) v-formula", v - v_claim)

# ---- (C5) flux formula at k = 1, 2, 3 ---------------------------------------
u = sp.diff(W, xi)
for kk in (1, 2, 3):
    avg = sp.integrate(sp.expand_trig((sp.diff(W, Z)*u**2).subs(k, kk)),
                       (xi, 0, 2*sp.pi))/(2*sp.pi)
    check(f"(C5) <W_Z u^2>_xi - (k^2/2) b_Z (P^2+Q^2) at k={kk}",
          avg - kk**2*sp.diff(b, Z)*(P**2 + Q**2)/2)

# ---- (C6) polar compatibility, SYMBOLIC k -----------------------------------
A = sp.Function('A')(Z)
phi0 = sp.Function('phi0')(Z)
Ppol = A*sp.sin(k*phi0)*sp.cos(k*beta) + A*sp.cos(k*phi0)*sp.sin(k*beta)
Qpol = -A*sp.sin(k*phi0)*sp.sin(k*beta) + A*sp.cos(k*phi0)*sp.cos(k*beta)
check("(C6) Cartesian composed form == polar form",
      Ppol*sp.cos(k*xi) + Qpol*sp.sin(k*xi)
      - A*sp.sin(k*xi + k*phi0 + k*beta))

# ---- (C7) parity: <u^3>_xi == 0 at k = 1, 2, 3 ------------------------------
for kk in (1, 2, 3):
    avg3 = sp.integrate(sp.expand_trig((u**3).subs(k, kk)),
                        (xi, 0, 2*sp.pi))/(2*sp.pi)
    check(f"(C7) <u^3>_xi at k={kk}", avg3)

if failures:
    print(f"CAS-CLASS-FAIL ({len(failures)} failed: {failures})")
    raise SystemExit(1)
print("CAS-CLASS-PASS: all classification/flux identities are exact CAS zeros")
raise SystemExit(0)
