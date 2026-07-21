#!/usr/bin/env python3
"""Grid check of the counterexample of Proposition 7.1 --- the
'flux costs unbounded v' witness

    W = b(Z) + A(Z) sin(xi + b'(Z) tau),
    A = 1 + cos Z,   b = eps * A * sin Z,   here eps = 0.3
    (Fourier datum (P0, Q0) = (0, A); S = A^2).

Design rule: NO hand-coded derivatives anywhere. W is sampled from its
closed form; every space derivative is computed by FFT collocation from the
SAMPLES (spectrally accurate for these smooth periodic fields), the time
derivative by centered finite differences of sampled snapshots, and the
predicted secular slope max_Z |A b''| by sympy differentiation of b, and the
flux/slope gates are ADDITIONALLY anchored to the literal reference
constants TARGET_FLUX = 0.1875 and SLOPE_REF = 1.390564, which are
independent of the EPS symbol -- so a wrong construction coefficient breaks
the residual, flux, or slope gate, and a drift of EPS itself is rejected by
the literal anchors.

Checks (floating-point spot check at eps = 0.3; nonzero exit on failure):
  (0) calibration control (calibrate-before-news): the identical flux
      pipeline with b' == 0 (b = 0, same A) returns |Phi_bar| < 1e-13;
  (1) PDE residual w_tau + {W,w} (w = W_xixi), all derivatives from
      sampled W only:  max < 1e-6  (measured ~1.5e-7; floor set by the
      tau finite-difference/roundoff balance, not by the solution);
  (2) time-mean flux <W_Z u^2> == 5 eps/8 AND == the literal 0.1875,
      relative error < 1e-10 each; the flux is also time-independent
      (spread < 1e-12), as the static class requires;
  (3) sup|v(tau)| grows linearly (fit-residual gate), fitted slope ==
      max_Z |A b''| (sympy) to relative error < 1e-2, and the sympy
      prediction itself matches the literal 1.390564 to 1e-5.
Confinement/topology is proved analytically in the paper (Prop. 7.1) and is
NOT tested here.
"""
import numpy as np
import sympy as sp

EPS = 0.3
# Literal reference constants -- independent of the EPS symbol (anti
# self-reference anchors): 5*0.3/8 and max_Z |A b''| at eps = 0.3.
TARGET_FLUX = 0.1875
SLOPE_REF = 1.390564
failures = []


def gate(name, ok, detail):
    print(f"{name}: {detail}  {'PASS' if ok else 'FAIL'}")
    if not ok:
        failures.append(name)


# ---- symbolic ground truth (W's closed form + slope prediction only) --------
Zs = sp.symbols('Z', real=True)
A_s = 1 + sp.cos(Zs)
b_s = EPS * A_s * sp.sin(Zs)
bp_s = sp.diff(b_s, Zs)            # b' (enters W's phase)
slope_s = sp.Abs(A_s * sp.diff(b_s, Zs, 2))   # |A b''| -- prediction only
A_f, b_f, bp_f = (sp.lambdify(Zs, e, 'numpy') for e in (A_s, b_s, bp_s))
zfine = np.linspace(0, 2*np.pi, 200001)
slope_pred = float(np.max(sp.lambdify(Zs, slope_s, 'numpy')(zfine)))


def make_grid(gxi, gz):
    xi = np.linspace(0, 2*np.pi, gxi, endpoint=False)
    z = np.linspace(0, 2*np.pi, gz, endpoint=False)
    XI, ZZ = np.meshgrid(xi, z, indexing='ij')
    kxi = np.fft.fftfreq(gxi, d=1.0/gxi)
    kz = np.fft.fftfreq(gz, d=1.0/gz)
    return XI, ZZ, kxi, kz


def d_xi(F, kxi):
    return np.real(np.fft.ifft(1j*kxi[:, None]*np.fft.fft(F, axis=0), axis=0))


def d_z(F, kz):
    return np.real(np.fft.ifft(1j*kz[None, :]*np.fft.fft(F, axis=1), axis=1))


def W_sample(XI, ZZ, tau, scale_b=1.0):
    return scale_b*b_f(ZZ) + A_f(ZZ)*np.sin(XI + scale_b*bp_f(ZZ)*tau)


# ---- (0) calibration control: b' == 0 => zero flux, SAME pipeline -----------
XI, ZZ, kxi, kz = make_grid(256, 256)


def flux_at(tau, scale_b=1.0):
    W = W_sample(XI, ZZ, tau, scale_b)
    return float(np.mean(d_z(W, kz) * d_xi(W, kxi)**2))


phibar0 = np.mean([flux_at(t, scale_b=0.0) for t in np.linspace(0, 100, 11)])
gate("(0) calibration b'==0 -> zero flux", abs(phibar0) < 1e-13,
     f"|Phi_bar| = {abs(phibar0):.2e}")


# ---- (1) PDE residual from sampled W only -----------------------------------
def residual_at(tau, dt=2e-4):
    W = W_sample(XI, ZZ, tau)
    w = d_xi(d_xi(W, kxi), kxi)
    w_plus = d_xi(d_xi(W_sample(XI, ZZ, tau+dt), kxi), kxi)
    w_minus = d_xi(d_xi(W_sample(XI, ZZ, tau-dt), kxi), kxi)
    w_tau = (w_plus - w_minus)/(2*dt)
    res = w_tau + d_xi(W, kxi)*d_z(w, kz) - d_z(W, kz)*d_xi(w, kxi)
    return float(np.max(np.abs(res)))


r = max(residual_at(t) for t in (0.0, 3.7, 25.0))
gate("(1) PDE residual (sampled-field check)", r < 1e-6, f"max {r:.2e}")

# ---- (2) time-mean flux == 5 eps/8 ------------------------------------------
taus = np.linspace(0, 100, 401)
phis = np.array([flux_at(t) for t in taus])
target = 5*EPS/8
rel = abs(phis.mean() - target)/abs(target)
gate("(2) time-mean flux == 5 eps/8", rel < 1e-10,
     f"{phis.mean():+.12f} vs {target:+.12f} (rel {rel:.1e})")
rel_lit = abs(phis.mean() - TARGET_FLUX)/TARGET_FLUX
gate("(2b) flux == literal 0.1875", rel_lit < 1e-10,
     f"rel {rel_lit:.1e} vs literal anchor")
spread = float(np.max(np.abs(phis - phis.mean())))
gate("(2c) flux time-independent", spread < 1e-12,
     f"max|Phi(tau) - mean| = {spread:.1e}")

# ---- (3) secular slope: fitted vs sympy prediction --------------------------
XI2, ZZ2, kxi2, kz2 = make_grid(256, 4096)
taus2 = np.linspace(0, 400, 401)
vmax = np.empty(len(taus2))
for i, t in enumerate(taus2):
    W = W_sample(XI2, ZZ2, t)
    vmax[i] = np.max(np.abs(-d_z(W, kz2)))
sel = taus2 >= 200
coef = np.polyfit(taus2[sel], vmax[sel], 1)
k_fit = coef[0]
fit_res = float(np.sqrt(np.mean((vmax[sel] - np.polyval(coef, taus2[sel]))**2))
                / np.mean(vmax[sel]))
rel_k = abs(k_fit - slope_pred)/slope_pred
gate("(3) sup|v| slope == max|A b''|", rel_k < 1e-2,
     f"fit {k_fit:.6f} vs predicted {slope_pred:.6f} (rel {rel_k:.1e})")
gate("(3b) sympy prediction == literal 1.390564",
     abs(slope_pred - SLOPE_REF)/SLOPE_REF < 1e-5,
     f"pred {slope_pred:.6f} vs literal anchor")
gate("(3c) growth is linear (fit residual)", fit_res < 1e-3,
     f"relative RMS fit residual {fit_res:.1e}")

if failures:
    print(f"COUNTEREXAMPLE-FAIL ({len(failures)} failed: {failures})")
    raise SystemExit(1)
print("COUNTEREXAMPLE-PASS: the witness carries flux 5eps/8 and violates the "
      "v-bound (linear growth at the predicted slope); topology not tested.")
raise SystemExit(0)
