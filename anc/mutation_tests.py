#!/usr/bin/env python3
"""Mutation battery for the ancillary checkers (corrupted-certificate test).

Each mutation applies a single mathematical corruption to a TEMPORARY copy of
one checker and asserts that the corrupted checker exits NONZERO. A checker
that accepts a corrupted certificate is not a verifier; this driver makes the
battery reproducible (`make verify-mutations`).

Mutations (checker: corruption -> expected catching gate):
  M1 cas: flux constant /2 -> /3            -> C5 residual nonzero
  M2 cas: sign-flipped rotation law          -> C1 residual nonzero
  M3 cx : literal flux target 0.1875 -> *1.01 -> gate (2b) literal anchor
  M4 cx : phase speed 1.01*b'                -> gate (1) residual (inconsistent
          field) and gate (3) slope
  M5 cx : slope prediction corrupted to the historical wrong-b'' formula
          -(EPS)*sin(Z)*(1+3*cos(Z))         -> gate (3) fit-vs-prediction
  M6 cx : EPS 0.3 -> 0.31 (parameter drift)  -> gates (2b)/(3b) literal anchors
  M7 vff: polar flux prefactor 1/2 -> 1/3    -> check (1) exact identity
  M8 cas: v-formula sign flip (+P_Z cos)     -> C4 residual nonzero
  M9 cx : calibration control corrupted (scale_b 0.0 -> 0.01: the control
          no longer zeroes the mean shear)   -> gate (0) calibration
The battery is representative (one probe per checker surface class), not an
exhaustive corruption of every gate.
Exit 0 iff every mutation is rejected (nonzero exit of the mutant).
"""
import pathlib
import subprocess
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
CAS = (HERE / "cas_class_exactness.py").read_text()
CX = (HERE / "counterexample_check.py").read_text()
VFF = (HERE / "verify_flux_formula.py").read_text()

MUTS = [
    ("M1-cas-flux-third", CAS,
     "kk**2*sp.diff(b, Z)*(P**2 + Q**2)/2",
     "kk**2*sp.diff(b, Z)*(P**2 + Q**2)/3"),
    ("M2-cas-law-sign", CAS,
     "- k**2*(sp.diff(Q, t) + k*sp.diff(b, Z)*P)*sp.sin(k*xi)",
     "+ k**2*(sp.diff(Q, t) + k*sp.diff(b, Z)*P)*sp.sin(k*xi)"),
    ("M3-cx-target-drift", CX, "TARGET_FLUX = 0.1875\nSLOPE_REF",
     "TARGET_FLUX = 0.1875*1.01\nSLOPE_REF"),
    ("M4-cx-phase-speed", CX,
     "np.sin(XI + scale_b*bp_f(ZZ)*tau)", "np.sin(XI + 1.01*scale_b*bp_f(ZZ)*tau)"),
    ("M5-cx-wrong-bpp-prediction", CX,
     "slope_s = sp.Abs(A_s * sp.diff(b_s, Zs, 2))",
     "slope_s = sp.Abs(A_s * (-EPS*sp.sin(Zs)*(1+3*sp.cos(Zs))))"),
    ("M6-cx-eps-drift", CX, "EPS = 0.3\n", "EPS = 0.31\n"),
    ("M7-vff-polar-third", VFF,
     "rhs = sp.Rational(1, 2) * sp.diff(b, Z) * A ** 2",
     "rhs = sp.Rational(1, 3) * sp.diff(b, Z) * A ** 2"),
    ("M8-cas-vformula-sign", CAS,
     "v_claim = -sp.diff(b, Z) - sp.diff(P, Z)*sp.cos(k*xi)",
     "v_claim = -sp.diff(b, Z) + sp.diff(P, Z)*sp.cos(k*xi)"),
    ("M9-cx-calibration-corrupt", CX,
     "phibar0 = np.mean([flux_at(t, scale_b=0.0) for t in np.linspace(0, 100, 11)])",
     "phibar0 = np.mean([flux_at(t, scale_b=0.01) for t in np.linspace(0, 100, 11)])"),
]

failures = []
with tempfile.TemporaryDirectory() as tmp:
    for name, src, old, new in MUTS:
        if old not in src:
            print(f"{name}: MUTATION SITE NOT FOUND (checker drifted)  FAIL")
            failures.append(name)
            continue
        mut = pathlib.Path(tmp) / f"{name}.py"
        mut.write_text(src.replace(old, new, 1))
        r = subprocess.run([sys.executable, str(mut)], capture_output=True,
                           text=True, timeout=900)
        rejected = r.returncode != 0
        tail = [l for l in r.stdout.strip().splitlines() if "FAIL" in l][-1:]
        print(f"{name}: {'REJECTED (exit %d)' % r.returncode if rejected else 'ACCEPTED -- BATTERY FAILURE'}"
              f"  {tail[0] if tail else ''}")
        if not rejected:
            failures.append(name)

if failures:
    print(f"MUTATION-BATTERY-FAIL ({len(failures)}: {failures})")
    raise SystemExit(1)
print("MUTATION-BATTERY-PASS: all 9 corrupted checkers rejected with nonzero exit")
raise SystemExit(0)
