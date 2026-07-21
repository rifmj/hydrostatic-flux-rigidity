# Ancillary verification scripts

Machine verification for "Flux rigidity for ancient single-mode solutions of
the two-dimensional hydrostatic Euler equations". None of these scripts is
load-bearing for the proofs (which are complete in the paper); they pin the
computational identities and make the counterexample independently checkable.

Run (Python >= 3.11 with sympy + numpy; reference environment: Python
3.13.13, SymPy 1.14.0, NumPy 2.4.6 -- see requirements.txt):

    make verify             # exact symbolic checks + gated numeric check
    make verify-mutations   # 9-mutation corruption battery (all must be REJECTED)

Every script asserts its checks and exits nonzero on any failure.

- cas_class_exactness.py -- exact symbolic (SymPy) verification of the
  single-mode classification (rotation law, SYMBOLIC harmonic k), conservation
  of S = P^2 + Q^2, the velocity formula, the flux formula and cubic parity
  (k = 1,2,3), and the Cartesian/polar compatibility.
- verify_flux_formula.py -- independent symbolic checks: flux formula in polar
  and Fourier form, parity cancellation, exact witness integral
  oint (cos Z + cos 2Z)(1 + cos Z)^2 dZ = 5*pi/2.
- counterexample_check.py -- float64 grid check of the flux witness
  (Proposition 7.1 of the paper, eps = 0.3) with no hand-coded derivatives
  (FFT collocation from sampled fields); 8 gates including a calibration
  control (b' == 0 => zero flux through the identical pipeline) and literal
  reference anchors (0.1875, 1.390564) independent of the eps symbol. Gate
  values are acceptance thresholds; exact residuals are platform-dependent.
- mutation_tests.py -- applies nine single-site mathematical corruptions to
  temporary copies of the checkers (including a corruption of the calibration
  control itself) and asserts each corrupted checker exits nonzero (a checker
  that accepts a corrupted certificate is not a verifier). The battery is
  representative -- one probe per checker surface class -- not an exhaustive
  corruption of every gate.
- SHA256SUMS -- checksums of the shipped files (verify: `shasum -c SHA256SUMS`).
