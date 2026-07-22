# Ancillary verification scripts

Machine verification for "Flux rigidity for ancient single-mode solutions of
the two-dimensional hydrostatic Euler equations". None of these scripts is
load-bearing for the proofs (which are complete in the paper); they pin the
computational identities and make the counterexample independently checkable.

Run (Python >= 3.11 with sympy + numpy; reference environment: Python
3.13.13, SymPy 1.14.0, NumPy 2.4.6 -- see requirements.txt):

    make verify             # exact symbolic checks + gated numeric check
    make verify-fast        # same, with counterexample_check.py --fast (~5 s total)
    make verify-mutations   # 9-mutation corruption battery (all must be REJECTED)

Every script asserts its checks and exits nonzero on any failure.

Timing and resources (reference environment above; the exact-symbolic
scripts are machine-independent, the float64 script is not):

    cas_class_exactness.py   ~3 s     symbolic, negligible memory
    verify_flux_formula.py   ~1 s     symbolic, negligible memory
    counterexample_check.py  ~15-40 s ~0.6 GB peak (one 256x4096 FFT loop)
    counterexample_check.py --fast  ~3 s   coarser slope grid; all gates pass
    make verify-mutations    ~2 min   (runs the counterexample script 7x)

On a slow machine the full `counterexample_check.py` can take several
minutes; use `--fast` (identical gates, slope grid 256x1024 instead of
256x4096) for a quick pass. Expected reference values (full run):
calibration |Phi_bar| < 1e-13 (env 7.8e-21); PDE residual ~1.5e-7;
time-mean flux = 5*eps/8 = 0.1875 (rel ~1e-16); fitted secular slope
1.390537 vs predicted 1.390564 (rel ~2e-5). The counterexample is proved
analytically in the paper; these numerics are auxiliary diagnostics.

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
  Accepts `--fast` (coarser slope grid; same gates) for a quick pass.
- mutation_tests.py -- applies nine single-site mathematical corruptions to
  temporary copies of the checkers (including a corruption of the calibration
  control itself) and asserts each corrupted checker exits nonzero (a checker
  that accepts a corrupted certificate is not a verifier). The battery is
  representative -- one probe per checker surface class -- not an exhaustive
  corruption of every gate.
- SHA256SUMS -- checksums of the shipped files. Verify with
  `sha256sum -c SHA256SUMS` (GNU/Linux) or `shasum -a 256 -c SHA256SUMS`
  (BSD/macOS); the explicit `-a 256` avoids the SHA-1 default of older
  `shasum`.
