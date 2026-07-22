# Flux rigidity for ancient single-mode solutions of 2D hydrostatic Euler

Paper and gated exact-verification package for

> R. Jumagulov, *Flux rigidity for ancient single-mode solutions of the two-dimensional
> hydrostatic Euler equations* (2026).

**Result.** On the periodic cell `T^2`, for each harmonic `k ≥ 1` the class of exact solutions
of the 2D hydrostatic Euler equations whose `ξ`-spectrum is supported on `{0, ±k}` is **exactly
parametrized** in Fourier form: an arbitrary smooth mean profile `b(Z,τ)` transports an
arbitrary initial fluctuation by the accumulated shear, and the squared modal amplitude `S = P² + Q²` is
time-independent. The main theorem is a Liouville-type rigidity statement:

- **Flux rigidity.** Every such solution which is *ancient*, has *uniformly bounded transverse
  velocity* `v` on `τ ≤ 0`, and for which, at almost every time, *no regular stream-function
  level component is freely homotopic to a horizontal circle* `{Z=const}` (integer class ±(1,0)), has vanishing backward
  Cesàro mean flux `Φ̄ = −⟨u²v⟩̄ = 0`.
- **Mechanism.** For a time-independent mean profile the mechanism is secular; in general it
  combines an annular separation lemma with an anchored estimate on the time-averaged shear.
- **Sharpness.** Each hypothesis is necessary at every fixed `k`: explicit witnesses carry flux
  `5k²ε/8` when the `v`-bound is dropped, and a compensated-gap shear when horizontal winding
  is allowed.
- **Cylinder analogue** for time-independent mean profiles; the multi-mode problem remains open.

> **Paper:** [`paper/main.pdf`](paper/main.pdf) (17 pp; source: `paper/main.tex`).
> The PDF is compiled from exactly the arXiv submission bundle.

## Verification

None of the scripts is load-bearing for the proofs (which are complete in the paper); they pin
the computational identities and make the counterexample independently checkable (`anc/`):

    cd anc
    python3 -m venv venv && source venv/bin/activate
    pip install -r requirements.txt   # pinned: sympy 1.14.0, numpy 2.4.6
    make verify                       # exact symbolic checks + gated numeric check
    make verify-mutations             # 9-mutation corruption battery (all must be REJECTED)

Every script asserts its own checks and exits nonzero on any failure
(see `anc/README.md` for details):

| script | verifies |
|---|---|
| `cas_class_exactness.py` | exact SymPy verification of the single-mode classification (rotation law, **symbolic** harmonic `k`), conservation of `S = P² + Q²`, the velocity formula, the flux formula and cubic parity (`k = 1,2,3`), Cartesian/polar compatibility |
| `verify_flux_formula.py` | independent symbolic checks: flux formula in polar and Fourier form, parity cancellation, the exact witness integral `∮ (cos Z + cos 2Z)(1 + cos Z)² dZ = 5π/2` |
| `counterexample_check.py` | float64 grid check of the flux witness (Prop. 7.1, `ε = 0.3`) with no hand-coded derivatives (FFT collocation); 8 gates including a calibration control (`b′ ≡ 0 ⟹` zero flux through the identical pipeline) and literal reference anchors independent of the `ε` symbol |
| `mutation_tests.py` | applies nine single-site mathematical corruptions to temporary copies of the checkers (including a corruption of the calibration control itself) and asserts each corrupted checker exits nonzero; representative, not exhaustive |

Integrity chain: the top-level `SHA256SUMS` pins the paper source and `anc/SHA256SUMS`,
which pins every ancillary file (verify: `shasum -c SHA256SUMS`).

## Provenance and disclosure

This work was carried out with substantial AI assistance in derivation, computation, and
drafting. Every claim is backed by the proofs in the paper and the gated verification scripts
here; the manuscript went through twenty-one numbered rounds of independent review
(130 tracked findings; all resolved except a few optional suggestions), including cross-model
audits and a final full-protocol fresh pass. The author takes full
responsibility for the contents.

## License

MIT for the code in this repository (see `LICENSE`). The paper text is © the author.
