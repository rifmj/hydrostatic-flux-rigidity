PY ?= python3
.PHONY: verify verify-mutations
verify:
	$(PY) cas_class_exactness.py
	$(PY) verify_flux_formula.py
	$(PY) counterexample_check.py
verify-mutations:
	$(PY) mutation_tests.py
