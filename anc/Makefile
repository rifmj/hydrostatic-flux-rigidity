PY ?= python3
.PHONY: verify verify-fast verify-mutations
verify:
	$(PY) cas_class_exactness.py
	$(PY) verify_flux_formula.py
	$(PY) counterexample_check.py
verify-fast:
	$(PY) cas_class_exactness.py
	$(PY) verify_flux_formula.py
	$(PY) counterexample_check.py --fast
verify-mutations:
	$(PY) mutation_tests.py
