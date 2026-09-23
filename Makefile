PYTHON ?= python3
export OPENBLAS_NUM_THREADS = 1
export OMP_NUM_THREADS = 1
.DEFAULT_GOAL := paper
.PHONY: discussion delivery patch-study physics reconcile figures paper test audit

discussion:
	$(PYTHON) python/run_discussion_v3.py

delivery:
	$(PYTHON) python/delivery_sensitivity.py

patch-study:
	$(PYTHON) python/dynamics_followup_analysis.py

physics:
	$(PYTHON) python/local_ordering.py
	$(PYTHON) python/physics_sweep.py

reconcile:
	$(PYTHON) python/workshop_reconciliation.py

figures: discussion delivery patch-study
	$(PYTHON) python/physics_sweep.py --plot
	$(PYTHON) python/plot_capacity.py
	$(PYTHON) python/draw_architectures.py

paper: figures
	$(PYTHON) python/build_manuscript.py
	$(PYTHON) python/presentation_manifest.py --write

test:
	$(PYTHON) -m pytest -q --junitxml=results/presentation_v1/tests.xml

audit:
	$(PYTHON) python/audit_refinement_v3.py
	$(PYTHON) python/check_discussion_v3.py
	$(PYTHON) python/delivery_sensitivity.py --check
	$(PYTHON) python/large_patch_analysis.py --check
	$(PYTHON) python/dynamics_followup_analysis.py --check
	$(PYTHON) python/physics_sweep.py --check
	$(PYTHON) python/workshop_reconciliation.py --check
	$(PYTHON) python/presentation_manifest.py
