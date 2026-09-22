PYTHON ?= python3
export OPENBLAS_NUM_THREADS = 1
export OMP_NUM_THREADS = 1
.DEFAULT_GOAL := paper
.PHONY: discussion delivery figures paper test audit

discussion:
	$(PYTHON) python/run_discussion_v3.py

delivery:
	$(PYTHON) python/delivery_sensitivity.py

figures: discussion delivery
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
	$(PYTHON) python/presentation_manifest.py
