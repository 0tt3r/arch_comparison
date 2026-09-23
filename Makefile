PYTHON ?= python3
CUDAQ_PYTHON ?= .venv-cudaq/bin/python
QEC_PYTHON ?= /home/ac.otten/soft/mambaforge/envs/stim/bin/python
export OPENBLAS_NUM_THREADS = 1
export OMP_NUM_THREADS = 1
.DEFAULT_GOAL := paper
.PHONY: discussion delivery patch-study physics reconcile figures paper test audit implementation implementation-figures implementation-check

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

implementation:
	$(PYTHON) python/nisq_estimator_study.py
	$(PYTHON) python/star_channel_study.py
	$(QEC_PYTHON) python/star_growth_probe.py
	$(CUDAQ_PYTHON) python/cudaq_hwp_stage.py
	$(PYTHON) python/locality_study.py
	$(PYTHON) python/implementation_analysis.py

implementation-figures:
	$(PYTHON) python/implementation_analysis.py

implementation-check:
	$(PYTHON) python/implementation_analysis.py --check
	$(QEC_PYTHON) python/star_growth_probe.py --check

figures: discussion delivery patch-study implementation-figures
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
	$(PYTHON) python/implementation_analysis.py --check
	$(QEC_PYTHON) python/star_growth_probe.py --check
	$(PYTHON) python/presentation_manifest.py
