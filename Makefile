PYTHON ?= python3
CUDAQ_PYTHON ?= .venv-cudaq/bin/python
QEC_PYTHON ?= /home/ac.otten/soft/mambaforge/envs/stim/bin/python
export OPENBLAS_NUM_THREADS = 1
export OMP_NUM_THREADS = 1
.DEFAULT_GOAL := paper
.PHONY: discussion delivery patch-study physics reconcile figures paper test audit implementation implementation-figures implementation-check frontier frontier-figures frontier-check integrated integrated-figures integrated-check published-benchmarks published-figures published-check willow-benchmark willow-figures willow-check willow-reference willow-reference-figures willow-reference-check

willow-reference:
	$(PYTHON) python/willow_reference_analysis.py
	$(PYTHON) python/plot_willow_reference.py

willow-reference-figures:
	$(PYTHON) python/plot_willow_reference.py

willow-reference-check:
	$(PYTHON) python/willow_reference_analysis.py --check

willow-benchmark:
	$(PYTHON) python/willow_benchmark.py
	$(PYTHON) python/plot_willow.py

willow-figures:
	$(PYTHON) python/plot_willow.py

willow-check:
	$(PYTHON) python/willow_benchmark.py --check

published-benchmarks:
	$(PYTHON) python/published_benchmarks.py
	$(PYTHON) python/plot_published_benchmarks.py

published-figures:
	$(PYTHON) python/plot_published_benchmarks.py

published-check:
	$(PYTHON) python/published_benchmarks.py --check

integrated:
	$(PYTHON) python/integrated_comparison.py
	$(PYTHON) python/plot_integrated.py

integrated-figures:
	$(PYTHON) python/plot_integrated.py

integrated-check:
	$(PYTHON) python/integrated_comparison.py --check

frontier:
	$(PYTHON) python/qubitization_study.py
	$(PYTHON) python/geometric_delivery.py
	$(PYTHON) python/exposure_diagnostics.py
	$(PYTHON) python/locality_certificate.py
	$(PYTHON) python/frontier_analysis.py

frontier-figures:
	$(PYTHON) python/frontier_analysis.py

frontier-check:
	$(PYTHON) python/frontier_analysis.py --check
	$(QEC_PYTHON) python/frontier_growth.py --check

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

figures: discussion delivery patch-study implementation-figures frontier-figures integrated-figures published-figures willow-figures willow-reference-figures
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
	$(PYTHON) python/frontier_analysis.py --check
	$(QEC_PYTHON) python/frontier_growth.py --check
	$(PYTHON) python/integrated_comparison.py --check
	$(PYTHON) python/published_benchmarks.py --check
	$(PYTHON) python/willow_benchmark.py --check
	$(PYTHON) python/willow_reference_analysis.py --check
	$(PYTHON) python/presentation_manifest.py
