# ASTRA brief: calibrate the Pinnacle logical-error law under our noise models

## Why

The manuscript's Pinnacle arm (qLDPC, generalized-bicycle codes) wins on all-to-all hardware by up to 20×.
That lead rests on a law we did not calibrate. It comes from Webster et al., arXiv:2602.11457v2 (local copy:
`refs/2602.11457v2.pdf`), Table II, "Log. Meas." row:

    p_L per logical qubit per logical cycle = (A/k) (p/B)^(d/2 + C),   A = 6.2, B = 0.0158, C = 0.47

It was fitted with most-likely-error (optimal) decoding, under uniform circuit noise, and only up to d = 16.
We then extrapolate it to the d = 24 code. We need the same quantity measured with a practical decoder
(BP-OSD) under our two noise models, over the range that can be sampled directly. That tells us whether the
paper's law is optimistic, and by how much.

This is a self-contained task. It does not depend on anything else running.

## Codes (Table I of the paper; also `python/pinnacle.py`)

Each code is a generalized-bicycle (GB) code of length n = 2ℓ. With the cyclic shift S of size ℓ, define

    A = Σ_{i∈a} S^i,   B = Σ_{j∈b} S^j,   H_X = [A | B],   H_Z = [Bᵀ | Aᵀ]

| name | ℓ | a | b | n | k | d (paper) | processing-block qubits |
|---|---|---|---|---|---|---|---|
| gb30 | 15 | (0,6,13) | (0,1,4) | 30 | 8 | 4 | 140 |
| gb62 | 31 | (0,6,15) | (0,5,7) | 62 | 10 | 6 | 244 |
| gb126 | 63 | (0,4,37) | (0,29,49) | 126 | 12 | 10 | 452 |
| gb254 | 127 | (0,32,100) | (0,28,49) | 254 | 14 | 16 | 860 |
| gb510 | 255 | (0,39,55) | (0,70,127) | 510 | 16 | 24 | 1620 |

A machine-readable copy, with logical bases, is in CUDA-Q Logical:
`/home/ac.otten/soft/cuda-quantum/preview/logical/python/cudaq/logical/codes/pinnacle.py`
(`PINNACLE_GB_INSTANCES`, `pinnacle_gb(name)`). Check k by rank before using any code.

## Noise models (match ours exactly)

- **SD6 (uniform):** depolarizing error p after every one- and two-qubit gate; DEPOLARIZE1(p) on every idle
  qubit in every layer, including the measurement and reset layers; X/Z flip with probability p after every
  reset and before every measurement. The reference implementation is `circuit()` in
  `python/surface_law_sd6.py`, which adds this noise layer by layer to a noiseless Stim circuit. Reuse that
  pattern.
- **Non-uniform:** the same, but one-qubit gates and one-qubit idles use p/10. Two-qubit gates, reset and
  measurement stay at p.
- Physical error rates: p ∈ {1e-3, 2e-3, 3e-3, 5e-3} as the main grid, plus 1e-4 wherever it can be resolved.

## What to measure

1. **Syndrome-extraction circuit** for each code with all-to-all connectivity (any qubit pair can interact).
   Use the paper's schedule if it gives one (check its appendices). Otherwise build a valid schedule from
   an edge colouring of the weight-6 checks, with one ancilla per check, and state which one you used.
   Verify it in two ways:
   - the noiseless circuit fires no detectors;
   - the circuit distance, from `stim.Circuit.search_for_undetectable_logical_errors`, or a bound from it.
     Report the circuit distance next to the paper's d.
2. **Memory experiment:** prepare in Z (and separately X), run d rounds (or r = d + 2, the paper's logical
   cycle), then measure out all k logical observables. Report:
   - p_L per logical qubit per cycle, from p_block = 1 − (1 − p_L)^(k·cycles);
   - the per-round rate as well.
3. **Decoder:** BP-OSD from `ldpc` (version 2.4.1 in the stim env: `ldpc.BpOsdDecoder` or
   `ldpc.SinterBpOsdDecoder` with sinter), OSD-CS order ≥ 7, max_iter ≈ 1000. Build it on the detector
   error model without graphlike decomposition. Record every decoder parameter. Optionally, add a BP-LSD
   comparison.
4. **Sampling:** at least 100 logical failures per point where affordable, with Wilson 95% intervals.
   - Direct sampling is expected to resolve gb30, gb62 and gb126 over the grid, and gb254 only at p ≥ 3e-3.
   - Do not try to resolve gb510 or gb254 at p = 1e-3 directly.
5. **Comparison:** fit the paper's functional form (A/k)(p/B)^(d/2 + C) to your points, per noise model.
   Report:
   - (A, B, C) with intervals;
   - the ratio yours/paper at each measured (code, p);
   - the extrapolated ratio at gb254 and gb510 for p = 1e-3 and 1e-4, with the caveat that it is
     extrapolated.

   The paper's law is for a logical *measurement* gadget, not memory. Also find the paper's memory-type
   row in Table II, if there is one, and compare like with like. Otherwise state the mismatch plainly.

## Environment

- Python: `/home/ac.otten/soft/mambaforge/envs/stim/bin/python` (stim 1.16.dev0, pymatching 2.3.1,
  ldpc 2.4.1, sinter). Install nothing into `qb_gsee` or `.venv-cudaq`. If you need a package, install it
  into the stim env and say so.
- Cluster: LCRC Improv, PBS. The template is `pbs/surface_law_sd6.pbs`: account `QC_NA_Sim`, queue
  `backfill` or `compute`, 128 cores per node, and `OMP_NUM_THREADS=1` with one worker process per core.
  Make the runs resumable, with one JSON per (code, p, noise, basis, rounds, seed), written atomically as
  in `python/surface_law_calibration.py`.
- Seeds: fixed and recorded in a manifest.

## Files

You may create:
- `python/pinnacle_calibration.py` (circuit construction, validation, sampling, fit; subcommands
  `validate`, `run`, `fit`, like `python/surface_law_sd6.py`);
- `pbs/pinnacle_calibration.pbs`;
- `results/pinnacle_calibration_v1/` (manifest, per-cell JSON, `fit.json`, `validation.json`);
- `docs/PINNACLE_CALIBRATION.md` (method, table, fit, comparison, caveats);
- `tests/test_pinnacle_calibration.py`, with small, fast checks: code k by rank, noiseless circuit
  silent, one tiny sample.

Do **not** modify:
- `python/pinnacle.py`, `python/optimization_round.py`, `python/quops_study.py`,
  `python/surface_law_*.py`, `python/growth*.py`;
- anything in `tex/`, `figures/`, `Makefile`, or any existing folder under `results/`;
- `results/presentation_v1/manifest.json`.

These are hash-pinned or feed the manuscript. Wiring the new law into the design search is our job
afterwards.

Do not commit or push. Leave the work uncommitted for review.

## Deliverable and stopping point

Stop when `docs/PINNACLE_CALIBRATION.md` contains:
1. the circuit used and its verified circuit distance per code;
2. a table of p_L per logical qubit per logical cycle, with intervals, for each measured (code, p, noise);
3. the fitted (A, B, C) per noise model, and the yours/paper ratios, measured and extrapolated;
4. a one-paragraph verdict: is the paper's law optimistic, and by roughly what factor at gb254/gb510 and
   p = 1e-3?

If a full grid would take more than about two node-days, run the SD6 model at p ∈ {2e-3, 3e-3, 5e-3} for
gb30–gb126 first. Report that, and then ask before going further.
