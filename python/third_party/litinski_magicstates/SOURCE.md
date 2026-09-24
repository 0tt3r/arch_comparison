# Vendored: Litinski, "Magic State Distillation: Not as Costly as You Think" (Quantum 3, 205, 2019)

Source: https://github.com/litinski/magicstates, commit 8fffa48c2e45d3bb0634559fce6b7b5b45f96966
(Python/ directory: definitions.py, magic_state_factory.py, onelevel15to1.py, twolevel15to1.py,
twolevel20to4.py), Apache License 2.0 (LICENSE). One modification: a stray `breakpoint()` in
`twolevel20to4.py` (line 569) is commented out, since it aborts non-interactive runs. `python/factory_calibration.py`
replaces each module's `plog` (the per-patch, per-cycle logical error 0.1(100p)^((d+1)/2)) with
the Stim-calibrated surface-code law at run time.
