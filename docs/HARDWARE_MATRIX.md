# Architectures × hardware types: plan and model definitions

Target structure of the manuscript: five architectures (Physical, STAR, Hybrid, Surface-code FT, Pinnacle)
on four hardware types (2D grid or all-to-all connectivity; uniform or non-uniform noise), for both
workloads (QUOPS random circuits, time-growing Hubbard). The hardware types change connectivity and
noise only; the clock stays at 240 ns everywhere, so the comparison isolates connectivity. A slower
syndrome cycle is, to a good approximation, a renormalization of the one-week budget (a 1 ms cycle ≈ a
week shrunk ~4000×, i.e. about ln 4000 ≈ 8 less affordable exposure), and can be added as a sensitivity.

## Noise models

| Name | Definition | Surface-code law | Status |
|---|---|---|---|
| Uniform (SD6) | p after every gate, on every idle qubit per layer, flips p on reset/measurement | Λ = 6.30, A = 0.103 (p = 10⁻³); Λ = 44 (10⁻⁴) — `results/surface_law_sd6_v1` | done |
| Non-uniform | as above but one-qubit gates p/10 (CZ schedule) | Λ = 7.41, A = 0.077; Λ = 89 — `results/surface_law_v1` | done (earlier model) |
| (optional) p/100 | one-qubit gates p/100 | needs a calibration run | not started |

Note: one-qubit noise does not enter STAR's injection error (its faults are two-qubit ZZ-rotation faults), so
non-uniform noise helps the surface code, and hence FT, as much as STAR. A STAR-favourable setting would
need angle-dependent two-qubit error; possible later variant.

## Cells

| | Physical | STAR | Hybrid | Surface-code FT | Pinnacle |
|---|---|---|---|---|---|
| 2D grid, uniform | done (v4/v5, quops_sd6) | done | done | done | not applicable (needs bounded-range non-local connectivity) |
| 2D grid, non-uniform | done (v1–v3, quops_v1) | done | done | done | not applicable |
| All-to-all, uniform | new | new | new | new | new (ledger) |
| All-to-all, non-uniform | new | new | new | new | new (ledger, uniform-noise law: conservative) |

## All-to-all model changes (declared)

- **Physical, QUOPS:** no routing; each benchmark layer runs as written (w/2 CNOTs, one layer deep).
- **Physical, Hubbard:** Jordan–Wigner still needs O(N·L) two-qubit gates per step without a local encoding
  (each vertical hop's Z string must be computed), so the gate count G is unchanged; the fSWAP network's
  O(L) depth becomes O(log L) parity trees. Exposure is unchanged; attempt time shrinks.
- **Surface code (STAR, Hybrid, FT):** logical CNOTs between patches become transversal (one code clock of
  d rounds instead of two for lattice surgery); no routing corridors, so lane and port tile charges drop to
  zero and ports are not limited by the block layout; T-state transfer from factories is a transversal
  teleportation (d rounds instead of 2·max(d, 15)). The surface-code memory law is unchanged.
- **Pinnacle:** ledger from Webster et al., arXiv:2602.11457v2 (vendored in CUDA-Q Logical):
  - codes [[30,8,4]] … [[510,16,24]], processing-block qubits 140 … 1620 (Table I);
  - logical error per logical qubit per logical cycle (A/k)(p/B)^(d/2+C), A = 6.2, B = 0.0158, C = 0.47
    (logical-measurement fit, uniform circuit noise, most-likely-error decoding: optimistic decoder);
  - logical cycle d + 2 rounds; one logical measurement (including T injection) per logical cycle per processing unit;
  - one processing unit of ⌈N_logical/k⌉ blocks plus one magic engine per sampler (their Eq. 18);
  - magic engines from CUDA-Q (p = 10⁻³: 4410 qubits, 26 cycles, output 6×10⁻¹⁰; or 5430 qubits, 1.2×10⁻¹¹;
    p = 10⁻⁴: 592/1807 qubits);
  - rotations synthesized as for FT (gridsynth/HWP T counts), each T one logical cycle, plus the repeat-until-success measurements where used;
  - cross-check: replay their Fermi–Hubbard ground-state Table IV (e.g. 58k qubits at L = 16, p = 10⁻³).

## Work plan

1. Pinnacle arm (Hubbard in `optimization_round.py`, QUOPS in `quops_study.py`), with the Table IV replay.
2. All-to-all switches for the existing arms (Hubbard and QUOPS).
3. Runs: all-to-all × {uniform, non-uniform} for both workloads.
4. Manuscript restructure: architecture sections, then hardware sections; small-multiple capability figures.
5. Optional: p/100 calibration; STAR-favourable angle-dependent noise; Pinnacle recalibration under our
   noise models (GB-code syndrome circuits + BP-OSD in Stim; an ASTRA brief).
