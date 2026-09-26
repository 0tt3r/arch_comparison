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
| All-to-all, uniform | done (v6, quops_a2a_sd6) | done | done | done | done (ledger) |
| All-to-all, non-uniform | done (v7, quops_a2a) | done | done | done | done (ledger, uniform-noise law: conservative) |

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

## Results (p = 10⁻³; figures/hardware_matrix_hubbard.pdf, figures/hardware_matrix_quops.pdf)

Hubbard, τ = L/4, ε = 0.01, full source menu (calibrated_factories):

| Hardware | Physical (eff.) | STAR | Hybrid | FT: L=16 / τ=8 | Pinnacle: L=16 / τ=8 |
|---|---|---|---|---|---|
| Grid, uniform (v5) | L=6 | L=10 | L=12 | 1.9e6 / 9.6e6 | n/a |
| Grid, non-uniform (v3) | L=6 | L=10 | L=14 | 1.5e6 / 7.8e6 | n/a |
| All-to-all, uniform (v6) | L=6 | L=10 | L=12 | 8.9e5 / 4.7e6 | 5.8e4 / 2.2e5 |
| All-to-all, non-uniform (v7) | L=6 | L=10 | L=14 | 7.2e5 / 3.8e6 | 5.8e4 / 2.2e5 |

At ε = 0.001 Pinnacle's lead narrows (p = 10⁻³: τ=8 at 2.9e7 vs FT 4.6e7 on all-to-all uniform) and reverses
at p = 10⁻⁴ (9.7e6 vs 4.7e6): one sampler has serial T gates (1.1e7 at L=32), so its cost is copies.

QUOPS, standard threshold (score at 10⁴ / 10⁶ / 10⁹ qubits; FT includes factories):

| Hardware | Physical | STAR (plateau) | FT | Pinnacle |
|---|---|---|---|---|
| Grid, uniform | 240 | 6.2e3 | 1.4e2 / 4.8e6 / 5.0e8 | n/a |
| Grid, non-uniform | 310 | 7.1e3 | 2.3e2 / 5.5e6 / 7.4e8 | n/a |
| All-to-all, uniform | 660 | 6.2e3 | 2.7e3 / 9.4e6 / 5.5e8 | 2.9e5 / 8.2e6 / 1.2e9 |
| All-to-all, non-uniform | 1.7e3 | 7.1e3 | 3.4e3 / 9.7e6 / 8.1e8 | 2.9e5 / 8.2e6 / 1.2e9 |

Pinnacle's jump near 10⁴ qubits is one sampler of six [[254,14,16]] blocks plus the 4410-qubit engine
(9570 qubits), allowed because the fit gives 3e-11 per logical cycle at d = 16; not a bug, but a ledger
consequence. The grid QUOPS panels are computed from costs.json on the 10–10⁹ grid (quops_v1's
capability.json starts at 10²).

## CUDA-Q Logical cross-check (results/cudaq_pinnacle_check_v1, `make pinnacle-check`)

- Codes (n, k, d, block qubits, cycle rounds) and all four engines (qubits, cycles per attempt, acceptance,
  output error) match python/pinnacle.py entry by entry.
- A 2×2 Hubbard plaquette step on a gb30 + p1e-3 engine device compiles to P2 via Clifford+T and PBC:
  744 T → 744 engine requests → 744 `inject` gadgets (one logical Pauli-product measurement each). This
  confirms the ledger's one-T-per-logical-measurement rule.
- CUDA-Q binds one extra processing block per sampler as RUS scratch (4690 vs 4550 qubits here); the
  ledger (paper Eq. 18) omits it. In our chosen designs: ≤3% for Hubbard L ≥ 16, +9% for the 9570-qubit
  QUOPS sampler, up to +43% for one-block samplers (gb126 + 592-qubit engine, QUOPS w = 12 at p = 1e-4).
- Not checkable in cudaq-logical 0.1.1.post1: P3 lowering of the GB codes ('fabric.resource_rotate_product'
  op has an invalid resource-assisted physical Pauli product), the analytical tier (WSC output_syndrome
  ops have no cost semantics), Pinnacle's RUS rotation route ('fabric.retry' op retry success semantics do
  not match an exact ordered OutcomeMap row), and the logical-error fit (CUDA-Q uses a generic
  0.1(p/0.01)^((d+1)/2) law). RUS would lower T counts, so our gridsynth counts are conservative for Pinnacle.

## Work plan

1. (done) Pinnacle arm (Hubbard in `optimization_round.py`, QUOPS in `quops_study.py`), with the Table IV replay.
2. (done) All-to-all switches for the existing arms (Hubbard and QUOPS).
3. (done) Runs: all-to-all × {uniform, non-uniform} for both workloads.
4. (done) Manuscript restructure: architecture sections, then hardware sections; small-multiple capability figures.
5. Optional: p/100 calibration; STAR-favourable angle-dependent noise; Pinnacle recalibration under our
   noise models (GB-code syndrome circuits + BP-OSD in Stim; an ASTRA brief).
