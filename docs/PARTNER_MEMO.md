# Notes for the partner efforts (crossover and IQUS)

24 September 2026. Every item below was re-checked against the current code in this repository
(`crossover_martin/`, and `fh_iqus_idears/` at commit `222131e`); line numbers refer to those copies.
The three efforts ask somewhat different questions, and most differences between them are choices,
not errors. The items marked **error** change a result in a direction the authors would likely not
intend; the rest are offered as context.

## crossover (Martin)

1. **Factory table (error).** `crossover/architectures.py:114–119`, `DEFAULT_FACTORIES`, sets
   `p_ref = 1e-3` for every row, but rows 1–4 are the p = 10⁻⁴ designs of Litinski (2019), Table 2
   (e.g. (15-to-1)₇,₃,₃: 4.4×10⁻⁸ at p = 10⁻⁴). Row 5 combines three designs: the name of the
   p = 10⁻³ small-footprint (15-to-1)₉,₅,₅×(15-to-1)₂₁,₉,₁₁ (7,780 qubits, 469 cycles, 6.1×10⁻¹⁰),
   the qubits and cycles of the p = 10⁻³ (15-to-1)⁶₁₁,₅,₅×(15-to-1)₂₅,₁₁,₁₁ (30,700, 82.5), and the
   output error of the p = 10⁻⁴ (15-to-1)⁴₉,₃,₃×(15-to-1)₂₅,₉,₉ (6.3×10⁻²⁵). At p = 10⁻³ the table is
   therefore optimistic by one to many orders of magnitude in output error; `error_at` then rescales
   from the wrong reference. The verified Tables 1–2 are in `configs/optimization_round_v3.json`
   (`litinski_factories`).
2. **FT logical error omits factory waiting (error).** `architectures.py:170–173` charges logical
   faults over `cycles_logic`, but the shot lasts `max(cycles_logic, cycles_fac)` (`:193`). When
   factories are the bottleneck, the data patches idle unprotected in the ledger.
3. **STAR per-rotation error.** `architectures.py:219, 233`: 0.2p per rotation, independent of angle,
   and T/Toffoli angles counted as analog rotations (`:231`). The transversal multi-rotation
   instrument of Toshio et al. gives an error that grows with the angle; for small Trotter angles
   0.2p is pessimistic, and for π/8 rotations it is strongly optimistic. In our runs, injecting a
   plaquette circuit's fixed π/8 rotations in analog costs 10¹⁷·⁴ qubits at L = 12 (p = 10⁻³,
   ε = 0.01), against 10¹⁰·⁴ for the termwise circuit.
4. **Trotter law.** `benchmarks/hubbard.py:5` uses an extensive bound, α·N·t³/r². For a local
   observable the error does not grow with N once the lattice contains its light cone; exact strips
   up to 8×2 and 6×3 confirm this for the U/J = 4 doublon at τ = 1 ([STRIP_DYNAMICS.md](STRIP_DYNAMICS.md)).
5. **QUOPS physical arm (about 150 against our 310).** `architectures.py:92` includes readout error
   (1% per measured qubit) and T₂ idling in the circuit fidelity. QUOPS scores process polarization,
   from which the mirror test removes state preparation and measurement; at w ≈ 17 readout alone
   uses about a third of the e^{-1/2} budget. This likely accounts for most of the factor of two.
6. **Physical-qubit mitigation model.** F = exp(−ΣpG) with variance 1/F² (`:92`) is a global
   depolarizing model; two-qubit PEC has exponent about 4pG for a depolarizing gate, and a local
   observable's attenuation is not the full-circuit fidelity.

## IQUS

1. **STAR logical-error criterion lacks the factor of two (error).** `fhcost/ftqc.py:446` tests
   q_L·rounds·p_L(d) ≤ ε_L, while the FT arm (`:387`) and the file's own convention (`:693`) double
   it because a flipped ±1 outcome biases the estimator by twice the fault probability.
2. **STAR timing and area.** `ftqc.py:445` charges 20 rounds per Trotter step whatever the number of
   rotations, and `:452` adds no injection area. Our delivery sensitivity finds STAR limited by
   exposure rather than delivery, so this mainly affects shot time, but a per-rotation service time
   with lanes would make STAR's runtime scale with the circuit.
3. **Default STAR law.** `budget.py:160, 189`: the default is the count law Γ² = exp(1.066 p N_rot),
   from the 2024 protocol (P = 2p/15 per rotation). The angle law of Toshio et al. (already
   implemented as `star_law = "angle"`) is the current reference and is much cheaper for small
   Trotter angles.
4. **Comparison preset.** `fhcost/presets.py:18–20` maps our model to ε_abs = 10⁻³ (1% of D* = 0.1).
   That is our strict case; our main case is ε = 0.01 absolute. Both are now reported.
5. **Surface-code law.** Both efforts use p_L = 0.1(p/0.01)^((d+1)/2). Under a circuit-level
   calibration with Stim (Λ = 7.4, A = 0.077 at p = 10⁻³) the distance rises by 4–6 at the operating
   points of this study.

## Results the partners may find useful

- **Calibrated surface-code law**: `results/surface_law_v1/fit.json` ([SURFACE_LAW.md](SURFACE_LAW.md)).
- **Lattice surgery**: isolated CNOTs fail at 0.73–1.15× the memory law at d merge rounds and take
  2d clocks plus 300 ns ([SURGERY_LAW.md](SURGERY_LAW.md)).
- **Factories re-priced under the calibrated law**: a Pareto menu of 26 (p = 10⁻³) and 34 (10⁻⁴)
  designs in `results/factory_calibration_v1/menu.json`, from Litinski's own code
  ([FACTORY_CALIBRATION.md](FACTORY_CALIBRATION.md)).
- **Trotter coefficients**: exact U/J = 4 doublon coefficients on strips to 8×2 and 6×3
  ([STRIP_DYNAMICS.md](STRIP_DYNAMICS.md)); the termwise value converges to 0.19 at τ = 1.
- **STAR angle law**: saturation s ≈ 0.41/(p⟨|θ|⟩) on random circuits
  ([QUOPS_SENSITIVITY.md](QUOPS_SENSITIVITY.md)); raw RUS versus PCEC reconciled with Toshio et al.
  ([RAW_RUS_NOTE.md](RAW_RUS_NOTE.md)).
