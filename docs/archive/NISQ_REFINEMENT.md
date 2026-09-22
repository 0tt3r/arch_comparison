# Native compilation, locality, and injected-state error

This second refinement replaces the physical routing proxy with fused native
circuits, restricts the workload and PEC noise locations to a sufficient local
observable cone, and computes spin-parity postselection acceptance under the
primitive noise model. It also separates distillation input error from physical
gate error. The derivations and figures are in [nisq_refinement.pdf](../../tex/archive/nisq_refinement.pdf),
with [LaTeX source](../../tex/archive/nisq_refinement.tex) and dedicated
[refinement_references.bib](../../tex/archive/refinement_references.bib), alongside the original
[references.bib](../../tex/archive/references.bib).

The task remains **1% relative error at provisional D*=0.1**, hence epsilon=0.001,
95% statistical confidence, U/J=4, tau=1, and r=4 unless specified. Neither the
signal scale nor the reserved Trotter bias is calibrated by these calculations.
All resource figures are conditional on those choices. No Hubbard dynamics
curve has been computed.

## What changes the conclusions

**Compilation changes the NISQ exponent substantially.** At L=4, the plaquette
circuit uses 3,149 native RZZ gates after fermionic routing, pair fusion, and
exact observable-cone restriction, versus 11,872 CNOTs in the first JW compilation.
The full fused circuit has 3,546 entanglers. The local-PEC exposure is

    log(Gamma^2) = 12.597 (2q) + 4.895 (1q) + 17.669 (idle) + 0.136 (SPAM)
                = 35.297.

With 22.07 microseconds per attempt, the week-cutoff area is 2.07e13 physical
qubits at nominal idle noise, 2.56e6 at one tenth of that rate, and 4.39e5 with
idle noise removed. These are sufficient-estimator costs, not NISQ lower bounds.
The nominal rate is p=0.001 per 100 ns, with p/10 per one-qubit gate. Both arbitrary
RZZ and one-qubit pulses take 10 ns. The old reference used native 10 ns CNOTs;
the new native gate assumption is explicit rather than a hardware calibration.

Moving preparation to first use and measurement to last use removes only 16
relevant idle qubit-layers from this L=4 circuit. Its idle exposure is almost
entirely internal. This boundary sensitivity assumes independent noise and no
measurement crosstalk; it is saved separately from the main curves.

**Exact circuit locality gives an explicit size plateau.** The sufficient
plaquette cone saturates at 36, 100, and 324 sites for r=1, 2, and 4 in the
interior cases checked. At r=1, L=8,12,16,20 all compile to the same 2,001
entanglers on 72 modes, with the same conditional resource cost. At r=4 the
L=4 cone still contains every mode, although late blocks can be omitted.
This restriction is available to all architectures. It is an exact statement
about the chosen circuit, without asserting that r=1 has sufficient Trotter
accuracy. The ordered termwise sweep has a much larger conservative circuit
dependency cone; that is not a Hamiltonian light-cone law.

**A 25% postselection yield is not evidence of an accurate observable.** We
propagate the three spin-parity characters through every native primitive and
noise location. For this model they remain Paulis, giving exact acceptance
without a many-body state simulation. At L=4 the refined plaquette acceptance
is 0.250505. The sufficient bias bound 1-P(no faults)/acceptance is 0.999726;
it does not certify epsilon=0.001. This is a loose upper bound, not a measured
bias or a lower bound. Exact spin-resolved particle-number selection and its
residual doublon bias remain open. No sampling advantage from postselection is
credited in the resource curves.

**The L=4 hybrid advantage reverses under conventional injection in the current
factory model.** Akahoshi et al. use p_in=46p/15 in their FT comparison; then
35 p_in^3 = 1.0094e-6, about 28.8 times 35p^3. The 576 fixed-angle states exceed
the strict magic budget, requiring two-level distillation. With the retained
12-tile one-level / 37-tile serial two-level factory reference:

| L=4, epsilon=0.001 | p_in=p | p_in=46p/15 |
|---|---:|---:|
| All-angle STAR qubits | 166,698 | 166,698 |
| Hybrid qubits | 66,150 | 669,714 |
| Hybrid distillation levels | 1 | 2 |
| FT qubits | 13,791,250 | 13,791,250 |

Thus the earlier 2.52-fold hybrid gain becomes a 4.02-fold penalty at this point.
The corrected hybrid still gains about 7.3-fold over all-angle STAR at L=8.
Relaxing replica integrality for the selected L=4 configurations gives about
129,701 qubits for STAR and 525,216 for the corrected hybrid, so rounding alone
does not explain the reversal. These encoded comparisons use the unchanged v1
workloads; we have not applied NISQ cone pruning to just one encoded arm.

At epsilon=0.01, one-level FT passes for p_in=p but still fails for p_in=46p/15;
the hybrid passes with one level for either input. This sensitivity keeps the
strict-precision FT words fixed to isolate sampling and factory switching.
Resynthesizing at looser precision could reduce the cost further.

The source injection coefficients are protocol- and noise-model-dependent;
they have not been recalibrated to our p_1q=p/10 convention. We retain p_in=p,
34p/15, 46p/15, and the improved prepared-state 2p/15 as explicit sensitivities.
The last is conditional on an injection supply whose acceptance and timing are
not yet compiled. Both FT and hybrid must be allowed improved injection in an
optimized comparison. These results identify a factory-switching mechanism,
not a universal hybrid penalty.

## Figures and saved evidence

- [Native compilation](../../results/nisq_v2/native_compilation.png): gate counts and
  durations for fused cones versus the v1 JW route. Full fused circuits and grid
  parity trees are also retained in the data table.
- [NISQ resource sensitivity](../../results/nisq_v2/nisq_resources.png): idle dependence
  and the L=4 exposure decomposition.
- [Locality and postselection](../../results/nisq_v2/locality_postselection.png): exact
  circuit-cone saturation, fixed-r compiled counts, parity acceptance, and bias bounds.
- [Factory reassessment](../../results/nisq_v2/factory_reassessment.png): strict/loose
  precision and injected-state quality together.

Each figure also has a vector PDF. [resource_summary.csv](../../results/nisq_v2/resource_summary.csv)
contains all 29 native cases; [factory_sensitivity.csv](../../results/nisq_v2/factory_sensitivity.csv)
contains 72 encoded sensitivity rows. Per-case directories contain the complete
native gate list, actual CUDA-Q P0 MLIR, counts, grid layout, live mode permutation,
selected lookahead and alternatives, noise exposure, and acceptance. The
[manifest](../../results/nisq_v2/manifest.json) hashes artifacts and sources, including
source snapshots matching the per-case hashes.

## Implementation and checks

Fermionic swaps occur only between adjacent modes in the live JW order. Adjacent
hopping and swap operations fuse into a two-mode Gaussian gate, implemented
with at most two RZZ gates and explicit local pulses. A separate grid-tree route
uses the square's extra edges while preserving the full JW string. This compiler
is a heuristic: the termwise fused route can have more gates than the original
JW route at larger sizes, and the plots retain that comparison.

The Phasecraft/Quantinuum reference has all-to-all connectivity, allowing ordinary
SWAP to be a free relabeling and a fermionic swap to cost one entangler. That
count cannot be transferred directly to this nearest-neighbor array; see
Sec. C4 of refs/2510.26300v2.pdf. Our fusion identities preserve the fermionic
signs rather than treating arbitrary geometric swaps as free.

There are **56 passing tests** in the qb_gsee environment. New independent
matrix/channel checks cover general two-mode Gaussian identities, four-mode
fermionic permutations, grid parity trees, occupied JW spectators under cone
compression, and parity acceptance with native noise and boundary scheduling.
These are small circuit/channel checks, not exact Hubbard dynamics or Trotter
error calibration. All 29 circuits compile with CPU CUDA-Q Logical; the saved
P0 action counts agree with the emitted native operations.

Reproduce from the repository root:

```bash
.venv-cudaq/bin/python run_nisq_refinement.py --sizes 2 4 6 8
.venv-cudaq/bin/python run_nisq_refinement.py --sizes 4 8 12 16 20 --steps 1 --algorithms plaquette --variants fused_cone
.venv-cudaq/bin/python run_nisq_refinement.py --cones-only
/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python boundary_timing.py
/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python factory_sensitivity.py
/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python plot_nisq_refinement.py
/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python -m pytest -q --junitxml=results/nisq_v2/tests.xml
tectonic nisq_refinement.tex
.venv-cudaq/bin/python run_nisq_refinement.py --manifest-only
```

The original analytical artifacts and first compilation results remain in their
versioned result directories. `compiled_model.encoded` now accepts an optional
injection error; its default remains p for reproducing the v1 model. Old stored
source hashes describe those earlier versions.

## Next architectural questions

The next NISQ refinements are better topology-aware circuit ordering and a
particle-number-conditioned observable estimator with a bias calculation.
STAR still needs explicit injection/growth/RUS schedules and a demonstration
that its residual logical channel is Pauli and branch-independent before C9's
unamplified bound becomes the default. Factory supply and consumption must then
be optimized together, including improved injection for both hybrid and FT.
The rejected encoded backend timing remains excluded from all curves; its
reduced reproducer is still pending. Trotter calibration and qubitization remain
later work in the agreed sequence.
