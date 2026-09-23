# Matched interacting reference for the Willow experiment

The [results](../results/willow_reference_v1/summary.md) separate circuit
approximation from observed experimental discrepancy for the main 4×4,
U=4, zero-flux center-hole case. Every saved nonzero time and all three
Trotter step counts are selected before inspecting their accuracy.
The 6×6 interacting experiment is outside this calculation.

## Same input and same circuit

The product occupations are read from the validated zero-time references
in the [preceding reduction](WILLOW_BENCHMARK.md). They have seven up spins,
eight down spins and a hole at (1,1). The solver uses the same open lattice,
J=1, U=4 and snake convention. Its centered onsite potential differs from
U times doublon number by a scalar within this fixed-particle sector.
That global phase affects neither observables nor the comparison.

The hopping order follows Appendix D2 of
[the pinned paper](../refs/2510.26845v1.pdf). A live fermionic permutation
tracks the alternating column swaps and applies each original nearest-neighbor
hop when its modes become adjacent. The reverse half undoes that permutation.
In the fixed fermionic frame, the circuit is the ordered hopping half sweep,
all onsite rotations, and the reverse hopping half sweep. Gate merging across
steps preserves that unitary. We simulate this ideal unitary; noise, pulse
details and physical routing costs are not inferred from it.

The noninteracting circuit reference is independently calculable with orbital
matrices. We compare all 288 saved points across both sizes, fluxes and the
three step counts. The archived ideal references include diagonal symmetry
averaging. Comparing unsymmetrized circuits to them would introduce a large
false discrepancy. The same symmetry averaging is applied to interacting
ideal observables before comparison with processed experimental estimates.
The zero-flux Hamiltonian and central-hole input are invariant under this
reflection, so averaging leaves the continuous target unchanged. A finite
ordered Trotter circuit can break the reflection symmetry.

## Numerical method

`cpp/willow_sector.cpp` reuses the frozen, previously validated matrix-free
Hubbard kernels without editing them. It constructs the (7,8) sector, which
has C(16,7)C(16,8)=147,232,800 complex amplitudes, and supports arbitrary
product occupations for small-system controls. Production is restricted to
U=4 and zero flux.

Continuous dynamics uses Chebyshev propagation in 0.05 increments through
Jt=1.2. The conservative spectral bound is 64. Each propagation has an
analytical Bessel-tail bound below 1e-13; the cumulative bound is recorded.
This is a truncation bound, not a floating-point certificate. Independent
full-Fock-space small-system tests compare complete complex states, including
fermionic signs and global-phase convention. A separate direct evolution
to Jt=1 at tighter truncation tolerance checks accumulated numerical error.

Ideal circuits are evaluated directly, using exact two-mode hopping and
onsite phase updates. Every output includes all local doublons, both spin
densities, norm, source/binary/input hashes and the PBS job identity. One
producer owns each series. Completed records are atomic and immutable;
an interrupted continuous series recomputes its prefix to recover the state,
checks the saved prefix and only writes missing time points.

## What is measured

For each observable and saved time, let E be the published processed value,
C the ideal circuit value and H the continuous Hamiltonian value. We retain

    E − H = (E − C) + (C − H).

The first term contains experimental sampling, mitigation and other processing
effects. It is an observed discrepancy, not the expectation-value bias of an
estimator. The second term is the deterministic circuit approximation error.
Their RMS values do not add: the squared total has a cross term that can be
negative. The analysis saves that term and checks the identity.

Finite-bond MPS references are compared to H separately. Published marginal
uncertainties remain statistical quantities; they are not replaced by these
observed differences. Temporal and observable covariances of the processed
estimates remain unavailable. Consequently, fractions of saved points within
a tolerance describe this data set, not a coverage guarantee.

The architecture curves continue to use their declared task and scenario
assumptions. This experiment has a central hole and a different hopping order
from the baseline Néel/plaquette problem; its error coefficient is not silently
substituted into that baseline.

## Reproduction

Offline analysis and figures use the committed sector results:

```bash
make willow-reference
make willow-reference-check
```

To rerun the expensive calculation, use the declared scientific environment,
compile with `python python/willow_sector.py --prepare`, and run the smoke
and production tasks in `pbs/willow_reference.pbs`. The wrapper verifies the
source and binary contract before executing. Validate the actual PBS smoke
before submitting production. The current run uses one whole Improv node,
16 OpenMP threads and NUMA interleaving; it does not change existing jobs.
Do not rebuild the production binary or change its pinned sources while a
writer is active. Tests compile a separate temporary binary.
