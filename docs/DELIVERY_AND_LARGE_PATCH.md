# Delivery tradeoffs and the 4 × 4 calibration

This pass tests two assumptions that can change the comparison: encoded
consumption throughput and transfer of the Trotter coefficient from small
patches. The original discussion-v3 curves remain the common reference.

## What is being changed in the delivery comparison

The [capacity figure](../figures/delivery_capacity.pdf) and
[shot-time breakdown](../figures/delivery_bottlenecks.pdf) retain the same
noise, precision, and finite-patch error envelope. Error allocations, distances,
source cleanup, source count, and integer replicas are reoptimized for each
layout. Both FT and hybrid get the same source menu, extended to 1024 sources
so a faster consumer is not artificially starved at the old 16-source cap.
Selected solutions must lie below the new source-count endpoint.

The previous convention was already **four effective consumption lanes** at
nine code clocks per T *per lane*. It is not the paper's serial compact block.
The paper's compact, intermediate, and fast blocks use respectively

| Block | Data tiles for n logical qubits | Code clocks per consumed T, whole block |
|---|---:|---:|
| Compact | ceil(1.5n + 3) | 9 |
| Intermediate | 2n + 4 | 5 |
| Fast | ceil(2n + sqrt(8n) + 1) | 1 |

These are the formulas in [Litinski, Quantum 3, 128 (2019), Sec. 2](https://quantum-journal.org/papers/q-2019-03-05-128/).
Here n = 2L² + 6, including the existing work-qubit allowance, and six
additional tiles per delivery/injection lane remain charged. The rounded
fast-block formula is an area estimate, not a placed rectangular floorplan.
In particular, the relative area increase at small L is larger than its
asymptotic 4/3 ratio.

The whole-block T rate is **not divided again by the number of analog lanes**.
An explicitly separate four-fast-port scenario performs that division and
charges the larger fast-block area plus four lane allowances. It is an
uncompiled throughput sensitivity; four independent arbitrary Pauli-product
measurements may require more area or a different dependency schedule.

STAR's rotation cost stays at 18 code clocks; it does not inherit the T-state
consumption speedup. We vary its effective analog-lane count over 1, 2, 4, 8,
and 16, with area charged in every case. The hybrid varies its analog lanes
and T delivery independently. Multi-lane rates remain assumptions; the
previous JW dependency calculation does not certify these higher rates.

The two overlaid menus are nested: the reference remains available, then
whole-block layouts and analog-lane choices are added, then four fast T ports.
Thus a larger menu cannot worsen cost. This is not a claim that every point
has a compiled implementation. Sparse sampled sizes are reported directly;
there is no interpolation to untested lattice sizes.

For source count f, T count N_T, analog count R, and distance d, the code uses

    t_supply = ceil(N_T/f) service(d) t_cycle
    t_consume = d [c_T N_T / b_T + (18 R + Clifford clocks)/b_A] t_cycle
    t_shot = max(t_supply, t_consume) + initialization + feedforward.

Here b_T is explicitly one for a paper whole-block rate. Both supplied-state
quality and logical failure volume are included in the error optimization.
The total-qubit objective includes integer replicas, so a faster shot need
not reduce total area when only one copy is needed.

At strict 4 × 4 and tau = 1, the single fast block lowers FT's reference area
from 1.43 million to 0.86 million qubits. The four-port scenario reaches
0.37 million and becomes supply-limited. STAR and hybrid prefer two lanes,
with about 61,000 and 63,000 qubits. At loose precision, one copy often fits
the week, and larger/faster blocks can lose to smaller ones. These comparisons
still use the 3 × 3 planning envelope; changing layout and changing physical
calibration are kept separate in the saved evidence.

The [complete menu](../results/delivery_v1/costs.json),
[selected budgets and bottlenecks](../results/delivery_v1/summary.csv), and
[assumptions](../configs/delivery_sensitivity.json) are reproducible with
`make delivery`.

## Matrix-free exact dynamics

The 4 × 4 half-filled fixed-spin sector has choose(16,8)² = 165,636,900
amplitudes. One complex128 vector uses 2.65 GB. We apply

    H Psi = h_up Psi + Psi h_down^T + V * Psi

without constructing the full Hamiltonian. Fixed-weight bit masks carry
fermionic signs, and the Neel state has eight particles of each spin. Exact
reference evolution uses a Chebyshev expansion with ||H|| <= 64 in hopping
units. For a = 64 tau,

    exp(-i H tau) = J_0(a) + 2 sum_{k>=1} (-i)^k J_k(a) T_k(H/64).

The retained degree is chosen using

    |J_k(a)| <= (a/2)^k / k! exp[a²/(4(k+1))]

and a geometric bound on the remaining terms, giving a state-norm truncation
bound below 1e-12. This bound excludes floating-point roundoff; norm checks
and independent full-state small-system comparisons test the implementation.
The local [0,1] observable error due to the polynomial truncation is at most
2 eta + eta² for state error eta. No dense diagonalization is required.

Plaquette stages use exact fermionic Hadamard transformations within disjoint
2 × 2 tiles; clipped two-site tiles use exact hopping gates. Termwise stages
apply exact two-state hopping rotations in precisely the saved circuit order.
Both orderings, and the Chebyshev reference, match independent Python full
state vectors on 2 × 2, 3 × 2, and 3 × 3 to better than 2e-11. All-site
observables are saved so the observed central coordinate can be checked
without rerunning the calculation or assuming a Trotter reflection symmetry.

The full-size topology smoke tested 8, 16, 32, and 64 OpenMP threads; extra
threads did not improve its short matrix-vector workload. Production uses 16.
A subsequent full-size Trotter benchmark found a fourfold speedup from
interleaving memory over all eight NUMA domains. The initial run was stopped
after its saved, converged r = 32 plaquette case; the remaining cases resumed
with interleaved placement. This intentional stop is preserved in the job
accounting, separately from successful completion of the resumed calculation.
The project allocation was exhausted, so the bounded jobs explicitly use
Improv's allocation-exempt backfill queue. Inputs, source and executable
hashes, per-case host/job IDs, CPU time, and cumulative child peak RSS accompany atomic
completion records. An interrupted case restarts from its initial state;
completed cases are reused only when their entire calculation contract matches.
The PBS scripts verify frozen input hashes before running. No state vectors
or executables are committed.

This initial calculation calibrates one finite 4 × 4 lattice at tau = 1.
The [follow-up](DIRECT_CHECKS_AND_TIME_EXTENSION.md) checks the selected
headline step counts and extends time to tau = 2,4. Neither calculation
establishes a bulk coefficient or validates a 4 × 4 transfer at tau = 8.

## Completed finite-lattice result

The [4 × 4 summary](../results/dynamics_4x4/summary.md) and
[figure](../figures/large_patch_calibration.pdf) now contain both orderings at
r = 8,16,32. The exact central doublon is 0.1452313819. The two finest
coefficients agree within 0.23% for plaquette and 0.35% for termwise. Their
r = 32 coefficients are 0.26183 and 0.18773, respectively: about 20% and
145% above the 3 × 3 values. The earlier 10% transfer margin does not cover
either. The new 4 × 3 tau = 1 check falls between the two square patches.

The interpolated 4 × 4 sensitivity uses coefficients 0.291324 and 0.210136,
including the 10% planning margin. At strict precision the reference STAR
and hybrid areas remain about 71,000 and 73,000 qubits, despite needing more
steps; one replica still fits the week. FT becomes about 1.88 million with
the reference delivery, 1.16 million with one fast block, and 0.385 million
with the uncompiled four-port assumption. Effective NISQ changes from 36
to 72 qubits, with mitigation bias still assumed affordable. These are
interpolated costs at one finite-lattice point, not a replacement capacity
curve or an accuracy validation of the NISQ estimator.

The directly measured step menu is also retained. It gives a strict NISQ
effective cost of 2,196 qubits because its available r = 32 overshoots the
interpolated optimum near 21. This is a step-menu penalty, not evidence of
a worse architecture. Reporting that table alone would introduce a new
comparison artifact. Exact r = 8,16,32 costs and interpolated optima are
therefore labeled separately.

The follow-up directly checks every selected headline r before reallocating
any error budget. All 16 saved allocations pass. Using their measured errors
then lowers strict reference FT to 1.66 million qubits and the single fast
block to 0.99 million; reference STAR, hybrid, and effective NISQ keep the same
areas. The earlier interpolation table remains frozen as the input to that
test. See the [new summary](../results/dynamics_followup_v1/summary.md).
