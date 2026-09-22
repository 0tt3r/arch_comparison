# Direct step checks and longer-time calibration

This calculation addresses two different questions. Does the finite 4 × 4
Trotter circuit satisfy the error allowance used in each headline resource
estimate? How much does the longer-time comparison change when its accuracy
model includes larger finite patches?

The [result summary](../results/dynamics_followup_v1/summary.md) records the
answers. The [capacity figure](../figures/patch_capacity.pdf) puts available
physical qubits on the horizontal axis. It isolates the physics calibration
from the separate [delivery comparison](DELIVERY_AND_LARGE_PATCH.md).

All 16 saved allocations pass. Reallocating their measured Trotter error lowers
strict reference FT from 1.88 to 1.66 million qubits and the single fast block
from 1.16 to 0.99 million; reference STAR, hybrid, and effective NISQ retain
their qubit counts. The final tau = 2,4 coefficients converge within 0.36%
between the finest pairs. At tau = 4 the plaquette coefficient is 4.193,
close to 4 × 3 and below 3 × 3, while termwise is 3.029, above both smaller
patches. A universal multiplicative patch correction would miss this behavior.

## Check the saved allocation first

The frozen input is `results/dynamics_4x4/interpolated.json`. We select the
effective and two PEC NISQ scenarios, reference STAR and hybrid, and reference,
single-fast-block, and four-fast-port FT. At both precision targets these give
16 resource estimates sharing nine distinct circuits. One circuit was already
measured; eight new circuits check the remaining selected step counts.

For each row, the audit compares

    |D_r - D_exact| + 1e-10 <= epsilon * saved_trotter_fraction.

It does this **before changing any budget fraction**. The saved qubit estimate
can therefore be assessed without silently moving error allowance from another
component. The precision convention stays fixed: epsilon = 0.001 means 1% of
the reference signal D* = 0.1, with epsilon = 0.01 as the loose sensitivity.
The newly calculated signal does not redefine either the target or shot model.

The next two tables answer a different optimization question. They charge the
measured Trotter error plus the numerical allowance and reoptimize the other
budget fractions, first at the saved r, then over all measured r values.
The latter remains a discrete candidate menu, not a global optimum. In
particular, it removes the earlier artificial restriction to r = 8,16,32.

These checks validate the deterministic Trotter component on this finite
lattice. They do not measure the residual mitigation bias. The NISQ discussion
model still transfers plaquette accuracy to an experimentally anchored gate
count proxy; the four-fast-port FT design remains a throughput assumption.

## Extend time with the same finite-lattice solver

The new runs use the validated matrix-free fixed-spin solver without changing
its Hamiltonian, initial state, ordering, or binary. Exact Chebyshev references
and both product formulas are evaluated at tau = 2 and 4. The observed site
remains (2,2). All-site doublons, norms, source/binary hashes, and execution
records are retained in a new result directory.

At each time, the initial step menu is 8 tau, 16 tau, and 32 tau. The driver
doubles the finest count if the last two absolute coefficients differ by more
than 5%, up to the declared cap of 128 tau. Failure at that cap stops the
analysis rather than labeling the coefficient converged. A converged pair
checks the leading r dependence on this finite patch; it does not bound an
infinite-lattice error or eliminate finite-time cancellations.

The three PBS jobs have disjoint output directories, and retain the previously
validated one-node, 16-thread, interleaved-memory topology. Existing completed
tau = 1 cases are reused only under matching contracts. The old result
directories remain unchanged. Terminal PBS accounting and a case census are
checked in addition to process exit status.

## Isolate the effect on capacity

The capacity comparison holds the reference delivery assumptions and expanded
shared source-count menu fixed. It reoptimizes allocations for two coefficient
policies:

- The earlier 3 × 3 running-maximum envelope, including its finite-step checks.
- An envelope that retains that reference and also includes every measured
  4 × 3 and 4 × 4 coefficient through the requested time, with the same 10%
  margin.

The second policy uses finite-r coefficients as well as converged tails.
Earlier peaks remain in the envelope, so a later local cancellation cannot
improve the planning cost. The L = 2 point retains its own measured calibration
in both policies. The same discrete lattice-size grid is used, with no
interpolation between sizes.

This is a **finite-patch sensitivity**, not a bulk bound. The larger patches
are sampled only at tau = 1,2,4; unsampled time peaks remain possible. Transfer
to L > 4 is still an assumption. No 4 × 4 coefficient is inferred at tau = 8,
and the original discussion-v3 curves remain available separately. Optimized
r values in this new capacity sweep are planning estimates; only the specified
tau = 1 headline allocations have the direct checks described above.

Run `make patch-study` to regenerate this analysis and its figures from the
saved dynamics. `python python/dynamics_followup_analysis.py --check` audits
the frozen allocations, raw case hashes, coefficients, budget/resource replay,
source-count endpoint, capacity inversion, and terminal job accounting.
