# Published benchmarks on their original tasks

The [result summary](../results/published_benchmarks_v1/summary.md) and two figures
answer a narrower question than an architecture ranking: can the common resource
accounting explain the published operating points, and which ingredients have
actually been checked? The integrated capacity curves retain their inputs.

The subsequent [Willow follow-up](WILLOW_BENCHMARK.md) now reduces the
superconducting data; the availability-only statement below records this
earlier pass's scope.

## STAR and FT: reconstruct the published ledger

We transcribe Table IV of the local, hashed
[Chung et al. paper](../refs/2603.13093v1.pdf) into
`configs/published_benchmarks.json`. Its task is ground-state energy estimation
at error 0.005 L². It assumes suitable ground-state overlap and uses QCELS for
STAR and phase estimation with plaquette Trotterization/Hamming-weight phasing
for FT. This is not our local-observable dynamics task.

The independently evaluated STAR footprint is (4L²+1)2d_f². Hamiltonian term
count is 8L(L−1)+L², and the reported deepest-circuit rotation count is
approximately twice this count times the maximum number of Trotter steps.
The QCELS time grid uses n=1,...,5 at each of J levels, with 100 samples per
quadrature and both real and imaginary parts. We retain the paper's
epsilon_QCELS=2 epsilon/3, epsilon_Trotter=epsilon/3 and delta=0.06.

For time fraction f=(n/5)2^(j−J), we take ceil(r_max f) steps. Appendix D gives
the control overhead 16+18r in addition to the supplied Trotter clocks. A
preparation/growth clock takes d_i+d_f+6 syndrome rounds. A 350 ns round
(six 25 ns gate layers, measurement and reset) reproduces the published
single-shot column, except for the flagged entry. Summing both quadratures,
all nonzero times, integer steps and control operations gives the PEC-free
campaign duration. The zero-time circuit is omitted from this evolution-time
ledger, as in the analytical sum; a physical campaign also prepares and reads
out those samples.

Distances and RUS clock counts are **supplied published inputs**. The paper's
growth law contains numerical fit parameters that are not tabulated, and its
predecessor's supplementary-code endpoint returned HTTP 403 during this pass.
We therefore do not call the ledger replay a reconstruction of the stochastic
scheduler or an independent distance optimization.

For the mitigation diagnostic we use Figure 10's d_i=9 small-angle slopes:
beta=0.4×5×6.96e−4 for target hardware and 0.4×5×3.90e−4 for desired hardware.
The angle convention is exp(−i theta P). The deepest Hadamard-test branch has
total angle Theta_max=T_max||H||_1/2, with
||H||_1=4L(L−1)+UL²/4. For each time point the exponent is

    chi(f) = 4 beta Theta_max f + 8 I_growth R(f),

where R(f)=2 N_terms ceil(r_max f). The factor eight follows from the two
growth contributions per successful RUS rotation and the PEC factor four.
We first infer the growth probability required to match the published
total/PEC-free ratio. We also fit only the U=4,L=4 anchor for each hardware
and transfer that single probability to the other rows without refitting.
This tests consistency across workloads; it cannot independently validate
the growth model used to produce the original table.

Two desired-hardware entries are retained as discrepancies. At U=8,L=8,
750 clocks predict 470 ms instead of the listed 297 ms. At U=8,L=6, the
published PEC multiplier is below the angle-only value under these stated
conventions. Agreement of a total runtime alone can hide compensating
differences in its PEC-free and mitigation factors.

For FT, the time-optimal throughput T_count d_f t_cycle reproduces all 16
reported time-optimal values within 2.8%. This checks one-port consumption,
not the private TopQAD factory/layout optimizer. The 16 space-optimal entries
are retained explicitly to show their throughput penalty. We do not claim
an independent reproduction of their optimized factory areas.

The fixed-workload sensitivity changes the cycle, angle coefficient or growth
contribution separately. Substituting our growth upper bound changes the noise
and protocol interpretation; it is an illustrative sensitivity, not a validated
hardware translation. It shows why resolving growth matters more than a small
clock correction at that operating point.

## NISQ: use the published measurements

We pin [Zenodo 17799843](https://zenodo.org/records/17799843), version 1.0.0,
with its archive MD5 and our SHA-256 hashes, rather than following a changing
latest-version link. The data are by Phasecraft, Brian Flynn, Ashley Montanaro,
Jan Lukas Bosse and Max Hunter-Gordon, under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Our committed
reduction selects the doublon observables, rescales sums to spatial means and
adds statistics reconstructed from the raw shots.
The 56-qubit experiment is a 7×4 periodic lattice with pi flux and a dimer
state containing a holon and a doublon, evolved for four Trotter steps.

`phasecraft_data.py` reduces the raw shot arrays and supplied processed
expectations. It retains all times and all single-site doublons, plus the
spatially averaged doublon. No TFLO fit or smoothing procedure is rerun.
The output identifies the exact method as well as the archive group: notably,
the group called `exact` also contains approximate Majorana propagation.
Selecting only by that group name would incorrectly label interacting
references as exact.

The raw mean is reconstructed directly from occupation bits using the supplied
qubit ordering. Its sample variance includes correlations between lattice
sites. Raw HDF5 arrays do not provide pseudo-twirling identifiers, so the
pooled independent-shot standard error is a diagnostic only. Figures and
comparison tables preserve the authors' reported errors; they do not replace
them with the pooled estimate.

We compare all 20 experimental times. Local comparisons include all 28 sites.
U=0 uses the separately labeled FLO and Exact references; it is a training
consistency check. U=4 uses TDVP at bond dimensions 1024 and 2048; it measures
reference agreement and includes Trotter and finite-bond errors. Small
differences between successive bond dimensions do not certify convergence.
Neither comparison establishes the interacting residual bias of mitigation.

The mean observable has much smaller shot variance than a single-site
Bernoulli variable. Training and smoothing also use information from other
circuits and times. A gatewise PEC sampling bound for one local observable
therefore does not predict the experiment's shot count. Conversely, the
published experiment does not certify our distinct 0.001 local-observable
contract. The new comparison supplies measured operating points without
inventing a transferable effective mitigation exponent.

Willow raw shots and processed estimates are also public in
[Zenodo 17881822](https://zenodo.org/records/17881822), version 1.0.1. We checked
the metadata and README; the 8.8 GB archive was not reduced in this pass.
`data_availability.json` records the URLs, sizes and checksums for both sources.

## Reproduction

`make published-benchmarks` replays the saved reduced data and table inputs and
regenerates the figures and summary. `make published-check` checks the source
hashes and numerical replay. Both work without network access.

For reduction from the original experimental archive, run
`python python/phasecraft_data.py` in the qb_gsee environment (NumPy and h5py
are required), followed by `make published-benchmarks`. The downloader verifies
the pinned 111.9 MB ZIP checksum and caches extracted HDF5 files under `build/`.
About 1.1 GB of uncompressed cache space is needed. Only the reduced data,
hashes and results are committed. The reader allows only the known NumPy
serialization types; redundant uncertainty objects are discarded in favor of
the explicit numeric value/error columns.
