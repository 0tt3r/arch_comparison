# Willow operating points and sampling design

The [summary](../results/willow_v1/summary.md) and two figures compare the public
superconducting experiment with our physical-qubit resource assumptions.
The gate counts, hardware footprint, attempted shots and processing can be
checked directly. Observable accuracy still depends on the available reference.
This pass leaves the integrated architecture comparison unchanged.

The subsequent [matched-reference calculation](WILLOW_REFERENCE.md) supplies
an independent interacting reference for the 4×4, U=4, zero-flux subset.
The archive limitations below still describe the supplied public data.

## Declared data selection and provenance

We use [Zenodo 17881822](https://zenodo.org/records/17881822), version 1.0.1,
and the pinned [paper](../refs/2510.26845v1.pdf). The declared selection contains
all center-hole cases at 4×4 and 6×6: zero and pi flux, U=0,4,8, and one, two
or three Trotter steps. Each case has 24 nonzero evolution times, 20 randomized
compilations per time, and 10,000 shots per compilation. Thus the selection
contains 36 cases, 864 time points and 172.8 million attempted shots.
The configuration was fixed before comparing observable performance.

`willow_archive.py` reads the ZIP directory and downloads only the 73 selected
members by HTTP byte range: 72 HDF5 files and the initial-state JSON. It checks
the returned range, member name, uncompressed size and ZIP CRC, then records
the SHA-256 of each extracted file. About 1.06 GB is transferred and 2.24 GB
is extracted under ignored `build/willow/`. The full published archive MD5 is
recorded but **not verified**, because the other archive members are not fetched.
The committed reduction uses deterministic gzip and supports offline analysis.

The data are by Phasecraft, Brian Flynn, Jan Lukas Bosse, Max Hunter-Gordon
and Ashley Montanaro, under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Our changes select Z and onsite opposite-spin ZZ, derive doublons and staggered
magnetization, and summarize the raw shot groups. The published trained
estimates are reused; TFLO, smoothing and their fit uncertainties are not refitted.

## What is actually supplied

Three distinctions prevent an incorrect accuracy claim:

- The U=0 `exact` group changes with Trotter step count. It is an ideal free
  **circuit** reference. We independently compute continuous free dynamics
  with the orbital Hamiltonian on the open lattice, and check that method
  against the existing many-body sector solver. Occupation-basis inputs make
  density observables invariant under the choice of magnetic gauge.
- The U=4 and U=8 `exact` groups are empty. The interacting exact comparisons
  in the paper therefore cannot be replayed from this pinned archive. TDVP
  at bond dimensions 1024 and 2048 supplies approximate reference comparisons;
  agreement between bond dimensions does not establish exact truth.
- The initial-state JSON repeats a singlet configuration under the center-hole
  labels, with out-of-range coordinates for 4×4. We preserve and flag it. We
  reconstruct the product-state occupations from the zero-time TDVP Z values,
  require every occupation to be 0 or 1, and check one missing particle. No
  interacting evolution is inferred from this initial-state check.

The paper's Table I counts 34 and 74 physical qubits, including two ancillas.
The headline 32 and 72 fermionic modes do not include those ancillas. Qubit
labels in the shot arrays identify the physical-to-fermionic mapping; ancillas
are excluded when computing occupation observables but included in the footprint.

## Shot grouping and the meaning of raw

The archive README specifies contiguous blocks of 10,000 shots for the same
randomized compilation. We retain these groups. For each group, we calculate
unfiltered estimates and number-selection variants allowing zero, one or
two total deviations from the initial up/down occupations. The acceptance
probability is charged relative to all attempted shots.

Published `raw` values already include selection and diagonal symmetry
averaging. To audit them, each accepted group's mean gets equal weight,
even when acceptance differs between groups. Averaging all retained shots
instead would overweight the compilations with higher acceptance. The
selection and symmetry reconstruction is compared to all published single-Z
and derived local-doublon means at every time. Any unmatched case remains
flagged; the comparison does not choose a filter to improve agreement with
classical truth.

The reconstructed thresholds are zero at 4×4, two at 6×6 with zero flux,
and one at 6×6 with pi flux. They were identified by matching the archived
raw estimates, then checked at every time, interaction and step count;
they were not optimized against an ideal observable.

For group means mu_k, conditional variances s_k² and retained counts n_k,
the measured standard error of the equal-weight group estimator is

    SE² = sample_variance(mu_k) / K.

This includes finite-shot fluctuations. Adding their variance a second time
would change the estimator's uncertainty. Our cluster standard errors are
reported separately from the authors' supplied marginal errors. Neither is
a bound on the remaining systematic error.

To illustrate a sampling-design choice, define

    w = mean(s_k² / n_k),
    a = max(0, sample_variance(mu_k) - w).

Then increasing shots within each existing circuit by h gives the planning
variance (a+w/h)/K, while increasing independent circuit groups by h gives
(a+w)/(Kh). The nonnegative component estimate uses only 20 observed groups.
It is uncertain and assumes exchangeable, independent new circuits with the
same acceptance and conditional variance. It is not a certified noise floor
and does not predict the accuracy of TFLO or smoothing. It makes explicit
which extra sampling the physical resource model would have to buy.

## Observables and reference comparisons

We derive each local doublon as

    D_i = (1 - Z_i,up - Z_i,down + Z_i,up Z_i,down) / 4,

and retain both every local D_i and the spatial mean. The staggered
magnetization uses S_i^z=(n_i,up-n_i,down)/2 and is normalized by L².
All RMS differences include the 24 measured nonzero times and every site;
zero-time constraints are excluded. U=0 remains training consistency, and
U>0 remains approximate-reference agreement, not measured residual bias.

The archive supplies marginal Z/ZZ errors but not their covariance after
processing. For a derived quantity sum_i c_i O_i, we report the bound
SE ≤ sum_i |c_i| SE_i, which holds for any compatible covariance matrix.
This avoids assuming independent trained observables or artificially gaining
a factor sqrt(L²). The figures show the derived central values without
presenting this loose bound as an experimental confidence interval.

The gatewise PEC column uses our specified two-qubit depolarizing channel
at p=0.001 and the paper's gate counts. It does not translate Willow's measured
gate infidelities into that channel, fit a mitigation exponent, or infer a
sampling penalty from the experimental reference discrepancy.

## Reproduction

```bash
make willow-benchmark PYTHON=/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python
make willow-check PYTHON=/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python
```

These commands analyze the committed reduction without network access.
To repeat the original retrieval and reduction, use:

```bash
python python/willow_archive.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python python/willow_data.py
make willow-benchmark
```

NumPy, SciPy, Matplotlib and h5py are provided by `qb_gsee`. The source HDF5
reader uses the restricted serialization reader from the preceding benchmark.
No exact many-body dynamics job or quantum hardware job is submitted.
