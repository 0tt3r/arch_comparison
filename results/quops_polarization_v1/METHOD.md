# Exact calculations and reproducibility

The requested individual-circuit sweep is complete for widths 2–6, with 20
circuits per width and every integer prefix depth for all three error rates.
At widths 7–8, individual-circuit checks cover only depths 1 and 2. The full-depth
mean at widths 2–8 is instead computed by exact analytical ensemble integration.
The deep individual-circuit requirement at widths 7–8 remains incomplete; the
angle-integrated samples must not be relabeled as concrete-circuit fidelities.

## Noise and circuits

Each layer is `w` simultaneous rotations, followed by a uniformly random
oriented matching of `floor(w/2)` CNOTs. The shuffled unmatched qubit idles for
10 ns at odd width. Every active gate lasts 10 ns. There is no routing, SPAM,
measurement, reset, feedforward, or additional waiting interval.
The order follows the brief; arXiv:2609.12146v1 instead places CNOTs before
rotations. Gates are `exp(-i theta P)`, with independent uniform axes and
`theta` uniform on `[0, 2*pi)`, consistent with PLAN.md C5. This gives the same
rotation ensemble up to global phase as uniform conventional Bloch angles.

For a k-qubit gate, the noise is

\[
\mathcal D_{k,p}(\rho)=(1-p)\rho+
 \frac{p}{4^k-1}\sum_{P\ne I}P\rho P.
\]

Thus `p` is process infidelity (nonidentity Pauli probability), rather than
average gate infidelity or the probability of replacing the state by the
maximally mixed state. Two-qubit gates use `p`; single-qubit gates use `p/10`.
The explicit brief's time-proportional idle convention gives `p/10` per 10 ns
idle. This uses linear probability scaling; it is not a continuous-time
semigroup conversion. PLAN.md C6's default `p` per idle layer is different and
is not used. The comparison with idle disabled is paired on the same circuits.

## Concrete circuits: full overlap propagation

Use orthonormal normalized Pauli operators. If `A` and `B` are respectively the
noisy and ideal Pauli transfer matrices, propagate `M=A B^T`, starting at the
identity matrix. An ideal gate with orthogonal transfer matrix `R`, followed by
noise with diagonal transfer matrix `L`, gives

\[
M' = L R M R^T,\qquad
F_{\rm proc}=\frac{\operatorname{Tr}M}{4^w},\qquad
\gamma=\frac{\operatorname{Tr}M-1}{4^w-1}.
\]

Noise acts on the left only. Rotations mix two Pauli indices; CNOTs are signed
permutations. No Kraus trajectories, stochastic traces, truncations, or pure
state-fidelity proxies are used. Arrays and accumulation use float64.
The density-matrix reference uses a normalized maximally entangled state,
Cirq 1.4 channels acting on one half, and the ideal final Choi vector, all in
complex128. Its agreement with the overlap implementation is recorded in
`validation.json`, including zero noise, idle on/off, and a larger positive
control error rate.

At D=1, M remains diagonal, so the D=2 trace needs only the diagonal of
`R M R^T`. Squared rotation matrix elements therefore give an exact shortcut
for the first two concrete-circuit prefixes. This is checked against the full
matrix at widths 2–5. It is deliberately prohibited for D>2.

## Full ensemble: exact analytical averaging

Starting from a diagonal M, averaging an independent uniform rotation angle
kills its off-diagonal entries. The two noncommuting Pauli entries are replaced
by their mean. Averaging the axis as well gives the local diagonal map

\[
v_I'=v_I,\qquad v_a'=\frac23v_a+
  \frac16\sum_{b\ne a,\ b\ne I}v_b \quad(a=X,Y,Z).
\]

Linearity and independence allow this average before each next layer. Noise
and CNOTs preserve diagonality. Full ensemble averaging is also invariant
under qubit permutations. A random oriented matching can be written as a
fixed matching conjugated by a uniform permutation. For permutation-invariant
input its exact average is therefore a fixed matching followed by projection
onto permutation orbits. These orbits are labeled by `(nX,nY,nZ)`; this is
**not** a Pauli-weight-only approximation. The odd-width idle qubit is included
in the permutation. All `4^w` diagonal entries are retained. The cost is
O(D w 4^w), rather than the full overlap matrix's O(D w 16^w).

Validation independently sums 4096 concrete angle-quadrature circuits and
explicitly enumerates all oriented matching permutations at widths 3 and 4.
The quadrature is exact for the first and second sine/cosine moments, not a
numerical approximation to the uniform integral. Twenty additional sampled
axis/matching sequences per width are integrated over angles only. Their raw
method label is `exact_angle_mean`; their spread is conditional-mean spread,
not the full concrete-circuit distribution.

## Products and fits

The local gate-polarization product is

\[
\gamma_{\rm local}=(1-16p/15)^{n_2}
 (1-4p/30)^{n_1+n_{\rm idle}}.
\]

Embedding a noisy k-qubit gate in w qubits leaves its process infidelity equal
to p, so its **w-qubit** polarization is `1-a_w p`, where
`a_w=4^w/(4^w-1)`, independently of k. The alternative product is therefore

\[
\gamma_{\rm embedded}=(1-a_wp)^{n_2}
 (1-a_wp/10)^{n_1+n_{\rm idle}}.
\]

Products of embedded gate polarizations are still a model: the simulation
retains return/cancellation correlations. The exact ensemble mean is computed
without either product or a fitted constant. Fit `gamma=exp(-c*x)` with
`x=p[n2+(n1+n_idle)/10]`, zero intercept, and equal total weight for each
width/error-rate cell. Train on p=0.001 and 0.01, with idle enabled; hold out
p=0.003. Both a shared c and one c per width are reported. Circuit depths are
correlated prefixes; no independence across depths or noise rates is assumed.
Sample intervals are pointwise t intervals from 20 circuits; they are not a
simultaneous confidence band or a hardware validation.

## Reproduction and files

Use `/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python -B` followed by
`python/quops_polarization.py` and one of these commands:

```
validate
ensemble
campaign --workers 8 --replicas 20 --exact-max-width 6
shallow
replay
analyze
```

Set `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1` for the campaign.
Run substantial work in a PBS allocation using the provided campaign script;
`submission.json` records the actual command and frozen source. Source
snapshots live in `sources/`; to execute a snapshot set `QUOPS_ROOT` to the
repository root. The campaign is restartable via completed hash-checked
receipts. It does not reuse a raw file without a completion receipt.

`ensemble.csv` is the exact mean, `comparison.csv` includes both products and
fits, `fit.json` records training/holdout results, and `endpoints.csv` records
the first crossing of 0.195 for every width/error rate. `raw/` and `circuits/`
contain circuit-level outputs and axes/angles/matchings. `sample_summary.csv`
keeps all calculation methods separate. `audit.json`, `validation.json`,
`replay.json`, terminal PBS accounting, and `MANIFEST.sha256` provide the
completion, correctness, and provenance checks. All result writes are confined
to this directory; the only other created paths are the standalone Python
entry point and `docs/QUOPS_POLARIZATION.md`.

Reference: [Proctor et al., arXiv:2609.12146v1](https://arxiv.org/abs/2609.12146v1),
especially the process-fidelity definition and QUOPS circuit ensemble in SI
Sections I B and II B. Gate-product predictions are tested here; the paper's
full compilation and experimental mirror protocol are not implemented.
