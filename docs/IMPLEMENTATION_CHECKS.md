# Estimator, injection, FT arithmetic, and locality checks

The four calculations in `results/implementation_v1` test assumptions that
could change the architectural comparison. They preserve the frozen capacity
curves. A measured small-patch bias, a generator-fault injection channel, a
logical dependency schedule, and a finite-reference difference have different
scopes; none is automatically a replacement for a complete architecture model.

## NISQ estimator calibration

The task is the doublon at coordinate (1,1) of a 2×2 open Hubbard lattice,
starting in the same Néel state. The simulator evolves the full eight-qubit
density matrix, including leakage between spin-number sectors. Native Rzz
circuits use the existing plaquette compiler and its nearest-neighbor routing.
The second-order step counts are 4, 8, 16, and 32. Gate noise is a uniform
nonidentity Pauli mixture with probabilities p₂=p and p₁=p/10. Preparation and
readout bit flips have probability p. Idle noise scales with time: p per
100 ns, including p/10 per unoccupied 10 ns layer. This convention scales with
the swept p; the configuration's `idle_p_per_100ns` is its p=0.001 anchor.

For each p and r, six U=0 circuits at tau=0.15,0.3,0.55,0.8,1.2,1.6 train an
affine map from noisy to ideal expectation. The U=0 interaction slots remain
present with zero angle and the same noise. Test circuits at U=2,4 and
tau=0.5,1,2 are held out. The mismatch test increases p₂ by 20% and all onsite
angles by 1%, leaving the training unchanged. Raw, spin-number selected,
trained-affine, and selected-plus-trained estimators are all reported.

The affine fit uses population training expectations to isolate residual
systematic bias. Finite training is charged: six million attempts per fit,
plus propagation of binomial uncertainty through the fitted intercept/slope.
If x is a training expectation and y its ideal label, X=(1,x), beta is the OLS
fit, e=y-X beta, and h=(1,x_test), its sensitivity is

```
d prediction / d x_i = h (Xᵀ X)^(-1) ([0,e_i]ᵀ - beta_slope X_iᵀ).
```

The selected-observable variance per physical attempt is D(1-D)/A, with A
the acceptance. The squared slope multiplies the test variance. Statistical
planning solves 1.96²(V_test/S + V_training) ≤ 0.0025². Training time is charged
to every displayed test campaign rather than amortized over many tests.
A separate 8000-replicate binomial training simulation checks this delta-method
variance. The total plan adds absolute mitigation bias, absolute Trotter bias,
and 0.0025; it cannot profit from cancellation of systematic errors. This is
Gaussian planning under a fully specified simulated channel, not a rigorous
hardware confidence statement. No 4×4 mitigation multiplier is fitted from it.

## STAR preparation and growth

The channel calculation follows the multirotation mechanism in Toshio et al.
(`refs/PhysRevX.15.021057.pdf`) while separating the faults it actually includes.
For odd distance d, k=(d+1)/2 physical factors use alpha=atan(|tan theta|^(1/k)).
The ideal projection acceptance is cos(alpha)^(2k)+sin(alpha)^(2k). Independent
rotation-generator faults have probabilities p/15 on weight-two factors and
p/30 on the final weight-one factor, consistent with p₁=p/10. Perfect syndrome
projection yields an explicitly enumerated mixture of logical angle errors.
Even-k logical Clifford corrections are charged. The other Pauli faults,
preparation faults, and syndrome-extraction faults are excluded from this
channel calculation; it is not a full injected-state error model.

For error-angle mixture weights w_j and angles delta_j, c=sum w_j exp(-2i delta_j).
Ideal probabilistic coherent-error cancellation (PCEC) composes the channel
with its conjugate, giving |c|². A one-level noisy-PCEC experiment implements
the sampled inverse angles with raw RUS channels instead. Their actual angles,
source acceptance, growth noise, retries, and delivery time are charged. The
remaining imaginary component is retained. After the dephasing inverse, its
observable bias is bounded by 0.5*(product_j (1+abs(Im c_j)/Re c_j)^R_j - 1),
which includes amplification by the inverse and composition of non-CPTP maps.
This is not an assumption that every noisy correction leaves a Pauli channel.

Main RUS branches succeed with probability 1/2. The workload chooses a retry
cap reserving 1% of epsilon for its omitted branches, with the configuration
value 12 serving only as a minimum; correction RUS has cap 32. Exact Clifford
tails terminate early. Independent growth dephasing is composed per attempted
resource. Protocol switching chooses the original-injection error reference
when it is smaller **after** the candidate's growth error. Its per-state error
2p/15 and acceptance 0.9 are declared reference assumptions, not new simulation
results. The conventional full-RUS 4p/15 is shown only as a guide.

The resource-time scenarios use

```
acceptance = projected_acceptance (1 - detector_factor*p)^(d²+6d)
trial rounds = (4 + even_k*d)/acceptance + growth_rounds + 9*d_data
```

with detector_factor 1 or 4, growth_rounds=d_data+1 for d<d_data, and
independent growth error

```
f(d,p) = min(0.25, f9*(p/0.001)^5*(d/9)²*(0.01/p)^(-(d-9)/2)).
```

The free f9 values 0, 1e-5, and 1e-4 are anchored at d=9, p=0.001.
The power-law extrapolation makes f scale as p^((d+1)/2) at each fixed
distance; these are scenarios rather than a measured p law. This suppression law is not
fitted to the small growth probes, and the detector factors are not calibrated.
Chung et al.'s reported fit uses a different noise model and is not transplanted.
Preparation is serial; supply, telemetry and PCEC costs are first moments, not
a parallel gadget schedule or a physical layout. The state source and logical
memory errors outside these components remain outside this new submodel.

The independent Clifford experiment uses Stim's four-layer nearest-neighbor
rotated-patch extraction circuit. An ideal d_i=3 or 5 patch is expanded directly
into d_f=5 or 7 by preparing the right strip in Z and bottom strip in X, then
running d_f noisy extraction rounds. No-growth controls use d_i=d_f. Exact
stabilizer flows determine transition checks; row reduction on the new checks
keeps the detector basis local enough for graph decomposition. Later checks
compare adjacent rounds. The experiment includes gate, preparation,
measurement/reset, and idle faults under the same p₂/p₁ convention. The
60 ns gate sequence, 100 ns measurement, and 80 ns reset/feedforward idle
allowance make the full 240 ns cycle. It measures
X/Z logical memory and decodes with PyMatching. Each of 20 cases has a fixed
500,000 shots, with exact binomial two-sided 95% intervals. Zero failures are
upper limits. Noiseless and deliberately inserted logical-flip controls check
both detectors and logical observables.

These probes include stabilization and readout, but exclude noisy small-patch
preparation. They neither isolate a growth-only floor by subtraction nor
establish a distance-25 arbitrary-angle injection channel. The full transition
circuits and seeds are saved for extension.

## FT arithmetic and delivery

The explicit test is one equal-angle onsite layer on 16 sites, with
exp(-i theta Z_up Z_down), theta=1/36. Parity CNOTs convert it to 16 equal Rz
rotations. A binary compressor network uses b-popcount(b) clean carries per
batch. Each temporary AND has four T gates and is an exact isometry on a clean
target. X measurement plus a conditional CZ uncomputes it without extra T.
Tests check all binary inputs through b=16 and arbitrary superpositions in
every measurement branch through b=4. CUDA-Q Logical compiles the complete
arithmetic, measurement and feedforward program, then GridSynth synthesizes
rotations at precision 1e-8. Full-program T counts agree with the per-angle
inventory plus explicit arithmetic.

A constructive router places logical patches on a square grid and replaces
remote two-patch gates with adjacent SWAP ladders. Dependencies, measurement
and feedforward are explicit. The declared macro costs are three d-round
clocks per adjacent CX/CZ, one per H/S or T consumption, two for X measurement
and reset, and one additional feedforward clock. Synthesis Clifford gates are
charged serially after each synthesized T sequence, a conservative scheduling
choice. This checks a logical schedule; lattice-surgery ancillas, decoder
failure rates and physical surgery occupancy are not simulated. Magic-state
ports constrain availability but their transport paths to the data are not
geometrically routed; these also remain outside the reported stage timing.

Source ports emit one state per 150 physical rounds, with d=25 and 240 ns
rounds. The bounded policy stores one state per port and pauses production
when occupied. Continuous production reports its required peak storage.
The patch plot includes packed data/routing patches and state buffers but
excludes source and surgery footprints. It is not a full area estimate.
This experiment diagnoses routing and buffering penalties; it does not certify
that the existing Hamming-weight full-campaign improvement survives an optimized
physical implementation, or that this greedy router is optimal.

## Finite-patch locality

The exact patch hierarchy is cross5, 2×2, 2×3, 3×3, 4×3, and 4×4 where data
exist, all with the same observed coordinate (1,1), Néel pattern, U=4 and OBC.
Only strict set inclusions are compared. Smaller patches are recomputed; saved
large exact results are reused. Spatial allowances are 0.00125 and 0.0025.
Differences are errors relative to the named finite reference, not infinite
lattice bounds. In particular, coordinate (1,1) retains a boundary at distance
one even in 4×4. These calculations test proposed truncations and show why
Trotter-coefficient convergence does not imply observable convergence.

## Reproduction

Science uses `qb_gsee` (including Numba), compilation uses `.venv-cudaq`, and
circuit-level growth uses the existing `stim` environment with Stim and
PyMatching. Versions are saved with results. No GPU or new source build is needed.

```bash
make implementation PYTHON=/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python \
  QEC_PYTHON=/home/ac.otten/soft/mambaforge/envs/stim/bin/python
make implementation-check QEC_PYTHON=/home/ac.otten/soft/mambaforge/envs/stim/bin/python \
  PYTHON=/home/ac.otten/soft/mambaforge/envs/qb_gsee/bin/python
```

`make implementation` regenerates the four bounded studies, compiler artifacts,
figures and summary; it can take several minutes. NISQ cases resume only with
matching source/configuration hashes. Use a fresh result directory after
changing its numerical protocol. `make implementation-figures` replays saved
results and redraws the four figures without simulating channels or circuits.
The ordinary paper build uses this saved-data path.
