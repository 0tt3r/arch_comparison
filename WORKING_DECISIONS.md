# Accepted conventions and sequence — 2026-09-22

This records the user's decisions after the merged PLAN.md. It supersedes that
document's precision default and implementation sequence. PLAN.md currently has an
editor lock, so it has been left intact.

C1–C3, C5–C8, and C10–C15 are accepted defaults, with C12 sequenced as below.
C4 is **1% relative error in the nonzero local doublon expectation**, with 95%
confidence on the statistical component and a budget for systematic bias.
Campbell's relative total-energy accuracy is fixed absolute energy-density accuracy;
local doublon density is already intensive. No factor of N should be introduced
into its allowed error. The first plots use the explicit provisional signal scale
D* = 0.1, and include sensitivity to D*; this is not a calculated dynamics value.

C9 uses a per-attempt logical-fault budget. For D in [0,1], worst-case unmitigated
bias is at most q_L. For the centered PEC estimator it is at most Gamma q_L if
q_L bounds each sampled branch. The reference includes this sufficient worst-case
amplification. A sharper q_L bound is retained as a sensitivity when fault location
and commutation permit cancellation without amplification; independence alone is
insufficient. See CUDA_Q_RESULTS.md for the condition and a counterexample.

The agreed sequence is:

1. Write the analytical forms, derivations, conventions, and sources in LaTeX,
   with a dedicated BibTeX file.
2. Plot those forms to understand scaling. Keep r as a shared circuit parameter;
   vary it without first claiming a calibrated Trotter-error law.
3. Use CUDA-Q Logical for the common workload and each architectural path. Replace
   gate, routing, timing, and supply proxies with compiler-supported counts.
4. Refine Trotter error using analysis and selected small-system dynamics.
5. Bring in qubitization as a separately identified algorithm comparison.

Campbell-style second-order plaquettes are the default algorithm candidate, adapted
to open boundaries. Termwise second order is retained as a useful comparison.
Their different compiled angle distributions matter: the plaquette construction
includes repeated fixed-angle basis changes. Hamming-weight phasing, routing,
cancellations, and hybrid distilled fixed angles remain refinements.

The analytical first pass is in analytical_forms.tex, references.bib, model.py,
assumptions.json, and results/analytical_v1/. Additional numerical choices in
assumptions.json are declared plotting assumptions, not additional user-approved
hardware specifications. Exact dynamics has not been run.

The first CUDA-Q compilation pass is now in cudaq_compile.py and results/cudaq_v1/.
It preserves the precision and shared-r decisions. Native physical gate counts,
two explicit NN realizations, signed angle classes, and FT synthesis counts are
compiled. Encoded timing, factories, and STAR injection remain reference models;
the complete Hubbard P3 macro schedule failed a timing sanity check and is not
used in the plots. This boundary is documented in CUDA_Q_RESULTS.md.

For the new compiled physical comparison, p per 100 ns of idle time is the default
model refinement, with p per 10 ns layer retained as a sensitivity. The original
assumptions.json and analytical_v1 artifacts are unchanged. Neither convention
has been calibrated against a common physical syndrome-extraction circuit.
The fixed-angle-distilled/analog-dynamic hybrid is now included in the reference
cost model. Trotter-error and dynamics calculations remain later work.

The second refinement is recorded in NISQ_REFINEMENT.md and nisq_refinement.tex,
with dedicated refinement_references.bib and results/nisq_v2/. It compiles native
RZZ circuits, exact fermionic observable cones, and spin-parity postselection
acceptance. Arbitrary RZZ and one-qubit pulses each take 10 ns in this explicit
native-gate scenario. The main idle rate remains p per 100 ns, with zero and
one-tenth-rate sensitivities; these are modeling choices, not new user hardware
specifications. Postselection has not met the bias target and receives no credit
in the PEC curves. Exact cone restriction is available to every architecture.

Distillation input error is now a separate parameter. The new sensitivity compares
p, 34p/15, 46p/15, and a conditional improved-state 2p/15. Literature coefficients
depend on their source circuit/noise model and are not a calibration to our
p_1q=p/10 convention. At L=4, conventional injection forces two-level hybrid
distillation and reverses the earlier gain under the retained serial factory
reference. The old model defaults and original v1 results remain reproducible.
The accepted precision target and the conditional C9 treatment are unchanged.
