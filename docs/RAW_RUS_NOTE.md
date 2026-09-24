# Raw RUS versus PCEC: what cancels, and what must be checked

**STAR discussion note, 2026-09-23.** Reproduced and checked against the paper text on the same day;
`make raw-rus-check` reruns every table entry and control (about 70 s). This is a statement about the restricted
[`star_channel_study.py`](../python/star_channel_study.py) instrument, not a
validation of complete noisy injection. Angles use `exp(-i theta Z)`.

**Derivation.** At retry index \(m\ge0\), let the accepted error channel multiply
off-diagonal density-matrix elements by \(c_m=1-a_m+i b_m\), with
\(a_m,b_m=O(p)\). Assume independent resource preparations/retry faults;
state- and history-independent success probability exactly \(1/2\); failure
multiplier \(c_m^*\); ideal feedforward and angle doubling
\(\theta_m=\operatorname{wrap}(2^m\theta)\), with Clifford wrapping preserving
this conjugation; and no further branch-dependent selection. Discard the retry
record after feedforward. Each branch implements the same ideal rotation, and
the averaged error multiplier is

\[
C=\sum_{n\ge0}2^{-n-1}c_n\prod_{m<n}c_m^*
 =1-\sum_{m\ge0}2^{-m}a_m+O(p^2),\qquad
[b_m]:\;2^{-m-1}-\sum_{n>m}2^{-n-1}=0.
\]

Thus **coherence cancels to first order in fault probability**, not exactly:
\(\operatorname{Im}C=O(p^2)\), while dephasing remains \(O(p)\). The expansion
requires uniformly small channel errors over the retries. Doubling ensures the
common ideal rotation; the cancellation algebra itself allows any sequence
\(c_m\). It describes the averaged quantum channel, not each recorded branch.
An additional fixed growth error \(f\) can leave mixed \(O(pf)\) terms.

**Numerical check.** Twenty-seven angles in \([10^{-3},\pi/8]\), including
0.01 and \(\pi/16\); 162 parameter points; no growth noise; RUS cap 48
(the implemented PCEC correction uses cap 32). Ranges and maxima below are
over sampled angles, not certified continuum bounds. \(q=(1-\operatorname{Re}C)/2\).

| \(d\) | \(p\) | Raw \(q\) range | Noisy-PCEC \(q\) range | max \(\lvert\operatorname{Im}C_{\rm raw}\rvert\) | max \(\lvert\operatorname{Im}C_{\rm PCEC}\rvert\) |
|---:|---:|---:|---:|---:|---:|
| 5 | 1e-3 | 8.19e-7–1.02e-4 | 1.64e-6–2.04e-4 | 7.80e-9 | 2.26e-12 |
| 5 | 1e-4 | 8.19e-8–1.02e-5 | 1.64e-7–2.04e-5 | 7.80e-11 | 2.26e-15 |
| 9 | 1e-3 | 1.04e-6–1.65e-4 | 2.08e-6–3.30e-4 | 2.66e-8 | 1.70e-11 |
| 9 | 1e-4 | 1.04e-7–1.65e-5 | 2.08e-7–3.30e-5 | 2.66e-10 | 1.70e-14 |
| 25 | 1e-3 | 2.50e-6–4.35e-4 | 5.02e-6–8.70e-4 | 2.32e-7 | 4.46e-10 |
| 25 | 1e-4 | 2.50e-7–4.36e-5 | 5.01e-7–8.72e-5 | 2.33e-9 | 4.48e-13 |

Reducing \(p\) tenfold reduces raw coherence approximately 100-fold. For
\(\chi=-2\log\operatorname{Re}C\), PCEC/raw ratios span 1.9992–2.0068.
The review's smaller ratios, 1.52 and 1.57, reproduce after composing its
growth upper-limit scenario at 4×4 and 16×16; they are not universal constants.

**Acceptance and failure conditions.** `preparation()` returns a mixture of
equatorial states, hence teleportation outcomes remain equiprobable even at
low preparation acceptance \(A\). At \(\theta=0.01,p=10^{-3}\),
\(A=(0.87269,0.47915,0.0054877)\) for \(d=(5,9,25)\): approximately
1.15, 2.09, and 182 projection attempts before additional detector rejection,
**not biased teleportation**.

- **Biased or nonconjugate branches.** If success is \(s_m=1/2+\eta_m\), with
  reach probability \(R_m=\prod_{j<m}(1-s_j)\), then
  \(\operatorname{Im}C^{(1)}=2\sum_m R_m\eta_m b_m\).
  At \(\theta=.01,d=9,p=.001\), a synthetic constant \(s=.501\) changes
  \(|\operatorname{Im}C|\) from 3.54e-9 to 8.60e-8. If \(\eta=O(p)\), this
  is still second order. Omitted preparation faults can also invalidate the
  scalar-channel assumption: in the standard teleportation circuit,
  \(s-1/2=z_Az_D/2\). An independent omitted fault of probability \(u\),
  always accepted and maximally polarizing the ancilla, permits
  \(|s-1/2|\le u/[2(u+(1-u)A)]\). Taking **hypothetically** \(u=p=.001\)
  gives 0.000573, 0.00104, and 0.0771 at the three acceptances above.
- **Correlations/systematics.** A fixed small systematic angle still cancels
  at first order if failure conjugates it. Same-sign success/failure errors
  do not: a 1e-4-rad error gives \(|\operatorname{Im}C|\simeq4e-4\).
  Rare correlated finite errors also evade the independent-fault expansion:
  a common 0.2-rad error on every retry with probability .001 gives
  \(|\operatorname{Im}C|=4.67e-5\), versus 6.15e-8 for independent occurrences.
- **Finite cap.** Conditional completion before \(K\) nonterminal trials leaves
  \(\operatorname{Im}C_K^{(1)}=2^{-K}\sum_{m<K}b_m/(1-2^{-K})\).
  Uncompleted probability is \(2^{-K}\); for a [0,1] observable the conditional
  versus eventual-completion bias before mitigation is at most this probability;
  inverse-channel weighting requires its own amplification allowance. Cap 48 versus
  64 changes the sampled raw multipliers by at most 1.12e-16. A wall-clock
  cutoff correlated with preparation outcomes needs a separate analysis.

**PCEC and the paper.** [Toshio et al.](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.15.021057),
§IV C, Eqs. (21)–(22), suppress the *single-trial* coherent channel to obtain
composable worst-case error control (Appendix D). Ideal PCEC maps each
\(c_m\mapsto|c_m|^2\), doubling dephasing to leading order without other noise.
For our phase channel, half diamond distance is
\(\sqrt{q^2+(\operatorname{Im}C/2)^2}\); after dephasing inversion, a [0,1]
expectation instead carries coherent-bias allowance
\(|\operatorname{Im}C|/(2\operatorname{Re}C)\), accumulated across gates.
Cancellation also reduces the **averaged channel's diamond distance**; it is
not justified simply by choosing an expectation value. Appendix E/Fig. 20
averages over RUS branches exactly as above, \(\tilde E=\sum_K 2^{-K}E_K\), yet reports a raw-RUS
prefactor \(\alpha'_{\rm RUS}\) "slightly larger" than the PCEC value of Fig. 11 in the small-angle
region, with a logarithmic dependence on the target angle. Their averaged coherent term therefore
stays first order, unlike this instrument, where raw RUS carries about half of PCEC's dephasing. The branch/sign/Clifford-wrapping conventions and channel
approximations therefore need reconciliation with the authors. **Use raw RUS
as a budgeted model option; do not infer that physical PCEC is unnecessary.**

**Reconciliation (2026-09-23).** [`python/rus_reconciliation.py`](../python/rus_reconciliation.py)
implements the paper's own Section IV model: a resource prepared at angle a > 0 is wrong with
probability P_L = a^{2/k}P_ud (P_ud = kp/15) and then over-rotates by Δ(a) = −a^{1−2/k} (their
Eqs. 11–16, stated for positive angles); completion at trial K composes
E_{2^{K−1}θ}∘E_{−2^{K−2}θ}∘⋯∘E_{−θ} (Eq. 19), with π/4 wrapping; the averaged channel is
Σ_K 2^{−K}E^K. The only input the paper leaves open is Δ at the negative angle of a failed trial.

- **Conjugate** (physical: a failed teleportation implements minus the resource's *actual* angle,
  so its error is −Δ): the averaged coherent term vanishes to first order (|ỹ|/x̃ ≈ 10⁻¹⁷ over
  θ = 10⁻⁴–π/8 and k = 2–6), and raw RUS carries exactly half of the PCEC dephasing. This is our
  instrument's result.
- **Same-sign** (Δ evaluated at |angle|): the coherent terms add over the ~log₂(1/θ) retries,
  and raw RUS exceeds PCEC at small angles with a logarithmic dependence on θ (for k = 4, α′/k
  rises from 0.17 at θ = 0.05 to 0.70 at 3.5×10⁻⁴). These are the two features reported in
  Appendix E and Fig. 20.

The disagreement is therefore a sign convention for the failure branch, and the physical sign
gives first-order cancellation. With our asymptotic P_L and Δ the model reproduces the paper's
α_RUS ≈ 0.4k (Fig. 11) at k = 3 but not uniformly in k, because the paper evaluates exact
expressions; the comparison between conventions is internal to one model and does not depend on
this. The hardware caveats above (biased branches, correlated or same-sign systematics, omitted
preparation faults) are not addressed by this reconciliation and still bound the raw-RUS gain.

**Consequence for the optimization round.** The `raw_rus` switch in
[OPTIMIZATION_ROUND.md](OPTIMIZATION_ROUND.md) changes the capacity leader in only one of 64 cells,
but it supplies most of the STAR/hybrid capacity gain, for example 16×16 STAR at τ = 1,
ε = 0.01 from 10¹⁰·¹ to 10⁸·⁰ qubits. Those gains hold only if the averaged coherent term cancels
at first order, as it does in this instrument. Under the paper's Appendix E behavior, raw RUS
would be no better than PCEC. The noisy-PCEC configurations are then the relevant STAR/hybrid
lines, and the FT-leads conclusions are unchanged.

Reproduce all rows, controls and input SHA-256 hashes without writing files:
`~/soft/mambaforge/envs/qb_gsee/bin/python -B python/raw_rus_check.py --json`.
