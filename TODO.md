# Open work

## Circuit-level simulation of STAR injection (deferred)

**Why.** Raw repeat-until-success supplies most of STAR's gain over PCEC (about 10² in qubits at 16×16,
τ = 1), and it rests on a first-order cancellation of the coherent error between successful and failed
teleportations ([docs/RAW_RUS_NOTE.md](docs/RAW_RUS_NOTE.md)). Our injection instrument models only
generator faults on the transversal rotations (p/15 per pair factor, p/30 for the last). Chan et al.
(arXiv:2609.17706) show that in T-state cultivation, hook errors propagate to *coherent* Clifford errors
after the final check and degrade the fault distance. An analogous hook-to-coherent path in STAR's
preparation would not be conjugated on failure and would break the cancellation.

**Expected impact.** STAR and hybrid only. The manuscript already brackets it with the `pcec_only`
configuration (round v3): if the cancellation fails, STAR falls to the PCEC line. Full fault tolerance's
crossovers do not depend on it.

**Plan.**
1. Build Toshio et al.'s transversal multi-rotation preparation (m = 2, d_i = 3 and 5) and the ZZ
   lattice-surgery teleportation, with the Clifford parts in Stim.
2. Treat the k = (d_i + 1)/2 small rotations as an explicit sum over 2^k Pauli branches (cos θ I + i sin θ ZZ),
   and enumerate every single fault (and every fault pair at d_i = 3) exactly, in the style of Wan and
   Zapirain (arXiv:2609.18922). Sampling cannot resolve coherent terms of order pθ ~ 10⁻⁵–10⁻⁷.
3. Per accepted outcome, record the logical resource channel: dephasing, coherent part by branch, and the
   teleportation-outcome bias s − 1/2. Add the surgery outcome-flip rate from results/surgery_v1.
4. Compose the RUS chain and compare raw RUS with PCEC; compare the preparation channel with our
   instrument and with Toshio et al., Fig. 7. Validate the harness on T-state injection at θ = π/8
   against the exact methods of Takada et al. (arXiv:2609.16929) where they overlap.

**Tools.** Stim for Clifford parts; PECOS (`pip install quantum-pecos`) or a direct state vector as an
independent check at d_i = 3. The efficient exact methods in the references target the third level of
the Clifford hierarchy (T, CCZ) and do not apply to arbitrary small angles directly.

**Scale.** Several days; a suitable bounded brief for ASTRA, reusing the Stim surgery setup.
