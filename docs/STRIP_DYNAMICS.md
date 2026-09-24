# Interacting Trotter coefficients beyond 4×4

The Hubbard step counts rest on a doublon coefficient C = |D_r − D| r² measured at U/J = 4 on
patches up to 4×4. The termwise value grew from 3×3 to 4×4 (0.077 → 0.19 at τ = 1), so it was
unclear whether larger lattices would raise it further. This note extends the exact dynamics.

## A tensor network does not reach far enough

`python/hubbard_mps.py` evolves the same termwise product formula as a fermionic MPS (one site
per lattice site in the solver's snake order, fermionic swaps for vertical bonds). Without
truncation it reproduces the exact sector solver at every site on 2×2 and 3×2 to 10 digits. It
is still not usable here: the quench from the Néel state at U/J = 4 builds near-volume-law
entanglement by τ = 1, and on 4×3 a bond dimension of 256 already discards 43% of the weight.
Reaching 6×6 with doublon errors near 10⁻⁵, which the coefficient requires, is out of reach.

## Exact strips (`python/strip_dynamics.py`, `results/strip_dynamics_v1/`)

The sector solver is exact up to 18 sites once its site limit is raised
(`cpp/hubbard_strip.cpp`, a copy of the hash-pinned `cpp/hubbard_sector.cpp` with only that
limit changed; 6×3 has 2.4×10⁹ amplitudes). We ran every strip of width 2 up to 8×2 and of width 3
up to 6×3, plus 4×4, at τ = 1 and 2, with Chebyshev exact references. The table gives the central
coefficient (the larger of the two central sites for even sizes); r = 16 at τ = 1 and 32 at τ = 2.

| Patch | 3×2 | 4×2 | 5×2 | 6×2 | 7×2 | 8×2 | 3×3 | 4×3 | 5×3 | 6×3 | 4×4 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| termwise, τ = 1 | .063 | .114 | .189 | .181 | .187 | .186 | .077 | .140 | .180 | .191 | .188 |
| plaquette, τ = 1 | .259 | .252 | .254 | .252 | .249 | .248 | .218 | .246 | .269 | .262 | .262 |
| termwise, τ = 2 | .662 | .557 | .672 | .643 | .678 | .684 | .382 | .702 | .713 | .755 | .775 |
| plaquette, τ = 2 | 1.35 | 1.10 | 1.04 | 1.03 | 1.06 | 1.09 | .935 | 1.14 | 1.06 | 1.08 | 1.10 |

Figure: [figures/strip_dynamics.pdf](../figures/strip_dynamics.pdf).

**At τ = 1 the termwise coefficient converges in length.** Along width-2 strips it rises to 0.19
at length 5 and stays within 0.181–0.189 to length 8; width 3 reaches 0.180–0.191 by length 5–6;
4×4 gives 0.188. The growth from 3×3 to 4×4 was a short-patch effect, not the start of a trend.
The plaquette coefficient is flat at 0.25–0.27 throughout.

**At τ = 2 length convergence holds, width convergence is not established.** At fixed width the
coefficients settle by length 5–7, but the termwise value rises slowly with width (0.68 at width 2,
0.76 at width 3 and length 6, 0.78 at 4×4).

**The planning coefficients cover every measured value.** They are 0.21 (termwise) and 0.29
(plaquette) at τ = 1 and 0.86 and 1.25 at τ = 2, that is 1.1 times the running maximum over the
earlier patches. The largest new value is 11% below the τ = 2 termwise planning value, so no step
count in the optimization changes.

## Scope

Strips probe convergence in one direction at a time. Because the free-fermion scan
([FREE_FERMION_TROTTER.md](FREE_FERMION_TROTTER.md)) shows convergence once a patch is wider than
the light cone (L ≳ 4τ + 2), the τ = 1 result is the more complete of the two: 2-wide strips
already converge in length, and widths 2–4 agree within 5%. At τ = 2 a converged two-dimensional
value would need patches near 10×10, beyond exact methods, and remains a transferred assumption with
the width trend above.
