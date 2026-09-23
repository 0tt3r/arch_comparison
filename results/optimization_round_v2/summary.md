# Optimization round v2: calibrated surface-code law

Identical to [round v1](../optimization_round_v1/summary.md) except that both encoded arms use the
Stim-calibrated surface-code law of [SURFACE_LAW.md](../../docs/SURFACE_LAW.md) (per-p fit; Λ = 7.41 at p = 10⁻³,
89 at p = 10⁻⁴; conservative X+Z amplitude). Physical-qubit rows replay round v1 exactly (2184 checks).
Improvement factors below are measured outputs of the declared model, not assumptions.

[Headline figure (τ = L/4)](../../figures/optimization_round_capacity_v2.pdf) · [Fixed τ = 1 benchmark](../../figures/optimization_round_fixed_tau_v2.pdf) · [Leader map](../../figures/optimization_round_leaders_v2.pdf) · [Methods](../../docs/OPTIMIZATION_ROUND.md)

## Encoded capacity leader by configuration

Entries give the leader and its lattice side at Q = 10⁶ and 10⁸. The time-growing family (τ = L/4) is the headline;
fixed τ = 1 is retained as an implementation benchmark.

### τ = 1

| p | ε | Q | v1 menu | + raw RUS | + lanes | + FT delivery | + separate distances | + Richardson | all five |
|---:|---:|---:|---|---|---|---|---|---|---|
| 0.001 | 0.01 | 10^6 | ft (10/12/14) | tie:hybrid+ft (10/14/14) | ft (10/12/14) | ft (10/12/14) | ft (10/12/14) | ft (10/12/14) | tie:hybrid+ft (10/14/14) |
| 0.001 | 0.01 | 10^8 | ft (12/16/48) | ft (14/20/48) | ft (12/18/48) | ft (12/16/48) | ft (12/16/48) | ft (12/16/48) | ft (16/24/48) |
| 0.001 | 0.001 | 10^6 | hybrid (4/6/4) | hybrid (4/6/4) | hybrid (4/6/4) | ft (4/6/8) | hybrid (4/6/4) | tie:hybrid+ft (4/6/6) | ft (6/8/10) |
| 0.001 | 0.001 | 10^8 | ft (6/10/16) | ft (8/10/16) | ft (8/10/16) | ft (6/10/44) | ft (6/10/16) | ft (8/10/18) | ft (8/14/48) |
| 0.0001 | 0.01 | 10^6 | ft (16/24/26) | ft (18/24/26) | ft (18/24/26) | ft (16/24/26) | ft (16/26/28) | ft (16/24/26) | ft (18/26/30) |
| 0.0001 | 0.01 | 10^8 | ft (22/36/48) | ft (28/38/48) | ft (24/38/48) | ft (22/38/48) | ft (22/36/48) | ft (22/36/48) | ft (32/46/48) |
| 0.0001 | 0.001 | 10^6 | tie:star+hybrid (8/8/6) | tie:star+hybrid (8/8/6) | hybrid (8/10/6) | tie:hybrid+ft (8/10/10) | hybrid (8/10/6) | hybrid (8/12/8) | tie:hybrid+ft (8/16/16) |
| 0.0001 | 0.001 | 10^8 | ft (12/16/34) | ft (14/18/34) | ft (14/20/34) | ft (12/16/48) | ft (12/16/34) | ft (12/20/38) | ft (18/26/48) |

Cells: leader (STAR/hybrid/FT sides).

### τ = L/4

| p | ε | Q | v1 menu | + raw RUS | + lanes | + FT delivery | + separate distances | + Richardson | all five |
|---:|---:|---:|---|---|---|---|---|---|---|
| 0.001 | 0.01 | 10^6 | ft (6/8/12) | ft (8/10/12) | ft (6/8/12) | ft (6/8/12) | ft (6/8/12) | ft (6/8/12) | ft (8/10/12) |
| 0.001 | 0.01 | 10^8 | ft (8/10/32) | ft (10/12/32) | ft (8/10/32) | ft (8/10/32) | ft (8/10/32) | ft (8/10/32) | ft (10/12/32) |
| 0.001 | 0.001 | 10^6 | tie:star+hybrid+ft (4/4/4) | tie:star+hybrid+ft (4/4/4) | tie:star+hybrid+ft (4/4/4) | ft (4/4/6) | tie:star+hybrid+ft (4/4/4) | hybrid (4/6/4) | ft (4/6/8) |
| 0.001 | 0.001 | 10^8 | ft (6/6/10) | ft (6/8/10) | ft (6/6/10) | ft (6/6/18) | ft (6/6/10) | ft (6/8/12) | ft (6/8/20) |
| 0.0001 | 0.01 | 10^6 | ft (10/14/18) | ft (10/14/18) | ft (10/14/18) | ft (10/14/24) | ft (10/14/20) | ft (10/14/18) | ft (12/16/24) |
| 0.0001 | 0.01 | 10^8 | ft (12/18/32) | ft (16/18/32) | ft (14/18/32) | ft (12/18/32) | ft (12/18/32) | ft (12/18/32) | ft (16/20/32) |
| 0.0001 | 0.001 | 10^6 | tie:star+hybrid (6/6/4) | tie:star+hybrid (6/6/4) | tie:star+hybrid (6/6/4) | tie:star+hybrid+ft (6/6/6) | tie:star+hybrid (6/6/4) | hybrid (6/8/6) | hybrid (6/10/8) |
| 0.0001 | 0.001 | 10^8 | ft (8/10/18) | ft (10/12/18) | ft (10/12/18) | ft (8/10/32) | ft (8/10/18) | ft (8/12/20) | ft (10/14/32) |

Cells: leader (STAR/hybrid/FT sides).

## How often each change moves the leader

Out of 64 (family, p, ε, Q) cells, relative to the v1 menu:

| Configuration | Cells with a different encoded leader |
|---|---:|
| + raw RUS | 2 |
| + lanes | 3 |
| + FT delivery | 8 |
| + separate distances | 4 |
| + Richardson | 6 |
| all five | 15 |

## What each change buys at representative points

log₁₀ Q_min for one week; — means no feasible design within the declared menu.

| Family | p | ε | L | Arm | v1 menu | + raw RUS | + lanes | + FT delivery | + separate distances | + Richardson | all five |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| τ=1 | 0.001 | 0.01 | 8 | star | 5.45 | 5.38 | 5.42 | 5.45 | 5.45 | 5.45 | 5.38 |
| τ=1 | 0.001 | 0.01 | 8 | hybrid | 5.39 | 5.39 | 5.39 | 5.39 | 5.39 | 5.39 | 5.36 |
| τ=1 | 0.001 | 0.01 | 8 | ft | 5.48 | 5.48 | 5.48 | 5.43 | 5.48 | 5.48 | 5.43 |
| τ=1 | 0.001 | 0.01 | 8 | nisq_effective | 8.29 | 8.29 | 8.29 | 8.29 | 8.29 | 8.29 | 8.29 |
| τ=1 | 0.001 | 0.001 | 8 | star | 8.01 | 7.63 | 7.87 | 8.01 | 8.01 | 7.40 | 6.88 |
| τ=1 | 0.001 | 0.001 | 8 | hybrid | 6.92 | 6.59 | 6.67 | 6.81 | 6.92 | 6.35 | 5.69 |
| τ=1 | 0.001 | 0.001 | 8 | ft | 6.97 | 6.97 | 6.97 | 5.86 | 6.97 | 6.64 | 5.74 |
| τ=1 | 0.001 | 0.001 | 8 | nisq_effective | 29.40 | 29.40 | 29.40 | 29.40 | 29.40 | 22.88 | 22.88 |
| τ=1 | 0.0001 | 0.01 | 8 | star | 4.89 | 4.89 | 4.85 | 4.89 | 4.89 | 4.89 | 4.85 |
| τ=1 | 0.0001 | 0.01 | 8 | hybrid | 5.02 | 5.02 | 4.98 | 5.02 | 4.89 | 5.02 | 4.86 |
| τ=1 | 0.0001 | 0.01 | 8 | ft | 5.02 | 5.02 | 5.02 | 4.97 | 4.75 | 5.02 | 4.70 |
| τ=1 | 0.0001 | 0.01 | 8 | nisq_effective | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 |
| τ=1 | 0.0001 | 0.001 | 8 | star | 5.93 | 5.89 | 5.85 | 5.93 | 5.93 | 5.36 | 5.36 |
| τ=1 | 0.0001 | 0.001 | 8 | hybrid | 5.48 | 5.48 | 5.21 | 5.22 | 5.35 | 5.18 | 5.03 |
| τ=1 | 0.0001 | 0.001 | 8 | ft | 6.26 | 6.26 | 6.26 | 5.35 | 6.21 | 5.84 | 5.10 |
| τ=1 | 0.0001 | 0.001 | 8 | nisq_effective | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 |
| τ=1 | 0.001 | 0.01 | 16 | star | 10.25 | 8.22 | 9.77 | 10.25 | 10.25 | 10.25 | 7.74 |
| τ=1 | 0.001 | 0.01 | 16 | hybrid | 7.25 | 6.22 | 6.94 | 7.13 | 7.25 | 7.25 | 6.17 |
| τ=1 | 0.001 | 0.01 | 16 | ft | 6.14 | 6.14 | 6.14 | 6.08 | 6.14 | 6.14 | 6.08 |
| τ=1 | 0.001 | 0.01 | 16 | nisq_effective | 97.64 | 97.64 | 97.64 | 97.64 | 97.64 | 97.64 | 97.64 |
| τ=1 | 0.001 | 0.001 | 16 | star | 15.49 | 14.39 | 15.03 | 15.49 | 15.49 | 14.29 | 12.28 |
| τ=1 | 0.001 | 0.001 | 16 | hybrid | 11.06 | 10.10 | 10.81 | 10.96 | 11.06 | 10.09 | 8.48 |
| τ=1 | 0.001 | 0.001 | 16 | ft | 7.85 | 7.85 | 7.85 | 7.16 | 7.85 | 7.68 | 6.44 |
| τ=1 | 0.001 | 0.001 | 16 | nisq_effective | 231.85 | 231.85 | 231.85 | 231.85 | 231.85 | 200.73 | 200.73 |
| τ=1 | 0.0001 | 0.01 | 16 | star | 5.43 | 5.43 | 5.43 | 5.43 | 5.43 | 5.43 | 5.43 |
| τ=1 | 0.0001 | 0.01 | 16 | hybrid | 5.56 | 5.56 | 5.55 | 5.56 | 5.44 | 5.56 | 5.43 |
| τ=1 | 0.0001 | 0.01 | 16 | ft | 5.56 | 5.56 | 5.56 | 5.55 | 5.44 | 5.56 | 5.42 |
| τ=1 | 0.0001 | 0.01 | 16 | nisq_effective | 6.62 | 6.62 | 6.62 | 6.62 | 6.62 | 6.62 | 6.62 |
| τ=1 | 0.0001 | 0.001 | 16 | star | 8.51 | 8.03 | 8.08 | 8.51 | 8.51 | 8.51 | 7.61 |
| τ=1 | 0.0001 | 0.001 | 16 | hybrid | 7.73 | 7.65 | 7.32 | 7.58 | 7.72 | 6.86 | 5.78 |
| τ=1 | 0.0001 | 0.001 | 16 | ft | 6.92 | 6.92 | 6.92 | 6.38 | 6.92 | 6.75 | 5.76 |
| τ=1 | 0.0001 | 0.001 | 16 | nisq_effective | 24.26 | 24.26 | 24.26 | 24.26 | 24.26 | 18.69 | 18.69 |
| τ=L/4 | 0.001 | 0.01 | 8 | star | 6.36 | 5.52 | 6.18 | 6.36 | 6.36 | 6.36 | 5.52 |
| τ=L/4 | 0.001 | 0.01 | 8 | hybrid | 5.46 | 5.46 | 5.44 | 5.46 | 5.46 | 5.46 | 5.42 |
| τ=L/4 | 0.001 | 0.01 | 8 | ft | 5.52 | 5.52 | 5.52 | 5.48 | 5.52 | 5.52 | 5.48 |
| τ=L/4 | 0.001 | 0.01 | 8 | nisq_effective | 21.30 | 21.30 | 21.30 | 21.30 | 21.30 | 21.30 | 21.30 |
| τ=L/4 | 0.001 | 0.001 | 8 | star | 10.29 | 9.68 | 10.17 | 10.29 | 10.29 | 9.43 | 8.53 |
| τ=L/4 | 0.001 | 0.001 | 8 | hybrid | 8.51 | 7.97 | 8.32 | 8.36 | 8.51 | 7.79 | 6.35 |
| τ=L/4 | 0.001 | 0.001 | 8 | ft | 7.33 | 7.33 | 7.33 | 6.88 | 7.33 | 7.20 | 5.92 |
| τ=L/4 | 0.001 | 0.001 | 8 | nisq_effective | 60.66 | 60.66 | 60.66 | 60.66 | 60.66 | 48.51 | 48.51 |
| τ=L/4 | 0.0001 | 0.01 | 8 | star | 4.89 | 4.89 | 4.85 | 4.89 | 4.89 | 4.89 | 4.85 |
| τ=L/4 | 0.0001 | 0.01 | 8 | hybrid | 5.02 | 5.02 | 4.98 | 5.02 | 4.89 | 5.02 | 4.86 |
| τ=L/4 | 0.0001 | 0.01 | 8 | ft | 5.02 | 5.02 | 5.02 | 4.97 | 4.77 | 5.02 | 4.71 |
| τ=L/4 | 0.0001 | 0.01 | 8 | nisq_effective | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 |
| τ=L/4 | 0.0001 | 0.001 | 8 | star | 7.12 | 6.89 | 7.00 | 7.12 | 7.12 | 6.36 | 6.19 |
| τ=L/4 | 0.0001 | 0.001 | 8 | hybrid | 6.36 | 6.30 | 6.13 | 6.24 | 6.14 | 5.48 | 5.11 |
| τ=L/4 | 0.0001 | 0.001 | 8 | ft | 6.56 | 6.56 | 6.56 | 6.28 | 6.51 | 6.48 | 5.28 |
| τ=L/4 | 0.0001 | 0.001 | 8 | nisq_effective | 5.74 | 5.74 | 5.74 | 5.74 | 5.74 | 3.06 | 3.06 |
| τ=L/4 | 0.001 | 0.01 | 16 | star | 27.89 | 19.86 | 27.44 | 27.89 | 27.89 | 27.89 | 19.40 |
| τ=L/4 | 0.001 | 0.01 | 16 | hybrid | 16.27 | 11.79 | 15.97 | 16.10 | 16.27 | 16.27 | 11.13 |
| τ=L/4 | 0.001 | 0.01 | 16 | ft | 6.29 | 6.29 | 6.29 | 6.19 | 6.29 | 6.29 | 6.19 |
| τ=L/4 | 0.001 | 0.01 | 16 | nisq_effective | 401.90 | 401.90 | 401.90 | 401.90 | 401.90 | 401.90 | 401.90 |
| τ=L/4 | 0.001 | 0.001 | 16 | star | 38.73 | 38.73 | 38.30 | 38.73 | 38.73 | 37.07 | 36.60 |
| τ=L/4 | 0.001 | 0.001 | 16 | hybrid | 23.10 | 23.10 | 22.78 | 22.93 | 23.10 | 20.63 | 19.99 |
| τ=L/4 | 0.001 | 0.001 | 16 | ft | 8.58 | 8.58 | 8.58 | 7.69 | 8.58 | 8.40 | 7.48 |
| τ=L/4 | 0.001 | 0.001 | 16 | nisq_effective | 1026.38 | 1026.38 | 1026.38 | 1026.38 | 1026.38 | 808.24 | 808.24 |
| τ=L/4 | 0.0001 | 0.01 | 16 | star | 9.75 | 7.84 | 9.32 | 9.75 | 9.75 | 9.75 | 7.43 |
| τ=L/4 | 0.0001 | 0.01 | 16 | hybrid | 6.82 | 6.56 | 6.64 | 6.68 | 6.71 | 6.82 | 5.69 |
| τ=L/4 | 0.0001 | 0.01 | 16 | ft | 5.68 | 5.68 | 5.68 | 5.57 | 5.56 | 5.68 | 5.45 |
| τ=L/4 | 0.0001 | 0.01 | 16 | nisq_effective | 37.81 | 37.81 | 37.81 | 37.81 | 37.81 | 37.81 | 37.81 |
| τ=L/4 | 0.0001 | 0.001 | 16 | star | 12.18 | 10.22 | 11.75 | 12.18 | 12.18 | 12.18 | 9.80 |
| τ=L/4 | 0.0001 | 0.001 | 16 | hybrid | 10.49 | 9.30 | 9.76 | 10.48 | 10.49 | 10.28 | 8.31 |
| τ=L/4 | 0.0001 | 0.001 | 16 | ft | 7.61 | 7.61 | 7.61 | 6.89 | 7.61 | 7.47 | 6.79 |
| τ=L/4 | 0.0001 | 0.001 | 16 | nisq_effective | 104.78 | 104.78 | 104.78 | 104.78 | 104.78 | 80.06 | 80.06 |

## Selected designs with all five changes

What the optimizer chose, which is the mechanism behind each shift.

| Family | p | ε | L | Arm | log₁₀ Q | r | extrapolation | d | RUS | lanes | ports | factories | HWP batch | qubits/copy | s/attempt | χ |
|---|---:|---:|---:|---|---:|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| τ=1 | 0.001 | 0.01 | 8 | star | 5.38 | 8 | single | 23 | raw | 4 | 1 | 0 | — | 241,224 | 0.287 | 1.89 |
| τ=1 | 0.001 | 0.01 | 8 | hybrid | 5.36 | 9 | single | 23 | raw | 2 | 1 | 1 | — | 230,610 | 0.481 | 0.89 |
| τ=1 | 0.001 | 0.01 | 8 | ft | 5.43 | 9 | single | 25 | — | — | 1 | 1 | 8 | 271,024 | 5.52 | 0.00 |
| τ=1 | 0.001 | 0.001 | 8 | star | 6.88 | 8 | richardson | 29 | raw | 8 | 1 | 0 | — | 588,700 | 0.333 | 3.54 |
| τ=1 | 0.001 | 0.001 | 8 | hybrid | 5.69 | 8 | richardson | 25 | raw | 8 | 16 | 8 | — | 493,192 | 0.0778 | 2.00 |
| τ=1 | 0.001 | 0.001 | 8 | ft | 5.74 | 8 | richardson | 25 | — | — | 8 | 32 | 32 | 550,268 | 0.123 | 0.90 |
| τ=1 | 0.0001 | 0.01 | 8 | star | 4.85 | 9 | single | 13 | raw | 1 | 1 | 0 | — | 70,980 | 0.609 | 0.98 |
| τ=1 | 0.0001 | 0.01 | 8 | hybrid | 4.86 | 11 | single | 13 | raw | 1 | 1 | 1 | — | 72,342 | 0.411 | 0.35 |
| τ=1 | 0.0001 | 0.01 | 8 | ft | 4.70 | 13 | single | 11 | — | — | 1 | 1 | 1 | 50,634 | 8.87 | 0.00 |
| τ=1 | 0.0001 | 0.001 | 8 | star | 5.36 | 8 | richardson | 13 | raw | 4 | 1 | 0 | — | 77,064 | 0.269 | 2.49 |
| τ=1 | 0.0001 | 0.001 | 8 | hybrid | 5.03 | 8 | richardson | 13 | raw | 1 | 2 | 2 | — | 107,504 | 0.315 | 1.36 |
| τ=1 | 0.0001 | 0.001 | 8 | ft | 5.10 | 8 | richardson | 11 | — | — | 2 | 16 | 64 | 126,010 | 0.15 | 0.90 |
| τ=L/4 | 0.001 | 0.01 | 8 | star | 5.52 | 16 | single | 27 | raw | 4 | 1 | 0 | — | 332,424 | 0.631 | 3.78 |
| τ=L/4 | 0.001 | 0.01 | 8 | hybrid | 5.42 | 23 | single | 25 | raw | 1 | 1 | 1 | — | 264,774 | 1.67 | 1.94 |
| τ=L/4 | 0.001 | 0.01 | 8 | ft | 5.48 | 26 | single | 27 | — | — | 1 | 1 | 1 | 299,914 | 34.1 | 0.00 |
| τ=L/4 | 0.001 | 0.001 | 8 | star | 8.53 | 16 | richardson | 31 | raw | 8 | 1 | 0 | — | 672,700 | 0.712 | 6.21 |
| τ=L/4 | 0.001 | 0.001 | 8 | hybrid | 6.35 | 16 | richardson | 27 | raw | 64 | 16 | 32 | — | 1,123,352 | 0.0303 | 3.10 |
| τ=L/4 | 0.001 | 0.001 | 8 | ft | 5.92 | 16 | richardson | 27 | — | — | 16 | 64 | 64 | 836,818 | 0.126 | 0.90 |
| τ=L/4 | 0.0001 | 0.01 | 8 | star | 4.85 | 16 | single | 13 | raw | 1 | 1 | 0 | — | 70,980 | 1.08 | 1.77 |
| τ=L/4 | 0.0001 | 0.01 | 8 | hybrid | 4.86 | 22 | single | 13 | raw | 1 | 1 | 1 | — | 72,342 | 0.817 | 0.69 |
| τ=L/4 | 0.0001 | 0.01 | 8 | ft | 4.71 | 17 | single | 11 | — | — | 1 | 1 | 4 | 51,844 | 9.2 | 0.00 |
| τ=L/4 | 0.0001 | 0.001 | 8 | star | 6.19 | 16 | richardson | 13 | raw | 8 | 1 | 0 | — | 118,300 | 0.269 | 4.12 |
| τ=L/4 | 0.0001 | 0.001 | 8 | hybrid | 5.11 | 16 | richardson | 13 | raw | 8 | 1 | 8 | — | 129,196 | 0.144 | 1.81 |
| τ=L/4 | 0.0001 | 0.001 | 8 | ft | 5.28 | 16 | richardson | 11 | — | — | 8 | 64 | 64 | 189,682 | 0.0737 | 0.90 |
| τ=L/4 | 0.001 | 0.01 | 16 | star | 19.40 | 32 | single | 45 | raw | 16 | 1 | 0 | — | 4,851,900 | 2.24 | 32.07 |
| τ=L/4 | 0.001 | 0.01 | 16 | hybrid | 11.13 | 40 | single | 35 | raw | 128 | 128 | 1024 | — | 38,866,776 | 0.0785 | 14.37 |
| τ=L/4 | 0.001 | 0.01 | 16 | ft | 6.19 | 38 | single | 25 | — | — | 8 | 32 | 32 | 1,547,768 | 1.02 | 0.00 |
| τ=L/4 | 0.001 | 0.001 | 16 | star | 36.60 | 32 | richardson | 65 | noisy_pcec | 16 | 1 | 0 | — | 10,123,100 | 6.35 | 66.28 |
| τ=L/4 | 0.001 | 0.001 | 16 | hybrid | 19.99 | 32 | richardson | 45 | noisy_pcec | 128 | 128 | 1024 | — | 63,564,376 | 0.155 | 29.91 |
| τ=L/4 | 0.001 | 0.001 | 16 | ft | 7.48 | 32 | richardson | 31 | — | — | 32 | 1024 | 256 | 29,878,576 | 0.464 | 0.90 |
| τ=L/4 | 0.0001 | 0.01 | 16 | star | 7.43 | 41 | single | 25 | raw | 16 | 1 | 0 | — | 1,497,500 | 16.8 | 4.43 |
| τ=L/4 | 0.0001 | 0.01 | 16 | hybrid | 5.69 | 34 | single | 13 | raw | 32 | 8 | 32 | — | 485,688 | 0.117 | 4.32 |
| τ=L/4 | 0.0001 | 0.01 | 16 | ft | 5.45 | 47 | single | 13 | — | — | 1 | 2 | 32 | 282,250 | 33.1 | 0.00 |
| τ=L/4 | 0.0001 | 0.001 | 16 | star | 9.80 | 186 | single | 25 | raw | 16 | 1 | 0 | — | 1,497,500 | 34.8 | 4.69 |
| τ=L/4 | 0.0001 | 0.001 | 16 | hybrid | 8.31 | 41 | richardson | 25 | raw | 128 | 64 | 128 | — | 4,706,072 | 1.2 | 3.67 |
| τ=L/4 | 0.0001 | 0.001 | 16 | ft | 6.79 | 32 | richardson | 13 | — | — | 16 | 1024 | 256 | 6,107,760 | 0.343 | 0.90 |

## Scope

Analytical capacity model with the integrated_v1 component instruments. Concurrency widths, per-lane and per-port tile charges, and the source-to-data transfer are declared assumptions, not compiled schedules. Richardson residuals are the measured finite-patch envelope of |(4E(2r)-E(r))/3| r^4 at r >= 8 tau. The time-growing family uses the running maximum of all measured patches through each time; no patch larger than 4x4 is measured. Interacting bulk convergence is not established. v2 uses the Stim-calibrated surface-code law for both encoded arms; its fit is memory-only at d=3-9 and is extrapolated to the distances used here.

Raw-RUS gains rest on first-order coherent cancellation in the restricted instrument; Toshio et al. App. E finds raw RUS slightly worse than PCEC. See [RAW_RUS_NOTE.md](../../docs/RAW_RUS_NOTE.md).

Both encoded arms use the Stim-calibrated surface-code law of [SURFACE_LAW.md](../../docs/SURFACE_LAW.md), fitted at d = 3–9 and extrapolated to the distances used here.

Validation: PASS: round-v1 NISQ replay, calibrated surface law, budgets, fault ledgers, switch domains, PCEC >= raw. Reproduce with `make optimization-round` (compute node recommended; `pbs/optimization_round.pbs`), redraw with `make optimization-round-figures`, audit with `make optimization-round-check`.
