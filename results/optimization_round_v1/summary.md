# Optimization round: which choices change the winner

This round implements the five changes agreed in [OPUS55_REVIEW.md](../../OPUS55_REVIEW.md) §0,
each alone and all together, against a `v1_menu` configuration that reproduces
`results/integrated_v1` and the NISQ rows of `results/physics_sweep_v1` (335 replays checked).
Improvement factors below are measured outputs of the declared model, not assumptions.

[Headline figure (τ = L/4)](../../figures/optimization_round_capacity.pdf) · [Fixed τ = 1 benchmark](../../figures/optimization_round_fixed_tau.pdf) · [Leader map](../../figures/optimization_round_leaders.pdf) · [Methods](../../docs/OPTIMIZATION_ROUND.md)

## Encoded capacity leader by configuration

Entries give the leader and its lattice side at Q = 10⁶ and 10⁸. The time-growing family (τ = L/4) is the headline;
fixed τ = 1 is retained as an implementation benchmark.

### τ = 1

| p | ε | Q | v1 menu | + raw RUS | + lanes | + FT delivery | + separate distances | + Richardson | all five |
|---:|---:|---:|---|---|---|---|---|---|---|
| 0.001 | 0.01 | 10^6 | ft (10/12/14) | tie:hybrid+ft (10/14/14) | tie:hybrid+ft (10/14/14) | ft (10/12/16) | ft (10/12/14) | ft (10/12/14) | ft (12/14/16) |
| 0.001 | 0.01 | 10^8 | ft (12/18/48) | ft (14/22/48) | ft (12/18/48) | ft (12/18/48) | ft (12/18/48) | ft (12/18/48) | ft (16/24/48) |
| 0.001 | 0.001 | 10^6 | hybrid (4/6/4) | hybrid (4/6/4) | hybrid (4/6/4) | ft (4/6/8) | hybrid (4/6/4) | tie:star+hybrid+ft (6/6/6) | tie:hybrid+ft (6/10/10) |
| 0.001 | 0.001 | 10^8 | ft (8/10/18) | ft (8/10/18) | ft (8/10/18) | ft (8/10/48) | ft (8/10/18) | ft (8/10/22) | ft (10/14/48) |
| 0.0001 | 0.01 | 10^6 | ft (16/24/26) | ft (18/24/26) | ft (18/24/26) | ft (16/24/26) | ft (16/26/28) | ft (16/24/26) | ft (18/26/30) |
| 0.0001 | 0.01 | 10^8 | ft (22/36/48) | ft (28/38/48) | ft (24/38/48) | ft (22/38/48) | ft (22/36/48) | ft (22/36/48) | ft (32/46/48) |
| 0.0001 | 0.001 | 10^6 | tie:star+hybrid (8/8/6) | tie:star+hybrid (8/8/6) | hybrid (8/10/6) | tie:hybrid+ft (8/10/10) | hybrid (8/10/6) | hybrid (8/12/8) | ft (8/16/18) |
| 0.0001 | 0.001 | 10^8 | ft (12/16/34) | ft (14/18/34) | ft (14/20/34) | ft (12/16/48) | ft (12/16/34) | ft (12/20/38) | ft (18/26/48) |

Cells: leader (STAR/hybrid/FT sides).

### τ = L/4

| p | ε | Q | v1 menu | + raw RUS | + lanes | + FT delivery | + separate distances | + Richardson | all five |
|---:|---:|---:|---|---|---|---|---|---|---|
| 0.001 | 0.01 | 10^6 | ft (6/8/12) | ft (8/10/12) | ft (8/8/12) | ft (6/8/14) | ft (6/8/12) | ft (6/8/12) | ft (8/10/14) |
| 0.001 | 0.01 | 10^8 | ft (8/10/32) | ft (10/12/32) | ft (8/10/32) | ft (8/10/32) | ft (8/10/32) | ft (8/10/32) | ft (10/12/32) |
| 0.001 | 0.001 | 10^6 | tie:star+hybrid+ft (4/4/4) | tie:star+hybrid+ft (4/4/4) | hybrid (4/6/4) | ft (4/4/6) | tie:star+hybrid+ft (4/4/4) | hybrid (4/6/4) | ft (4/6/8) |
| 0.001 | 0.001 | 10^8 | ft (6/6/12) | ft (6/8/12) | ft (6/6/12) | ft (6/6/20) | ft (6/6/12) | ft (6/8/12) | ft (6/8/22) |
| 0.0001 | 0.01 | 10^6 | ft (10/14/18) | ft (10/14/18) | ft (10/14/18) | ft (10/14/24) | ft (10/14/20) | ft (10/14/18) | ft (12/16/26) |
| 0.0001 | 0.01 | 10^8 | ft (12/18/32) | ft (16/18/32) | ft (14/18/32) | ft (12/18/32) | ft (12/18/32) | ft (12/18/32) | ft (16/20/32) |
| 0.0001 | 0.001 | 10^6 | tie:star+hybrid (6/6/4) | tie:star+hybrid (6/6/4) | tie:star+hybrid (6/6/4) | tie:star+hybrid+ft (6/6/6) | tie:star+hybrid (6/6/4) | hybrid (6/8/6) | hybrid (6/10/8) |
| 0.0001 | 0.001 | 10^8 | ft (8/10/18) | ft (10/12/18) | ft (10/12/18) | ft (8/10/32) | ft (8/10/18) | ft (8/12/20) | ft (10/14/32) |

Cells: leader (STAR/hybrid/FT sides).

## How often each change moves the leader

Out of 64 (family, p, ε, Q) cells, relative to the v1 menu:

| Configuration | Cells with a different encoded leader |
|---|---:|
| + raw RUS | 1 |
| + lanes | 6 |
| + FT delivery | 8 |
| + separate distances | 4 |
| + Richardson | 6 |
| all five | 14 |

## What each change buys at representative points

log₁₀ Q_min for one week; — means no feasible design within the declared menu.

| Family | p | ε | L | Arm | v1 menu | + raw RUS | + lanes | + FT delivery | + separate distances | + Richardson | all five |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| τ=1 | 0.001 | 0.01 | 8 | star | 5.30 | 5.30 | 5.30 | 5.30 | 5.30 | 5.30 | 5.28 |
| τ=1 | 0.001 | 0.01 | 8 | hybrid | 5.31 | 5.31 | 5.27 | 5.31 | 5.31 | 5.31 | 5.27 |
| τ=1 | 0.001 | 0.01 | 8 | ft | 5.39 | 5.39 | 5.39 | 5.34 | 5.39 | 5.39 | 5.34 |
| τ=1 | 0.001 | 0.01 | 8 | nisq_effective | 8.29 | 8.29 | 8.29 | 8.29 | 8.29 | 8.29 | 8.29 |
| τ=1 | 0.001 | 0.001 | 8 | star | 7.84 | 7.40 | 7.65 | 7.84 | 7.84 | 7.17 | 6.64 |
| τ=1 | 0.001 | 0.001 | 8 | hybrid | 6.69 | 6.39 | 6.48 | 6.58 | 6.69 | 6.09 | 5.58 |
| τ=1 | 0.001 | 0.001 | 8 | ft | 6.86 | 6.86 | 6.86 | 5.73 | 6.86 | 6.51 | 5.66 |
| τ=1 | 0.001 | 0.001 | 8 | nisq_effective | 29.40 | 29.40 | 29.40 | 29.40 | 29.40 | 22.88 | 22.88 |
| τ=1 | 0.0001 | 0.01 | 8 | star | 4.89 | 4.89 | 4.85 | 4.89 | 4.89 | 4.89 | 4.85 |
| τ=1 | 0.0001 | 0.01 | 8 | hybrid | 5.02 | 5.02 | 4.98 | 5.02 | 4.89 | 5.02 | 4.86 |
| τ=1 | 0.0001 | 0.01 | 8 | ft | 5.02 | 5.02 | 5.02 | 4.97 | 4.75 | 5.02 | 4.70 |
| τ=1 | 0.0001 | 0.01 | 8 | nisq_effective | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 |
| τ=1 | 0.0001 | 0.001 | 8 | star | 5.93 | 5.89 | 5.85 | 5.93 | 5.93 | 5.36 | 5.36 |
| τ=1 | 0.0001 | 0.001 | 8 | hybrid | 5.48 | 5.48 | 5.21 | 5.22 | 5.35 | 5.18 | 5.03 |
| τ=1 | 0.0001 | 0.001 | 8 | ft | 6.26 | 6.26 | 6.26 | 5.35 | 6.15 | 5.84 | 5.04 |
| τ=1 | 0.0001 | 0.001 | 8 | nisq_effective | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 |
| τ=1 | 0.001 | 0.01 | 16 | star | 10.08 | 8.04 | 9.61 | 10.08 | 10.08 | 10.08 | 7.56 |
| τ=1 | 0.001 | 0.01 | 16 | hybrid | 7.06 | 6.08 | 6.75 | 6.93 | 7.06 | 7.06 | 6.07 |
| τ=1 | 0.001 | 0.01 | 16 | ft | 6.01 | 6.01 | 6.01 | 6.00 | 6.01 | 6.01 | 6.00 |
| τ=1 | 0.001 | 0.01 | 16 | nisq_effective | 97.64 | 97.64 | 97.64 | 97.64 | 97.64 | 97.64 | 97.64 |
| τ=1 | 0.001 | 0.001 | 16 | star | 15.32 | 14.25 | 14.85 | 15.32 | 15.32 | 14.12 | 12.11 |
| τ=1 | 0.001 | 0.001 | 16 | hybrid | 10.90 | 9.93 | 10.65 | 10.84 | 10.90 | 9.87 | 8.32 |
| τ=1 | 0.001 | 0.001 | 16 | ft | 7.71 | 7.71 | 7.71 | 6.88 | 7.71 | 7.53 | 6.31 |
| τ=1 | 0.001 | 0.001 | 16 | nisq_effective | 231.85 | 231.85 | 231.85 | 231.85 | 231.85 | 200.73 | 200.73 |
| τ=1 | 0.0001 | 0.01 | 16 | star | 5.43 | 5.43 | 5.43 | 5.43 | 5.43 | 5.43 | 5.43 |
| τ=1 | 0.0001 | 0.01 | 16 | hybrid | 5.56 | 5.56 | 5.55 | 5.56 | 5.44 | 5.56 | 5.43 |
| τ=1 | 0.0001 | 0.01 | 16 | ft | 5.56 | 5.56 | 5.56 | 5.55 | 5.44 | 5.56 | 5.30 |
| τ=1 | 0.0001 | 0.01 | 16 | nisq_effective | 6.62 | 6.62 | 6.62 | 6.62 | 6.62 | 6.62 | 6.62 |
| τ=1 | 0.0001 | 0.001 | 16 | star | 8.51 | 8.03 | 8.08 | 8.51 | 8.51 | 8.51 | 7.61 |
| τ=1 | 0.0001 | 0.001 | 16 | hybrid | 7.73 | 7.65 | 7.32 | 7.58 | 7.62 | 6.86 | 5.78 |
| τ=1 | 0.0001 | 0.001 | 16 | ft | 6.92 | 6.92 | 6.92 | 6.38 | 6.92 | 6.75 | 5.75 |
| τ=1 | 0.0001 | 0.001 | 16 | nisq_effective | 24.26 | 24.26 | 24.26 | 24.26 | 24.26 | 18.69 | 18.69 |
| τ=L/4 | 0.001 | 0.01 | 8 | star | 6.06 | 5.38 | 5.94 | 6.06 | 6.06 | 6.06 | 5.38 |
| τ=L/4 | 0.001 | 0.01 | 8 | hybrid | 5.39 | 5.31 | 5.35 | 5.39 | 5.39 | 5.39 | 5.28 |
| τ=L/4 | 0.001 | 0.01 | 8 | ft | 5.40 | 5.40 | 5.40 | 5.34 | 5.40 | 5.40 | 5.34 |
| τ=L/4 | 0.001 | 0.01 | 8 | nisq_effective | 21.30 | 21.30 | 21.30 | 21.30 | 21.30 | 21.30 | 21.30 |
| τ=L/4 | 0.001 | 0.001 | 8 | star | 10.12 | 9.47 | 10.01 | 10.12 | 10.12 | 9.25 | 8.35 |
| τ=L/4 | 0.001 | 0.001 | 8 | hybrid | 8.33 | 7.80 | 8.17 | 8.22 | 8.33 | 7.62 | 6.06 |
| τ=L/4 | 0.001 | 0.001 | 8 | ft | 7.21 | 7.21 | 7.21 | 6.76 | 7.21 | 7.03 | 5.80 |
| τ=L/4 | 0.001 | 0.001 | 8 | nisq_effective | 60.66 | 60.66 | 60.66 | 60.66 | 60.66 | 48.51 | 48.51 |
| τ=L/4 | 0.0001 | 0.01 | 8 | star | 4.89 | 4.89 | 4.85 | 4.89 | 4.89 | 4.89 | 4.85 |
| τ=L/4 | 0.0001 | 0.01 | 8 | hybrid | 5.02 | 5.02 | 4.98 | 5.02 | 4.89 | 5.02 | 4.86 |
| τ=L/4 | 0.0001 | 0.01 | 8 | ft | 5.02 | 5.02 | 5.02 | 4.97 | 4.75 | 5.02 | 4.70 |
| τ=L/4 | 0.0001 | 0.01 | 8 | nisq_effective | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 |
| τ=L/4 | 0.0001 | 0.001 | 8 | star | 7.12 | 6.89 | 7.00 | 7.12 | 7.12 | 6.36 | 6.19 |
| τ=L/4 | 0.0001 | 0.001 | 8 | hybrid | 6.36 | 6.30 | 6.13 | 6.24 | 6.14 | 5.48 | 5.11 |
| τ=L/4 | 0.0001 | 0.001 | 8 | ft | 6.56 | 6.56 | 6.56 | 6.28 | 6.51 | 6.48 | 5.28 |
| τ=L/4 | 0.0001 | 0.001 | 8 | nisq_effective | 5.74 | 5.74 | 5.74 | 5.74 | 5.74 | 3.06 | 3.06 |
| τ=L/4 | 0.001 | 0.01 | 16 | star | 27.71 | 19.68 | 27.25 | 27.71 | 27.71 | 27.71 | 19.22 |
| τ=L/4 | 0.001 | 0.01 | 16 | hybrid | 16.08 | 11.60 | 15.81 | 15.97 | 16.08 | 16.08 | 11.01 |
| τ=L/4 | 0.001 | 0.01 | 16 | ft | 6.17 | 6.17 | 6.17 | 6.09 | 6.17 | 6.17 | 6.09 |
| τ=L/4 | 0.001 | 0.01 | 16 | nisq_effective | 401.90 | 401.90 | 401.90 | 401.90 | 401.90 | 401.90 | 401.90 |
| τ=L/4 | 0.001 | 0.001 | 16 | star | 38.56 | 38.56 | 38.13 | 38.56 | 38.56 | 36.89 | 36.43 |
| τ=L/4 | 0.001 | 0.001 | 16 | hybrid | 22.90 | 22.90 | 22.62 | 22.74 | 22.90 | 20.45 | 19.81 |
| τ=L/4 | 0.001 | 0.001 | 16 | ft | 8.44 | 8.44 | 8.44 | 7.58 | 8.44 | 8.27 | 7.36 |
| τ=L/4 | 0.001 | 0.001 | 16 | nisq_effective | 1026.38 | 1026.38 | 1026.38 | 1026.38 | 1026.38 | 808.24 | 808.24 |
| τ=L/4 | 0.0001 | 0.01 | 16 | star | 9.75 | 7.84 | 9.32 | 9.75 | 9.75 | 9.75 | 7.43 |
| τ=L/4 | 0.0001 | 0.01 | 16 | hybrid | 6.82 | 6.56 | 6.64 | 6.68 | 6.64 | 6.82 | 5.69 |
| τ=L/4 | 0.0001 | 0.01 | 16 | ft | 5.68 | 5.68 | 5.68 | 5.57 | 5.56 | 5.68 | 5.45 |
| τ=L/4 | 0.0001 | 0.01 | 16 | nisq_effective | 37.81 | 37.81 | 37.81 | 37.81 | 37.81 | 37.81 | 37.81 |
| τ=L/4 | 0.0001 | 0.001 | 16 | star | 12.18 | 10.22 | 11.75 | 12.18 | 12.18 | 12.18 | 9.80 |
| τ=L/4 | 0.0001 | 0.001 | 16 | hybrid | 10.49 | 9.30 | 9.76 | 10.48 | 10.49 | 10.19 | 8.31 |
| τ=L/4 | 0.0001 | 0.001 | 16 | ft | 7.59 | 7.59 | 7.59 | 6.89 | 7.59 | 7.45 | 6.76 |
| τ=L/4 | 0.0001 | 0.001 | 16 | nisq_effective | 104.78 | 104.78 | 104.78 | 104.78 | 104.78 | 80.06 | 80.06 |

## Selected designs with all five changes

What the optimizer chose, which is the mechanism behind each shift.

| Family | p | ε | L | Arm | log₁₀ Q | r | extrapolation | d | RUS | lanes | ports | factories | HWP batch | qubits/copy | s/attempt | χ |
|---|---:|---:|---:|---|---:|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| τ=1 | 0.001 | 0.01 | 8 | star | 5.28 | 10 | single | 21 | raw | 2 | 1 | 0 | — | 190,512 | 0.646 | 2.06 |
| τ=1 | 0.001 | 0.01 | 8 | hybrid | 5.27 | 10 | single | 21 | raw | 1 | 1 | 1 | — | 187,126 | 0.671 | 0.91 |
| τ=1 | 0.001 | 0.01 | 8 | ft | 5.34 | 10 | single | 23 | — | — | 1 | 1 | 1 | 217,914 | 10.8 | 0.00 |
| τ=1 | 0.001 | 0.001 | 8 | star | 6.64 | 8 | richardson | 25 | raw | 8 | 1 | 0 | — | 437,500 | 0.263 | 3.54 |
| τ=1 | 0.001 | 0.001 | 8 | hybrid | 5.58 | 8 | richardson | 23 | raw | 8 | 2 | 4 | — | 380,744 | 0.109 | 2.00 |
| τ=1 | 0.001 | 0.001 | 8 | ft | 5.66 | 8 | richardson | 23 | — | — | 4 | 32 | 32 | 462,316 | 0.172 | 0.90 |
| τ=1 | 0.0001 | 0.01 | 8 | star | 4.85 | 9 | single | 13 | raw | 1 | 1 | 0 | — | 70,980 | 0.609 | 0.98 |
| τ=1 | 0.0001 | 0.01 | 8 | hybrid | 4.86 | 11 | single | 13 | raw | 1 | 1 | 1 | — | 72,342 | 0.411 | 0.35 |
| τ=1 | 0.0001 | 0.01 | 8 | ft | 4.70 | 13 | single | 11 | — | — | 1 | 1 | 1 | 50,634 | 8.87 | 0.00 |
| τ=1 | 0.0001 | 0.001 | 8 | star | 5.36 | 8 | richardson | 13 | raw | 4 | 1 | 0 | — | 77,064 | 0.269 | 2.49 |
| τ=1 | 0.0001 | 0.001 | 8 | hybrid | 5.03 | 8 | richardson | 13 | raw | 1 | 2 | 2 | — | 107,504 | 0.315 | 1.36 |
| τ=1 | 0.0001 | 0.001 | 8 | ft | 5.04 | 8 | richardson | 11 | — | — | 2 | 16 | 32 | 109,796 | 0.21 | 0.90 |
| τ=L/4 | 0.001 | 0.01 | 8 | star | 5.38 | 16 | single | 23 | raw | 4 | 1 | 0 | — | 241,224 | 0.575 | 3.78 |
| τ=L/4 | 0.001 | 0.01 | 8 | hybrid | 5.28 | 18 | single | 21 | raw | 2 | 1 | 1 | — | 192,418 | 0.874 | 1.77 |
| τ=L/4 | 0.001 | 0.01 | 8 | ft | 5.34 | 21 | single | 23 | — | — | 1 | 1 | 1 | 217,914 | 23.3 | 0.00 |
| τ=L/4 | 0.001 | 0.001 | 8 | star | 8.35 | 16 | richardson | 27 | raw | 8 | 1 | 0 | — | 510,300 | 0.62 | 6.21 |
| τ=L/4 | 0.001 | 0.001 | 8 | hybrid | 6.06 | 16 | richardson | 25 | raw | 64 | 64 | 64 | — | 1,160,536 | 0.0203 | 3.10 |
| τ=L/4 | 0.001 | 0.001 | 8 | ft | 5.80 | 16 | richardson | 23 | — | — | 16 | 64 | 64 | 625,218 | 0.108 | 0.90 |
| τ=L/4 | 0.0001 | 0.01 | 8 | star | 4.85 | 16 | single | 13 | raw | 1 | 1 | 0 | — | 70,980 | 1.08 | 1.77 |
| τ=L/4 | 0.0001 | 0.01 | 8 | hybrid | 4.86 | 22 | single | 13 | raw | 1 | 1 | 1 | — | 72,342 | 0.817 | 0.69 |
| τ=L/4 | 0.0001 | 0.01 | 8 | ft | 4.70 | 26 | single | 11 | — | — | 1 | 1 | 1 | 50,634 | 18.3 | 0.00 |
| τ=L/4 | 0.0001 | 0.001 | 8 | star | 6.19 | 16 | richardson | 13 | raw | 8 | 1 | 0 | — | 118,300 | 0.269 | 4.12 |
| τ=L/4 | 0.0001 | 0.001 | 8 | hybrid | 5.11 | 16 | richardson | 13 | raw | 8 | 1 | 8 | — | 129,196 | 0.144 | 1.81 |
| τ=L/4 | 0.0001 | 0.001 | 8 | ft | 5.28 | 16 | richardson | 11 | — | — | 8 | 64 | 64 | 189,682 | 0.0737 | 0.90 |
| τ=L/4 | 0.001 | 0.01 | 16 | star | 19.22 | 32 | single | 39 | raw | 16 | 1 | 0 | — | 3,644,316 | 1.94 | 32.07 |
| τ=L/4 | 0.001 | 0.01 | 16 | hybrid | 11.01 | 33 | single | 31 | raw | 128 | 128 | 1024 | — | 30,716,568 | 0.0614 | 13.77 |
| τ=L/4 | 0.001 | 0.01 | 16 | ft | 6.09 | 39 | single | 21 | — | — | 16 | 64 | 64 | 1,235,068 | 0.361 | 0.00 |
| τ=L/4 | 0.001 | 0.001 | 16 | star | 36.43 | 32 | richardson | 57 | noisy_pcec | 16 | 1 | 0 | — | 7,784,604 | 5.57 | 66.28 |
| τ=L/4 | 0.001 | 0.001 | 16 | hybrid | 19.81 | 32 | richardson | 39 | noisy_pcec | 128 | 128 | 1024 | — | 48,004,888 | 0.135 | 29.91 |
| τ=L/4 | 0.001 | 0.001 | 16 | ft | 7.36 | 32 | richardson | 27 | — | — | 32 | 1024 | 256 | 22,918,576 | 0.437 | 0.90 |
| τ=L/4 | 0.0001 | 0.01 | 16 | star | 7.43 | 41 | single | 25 | raw | 16 | 1 | 0 | — | 1,497,500 | 16.8 | 4.43 |
| τ=L/4 | 0.0001 | 0.01 | 16 | hybrid | 5.69 | 34 | single | 13 | raw | 32 | 8 | 32 | — | 485,688 | 0.117 | 4.32 |
| τ=L/4 | 0.0001 | 0.01 | 16 | ft | 5.45 | 47 | single | 13 | — | — | 1 | 1 | 32 | 280,888 | 36 | 0.00 |
| τ=L/4 | 0.0001 | 0.001 | 16 | star | 9.80 | 186 | single | 25 | raw | 16 | 1 | 0 | — | 1,497,500 | 34.8 | 4.69 |
| τ=L/4 | 0.0001 | 0.001 | 16 | hybrid | 8.31 | 41 | richardson | 25 | raw | 128 | 64 | 128 | — | 4,706,072 | 1.2 | 3.67 |
| τ=L/4 | 0.0001 | 0.001 | 16 | ft | 6.76 | 32 | richardson | 13 | — | — | 4 | 256 | 256 | 1,938,624 | 1.41 | 0.90 |

## Scope

Analytical capacity model with the integrated_v1 component instruments. Concurrency widths, per-lane and per-port tile charges, and the source-to-data transfer are declared assumptions, not compiled schedules. Richardson residuals are the measured finite-patch envelope of |(4E(2r)-E(r))/3| r^4 at r >= 8 tau. The time-growing family uses the running maximum of all measured patches through each time; no patch larger than 4x4 is measured. Interacting bulk convergence is not established.

Raw-RUS gains rest on first-order coherent cancellation in the restricted instrument; Toshio et al. App. E finds raw RUS slightly worse than PCEC. See [RAW_RUS_NOTE.md](../../docs/RAW_RUS_NOTE.md).

Both encoded arms use the reference surface-code law (A = 0.1, Λ = 10 at p = 10⁻³). The Stim calibration in [SURFACE_LAW.md](../../docs/SURFACE_LAW.md) gives Λ = 7.4, about +4–6 in distance; see [round v2](../optimization_round_v2/summary.md) for the propagated result.

Validation: PASS: integrated_v1/physics_sweep_v1 replay, budgets, fault ledgers, switch domains, PCEC >= raw. Reproduce with `make optimization-round` (compute node recommended; `pbs/optimization_round.pbs`), redraw with `make optimization-round-figures`, audit with `make optimization-round-check`.
