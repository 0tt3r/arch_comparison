# Optimization round v2: calibrated surface-code law

Identical to [round v1](../optimization_round_v1/summary.md) except that both encoded arms use the
Stim-calibrated surface-code law of [SURFACE_LAW.md](../../docs/SURFACE_LAW.md) (per-p fit; Λ = 7.41 at p = 10⁻³,
89 at p = 10⁻⁴; conservative X+Z amplitude). Physical-qubit rows replay round v1 exactly (1092 checks).
Improvement factors below are measured outputs of the declared model, not assumptions.

[Headline figure (τ = L/4)](../../figures/optimization_round_capacity_v4.pdf) · [Fixed τ = 1 benchmark](../../figures/optimization_round_fixed_tau_v4.pdf) · [Leader map](../../figures/optimization_round_leaders_v4.pdf) · [Methods](../../docs/OPTIMIZATION_ROUND.md)

## Encoded capacity leader by configuration

Entries give the leader and its lattice side at Q = 10⁶ and 10⁸. The time-growing family (τ = L/4) is the headline;
fixed τ = 1 is retained as an implementation benchmark.

### τ = 1

| p | ε | Q | v1 menu | + raw RUS | + lanes | + FT delivery | + separate distances | + Richardson | all five |
|---:|---:|---:|---|---|---|---|---|---|---|
| 0.001 | 0.01 | 10^6 | ft (8/10/12) | tie:hybrid+ft (10/12/12) | ft (8/10/12) | ft (8/10/12) | ft (8/10/12) | ft (8/10/12) | tie:hybrid+ft (10/12/12) |
| 0.001 | 0.01 | 10^8 | ft (10/16/48) | ft (12/18/48) | ft (12/16/48) | ft (10/16/48) | ft (10/16/48) | ft (10/16/48) | ft (12/20/48) |
| 0.001 | 0.001 | 10^6 | tie:star+hybrid+ft (4/4/4) | tie:star+hybrid+ft (4/4/4) | hybrid (4/6/4) | ft (4/4/8) | tie:star+hybrid+ft (4/4/4) | hybrid (4/6/4) | tie:hybrid+ft (4/8/8) |
| 0.001 | 0.001 | 10^8 | ft (6/8/16) | ft (6/8/16) | ft (6/8/16) | ft (6/8/38) | ft (6/8/16) | ft (6/10/18) | ft (8/12/40) |
| 0.0001 | 0.01 | 10^6 | tie:hybrid+ft (16/22/22) | tie:hybrid+ft (16/22/22) | tie:hybrid+ft (18/22/22) | ft (16/22/26) | tie:hybrid+ft (16/22/22) | tie:hybrid+ft (16/22/22) | ft (18/24/26) |
| 0.0001 | 0.01 | 10^8 | ft (22/36/48) | ft (28/38/48) | ft (24/38/48) | ft (22/36/48) | ft (22/36/48) | ft (22/36/48) | ft (32/46/48) |
| 0.0001 | 0.001 | 10^6 | hybrid (6/8/6) | hybrid (6/8/6) | hybrid (6/10/6) | tie:hybrid+ft (6/10/10) | hybrid (6/8/6) | hybrid (8/12/8) | tie:hybrid+ft (8/16/16) |
| 0.0001 | 0.001 | 10^8 | ft (12/16/28) | ft (14/18/28) | ft (14/20/28) | ft (12/16/48) | ft (12/16/28) | ft (12/20/32) | ft (18/26/48) |

Cells: leader (STAR/hybrid/FT sides).

### τ = L/4

| p | ε | Q | v1 menu | + raw RUS | + lanes | + FT delivery | + separate distances | + Richardson | all five |
|---:|---:|---:|---|---|---|---|---|---|---|
| 0.001 | 0.01 | 10^6 | ft (6/8/10) | ft (6/8/10) | ft (6/8/10) | ft (6/8/12) | ft (6/8/10) | ft (6/8/10) | ft (6/10/12) |
| 0.001 | 0.01 | 10^8 | ft (8/10/32) | ft (8/10/32) | ft (8/10/32) | ft (8/10/32) | ft (8/10/32) | ft (8/10/32) | ft (8/12/32) |
| 0.001 | 0.001 | 10^6 | tie:star+hybrid+ft (4/4/4) | tie:star+hybrid+ft (4/4/4) | tie:star+hybrid+ft (4/4/4) | tie:star+hybrid+ft (4/4/4) | tie:star+hybrid+ft (4/4/4) | tie:star+hybrid+ft (4/4/4) | tie:hybrid+ft (4/6/6) |
| 0.001 | 0.001 | 10^8 | ft (6/6/10) | ft (6/6/10) | ft (6/6/10) | ft (6/6/18) | ft (6/6/10) | ft (6/6/12) | ft (6/8/18) |
| 0.0001 | 0.01 | 10^6 | ft (10/14/18) | ft (10/14/18) | ft (10/14/18) | ft (10/14/20) | ft (10/14/18) | ft (10/14/18) | ft (12/16/20) |
| 0.0001 | 0.01 | 10^8 | ft (12/18/32) | ft (16/18/32) | ft (14/18/32) | ft (12/18/32) | ft (12/18/32) | ft (12/18/32) | ft (16/20/32) |
| 0.0001 | 0.001 | 10^6 | tie:star+hybrid (6/6/4) | tie:star+hybrid (6/6/4) | tie:star+hybrid (6/6/4) | tie:star+hybrid+ft (6/6/6) | tie:star+hybrid (6/6/4) | hybrid (6/8/6) | hybrid (6/10/8) |
| 0.0001 | 0.001 | 10^8 | ft (8/10/18) | ft (10/12/18) | ft (10/12/18) | ft (8/10/32) | ft (8/10/18) | ft (8/12/18) | ft (10/14/32) |

Cells: leader (STAR/hybrid/FT sides).

## How often each change moves the leader

Out of 64 (family, p, ε, Q) cells, relative to the v1 menu:

| Configuration | Cells with a different encoded leader |
|---|---:|
| + raw RUS | 2 |
| + lanes | 4 |
| + FT delivery | 10 |
| + separate distances | 2 |
| + Richardson | 4 |
| all five | 15 |

## What each change buys at representative points

log₁₀ Q_min for one week; — means no feasible design within the declared menu.

| Family | p | ε | L | Arm | v1 menu | + raw RUS | + lanes | + FT delivery | + separate distances | + Richardson | all five |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| τ=1 | 0.001 | 0.01 | 8 | star | 5.52 | 5.52 | 5.52 | 5.52 | 5.52 | 5.52 | 5.50 |
| τ=1 | 0.001 | 0.01 | 8 | hybrid | 5.52 | 5.52 | 5.49 | 5.52 | 5.52 | 5.52 | 5.49 |
| τ=1 | 0.001 | 0.01 | 8 | ft | 5.59 | 5.59 | 5.59 | 5.54 | 5.59 | 5.59 | 5.54 |
| τ=1 | 0.001 | 0.01 | 8 | nisq_effective | 8.29 | 8.29 | 8.29 | 8.29 | 8.29 | 8.29 | 8.29 |
| τ=1 | 0.001 | 0.001 | 8 | star | 8.96 | 8.64 | 8.84 | 8.96 | 8.96 | 8.09 | 7.59 |
| τ=1 | 0.001 | 0.001 | 8 | hybrid | 7.47 | 7.14 | 7.18 | 7.36 | 7.47 | 6.65 | 5.85 |
| τ=1 | 0.001 | 0.001 | 8 | ft | 7.08 | 7.08 | 7.08 | 5.93 | 7.08 | 6.76 | 5.86 |
| τ=1 | 0.001 | 0.001 | 8 | nisq_effective | 29.40 | 29.40 | 29.40 | 29.40 | 29.40 | 22.88 | 22.88 |
| τ=1 | 0.0001 | 0.01 | 8 | star | 4.89 | 4.89 | 4.85 | 4.89 | 4.89 | 4.89 | 4.85 |
| τ=1 | 0.0001 | 0.01 | 8 | hybrid | 5.02 | 5.02 | 4.98 | 5.02 | 4.89 | 5.02 | 4.86 |
| τ=1 | 0.0001 | 0.01 | 8 | ft | 5.02 | 5.02 | 5.02 | 4.97 | 4.91 | 5.02 | 4.87 |
| τ=1 | 0.0001 | 0.01 | 8 | nisq_effective | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 |
| τ=1 | 0.0001 | 0.001 | 8 | star | 6.22 | 6.19 | 6.10 | 6.22 | 6.22 | 5.61 | 5.50 |
| τ=1 | 0.0001 | 0.001 | 8 | hybrid | 5.48 | 5.48 | 5.21 | 5.48 | 5.48 | 5.18 | 5.07 |
| τ=1 | 0.0001 | 0.001 | 8 | ft | 6.35 | 6.35 | 6.35 | 5.39 | 6.35 | 5.84 | 5.23 |
| τ=1 | 0.0001 | 0.001 | 8 | nisq_effective | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 |
| τ=1 | 0.001 | 0.01 | 16 | star | 11.67 | 9.66 | 11.23 | 11.67 | 11.67 | 11.67 | 9.21 |
| τ=1 | 0.001 | 0.01 | 16 | hybrid | 7.76 | 6.64 | 7.50 | 7.65 | 7.76 | 7.76 | 6.29 |
| τ=1 | 0.001 | 0.01 | 16 | ft | 6.20 | 6.20 | 6.20 | 6.19 | 6.20 | 6.20 | 6.19 |
| τ=1 | 0.001 | 0.01 | 16 | nisq_effective | 97.64 | 97.64 | 97.64 | 97.64 | 97.64 | 97.64 | 97.64 |
| τ=1 | 0.001 | 0.001 | 16 | star | 18.54 | 18.05 | 18.06 | 18.54 | 18.54 | 16.96 | 15.08 |
| τ=1 | 0.001 | 0.001 | 16 | hybrid | 12.15 | 11.39 | 11.85 | 11.99 | 12.15 | 10.84 | 9.23 |
| τ=1 | 0.001 | 0.001 | 16 | ft | 7.97 | 7.97 | 7.97 | 7.26 | 7.97 | 7.77 | 6.52 |
| τ=1 | 0.001 | 0.001 | 16 | nisq_effective | 231.85 | 231.85 | 231.85 | 231.85 | 231.85 | 200.73 | 200.73 |
| τ=1 | 0.0001 | 0.01 | 16 | star | 5.56 | 5.56 | 5.56 | 5.56 | 5.56 | 5.56 | 5.56 |
| τ=1 | 0.0001 | 0.01 | 16 | hybrid | 5.56 | 5.56 | 5.55 | 5.56 | 5.56 | 5.56 | 5.55 |
| τ=1 | 0.0001 | 0.01 | 16 | ft | 5.56 | 5.56 | 5.56 | 5.55 | 5.56 | 5.56 | 5.55 |
| τ=1 | 0.0001 | 0.01 | 16 | nisq_effective | 6.62 | 6.62 | 6.62 | 6.62 | 6.62 | 6.62 | 6.62 |
| τ=1 | 0.0001 | 0.001 | 16 | star | 8.51 | 8.03 | 8.08 | 8.51 | 8.51 | 8.51 | 7.61 |
| τ=1 | 0.0001 | 0.001 | 16 | hybrid | 7.90 | 7.68 | 7.32 | 7.76 | 7.90 | 7.01 | 5.90 |
| τ=1 | 0.0001 | 0.001 | 16 | ft | 7.09 | 7.09 | 7.09 | 6.50 | 7.09 | 6.94 | 5.91 |
| τ=1 | 0.0001 | 0.001 | 16 | nisq_effective | 24.26 | 24.26 | 24.26 | 24.26 | 24.26 | 18.69 | 18.69 |
| τ=L/4 | 0.001 | 0.01 | 8 | star | 7.04 | 6.12 | 6.94 | 7.04 | 7.04 | 7.04 | 6.07 |
| τ=L/4 | 0.001 | 0.01 | 8 | hybrid | 5.59 | 5.52 | 5.59 | 5.59 | 5.59 | 5.59 | 5.50 |
| τ=L/4 | 0.001 | 0.01 | 8 | ft | 5.61 | 5.61 | 5.61 | 5.55 | 5.61 | 5.61 | 5.55 |
| τ=L/4 | 0.001 | 0.01 | 8 | nisq_effective | 21.30 | 21.30 | 21.30 | 21.30 | 21.30 | 21.30 | 21.30 |
| τ=L/4 | 0.001 | 0.001 | 8 | star | 11.99 | 11.41 | 11.86 | 11.99 | 11.99 | 10.76 | 9.90 |
| τ=L/4 | 0.001 | 0.001 | 8 | hybrid | 9.18 | 8.73 | 8.98 | 9.04 | 9.18 | 8.26 | 6.88 |
| τ=L/4 | 0.001 | 0.001 | 8 | ft | 7.43 | 7.43 | 7.43 | 6.98 | 7.43 | 7.25 | 6.00 |
| τ=L/4 | 0.001 | 0.001 | 8 | nisq_effective | 60.66 | 60.66 | 60.66 | 60.66 | 60.66 | 48.51 | 48.51 |
| τ=L/4 | 0.0001 | 0.01 | 8 | star | 4.89 | 4.89 | 4.85 | 4.89 | 4.89 | 4.89 | 4.85 |
| τ=L/4 | 0.0001 | 0.01 | 8 | hybrid | 5.02 | 5.02 | 4.98 | 5.02 | 4.89 | 5.02 | 4.86 |
| τ=L/4 | 0.0001 | 0.01 | 8 | ft | 5.02 | 5.02 | 5.02 | 4.97 | 5.02 | 5.02 | 4.94 |
| τ=L/4 | 0.0001 | 0.01 | 8 | nisq_effective | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 | 2.16 |
| τ=L/4 | 0.0001 | 0.001 | 8 | star | 7.12 | 6.89 | 7.00 | 7.12 | 7.12 | 6.67 | 6.50 |
| τ=L/4 | 0.0001 | 0.001 | 8 | hybrid | 6.44 | 6.36 | 6.18 | 6.28 | 6.44 | 5.65 | 5.20 |
| τ=L/4 | 0.0001 | 0.001 | 8 | ft | 6.72 | 6.72 | 6.72 | 6.38 | 6.72 | 6.58 | 5.40 |
| τ=L/4 | 0.0001 | 0.001 | 8 | nisq_effective | 5.74 | 5.74 | 5.74 | 5.74 | 5.74 | 3.06 | 3.06 |
| τ=L/4 | 0.001 | 0.01 | 16 | star | 33.07 | 25.09 | 32.60 | 33.07 | 33.07 | 33.07 | 24.62 |
| τ=L/4 | 0.001 | 0.01 | 16 | hybrid | 17.76 | 13.28 | 17.42 | 17.64 | 17.76 | 17.76 | 12.61 |
| τ=L/4 | 0.001 | 0.01 | 16 | ft | 6.35 | 6.35 | 6.35 | 6.28 | 6.35 | 6.35 | 6.28 |
| τ=L/4 | 0.001 | 0.01 | 16 | nisq_effective | 401.90 | 401.90 | 401.90 | 401.90 | 401.90 | 401.90 | 401.90 |
| τ=L/4 | 0.001 | 0.001 | 16 | star | 49.03 | 49.03 | 48.57 | 49.03 | 49.03 | 47.19 | 46.75 |
| τ=L/4 | 0.001 | 0.001 | 16 | hybrid | 26.63 | 26.63 | 26.34 | 26.50 | 26.63 | 23.28 | 22.57 |
| τ=L/4 | 0.001 | 0.001 | 16 | ft | 8.67 | 8.67 | 8.67 | 7.80 | 8.67 | 8.52 | 7.58 |
| τ=L/4 | 0.001 | 0.001 | 16 | nisq_effective | 1026.38 | 1026.38 | 1026.38 | 1026.38 | 1026.38 | 808.24 | 808.24 |
| τ=L/4 | 0.0001 | 0.01 | 16 | star | 9.75 | 7.84 | 9.32 | 9.75 | 9.75 | 9.75 | 7.43 |
| τ=L/4 | 0.0001 | 0.01 | 16 | hybrid | 6.92 | 6.66 | 6.74 | 6.79 | 6.92 | 6.92 | 5.84 |
| τ=L/4 | 0.0001 | 0.01 | 16 | ft | 5.72 | 5.72 | 5.72 | 5.64 | 5.72 | 5.72 | 5.64 |
| τ=L/4 | 0.0001 | 0.01 | 16 | nisq_effective | 37.81 | 37.81 | 37.81 | 37.81 | 37.81 | 37.81 | 37.81 |
| τ=L/4 | 0.0001 | 0.001 | 16 | star | 12.18 | 10.22 | 11.75 | 12.18 | 12.18 | 12.18 | 9.80 |
| τ=L/4 | 0.0001 | 0.001 | 16 | hybrid | 10.49 | 9.30 | 9.76 | 10.48 | 10.49 | 10.39 | 8.31 |
| τ=L/4 | 0.0001 | 0.001 | 16 | ft | 7.84 | 7.84 | 7.84 | 6.99 | 7.84 | 7.66 | 6.97 |
| τ=L/4 | 0.0001 | 0.001 | 16 | nisq_effective | 104.78 | 104.78 | 104.78 | 104.78 | 104.78 | 80.06 | 80.06 |

## Selected designs with all five changes

What the optimizer chose, which is the mechanism behind each shift.

| Family | p | ε | L | Arm | log₁₀ Q | r | extrapolation | d | RUS | lanes | ports | factories | HWP batch | qubits/copy | s/attempt | χ |
|---|---:|---:|---:|---|---:|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| τ=1 | 0.001 | 0.01 | 8 | star | 5.50 | 9 | single | 27 | raw | 2 | 1 | 0 | — | 314,928 | 0.754 | 2.73 |
| τ=1 | 0.001 | 0.01 | 8 | hybrid | 5.49 | 10 | single | 27 | raw | 1 | 1 | 1 | — | 308,662 | 0.862 | 1.13 |
| τ=1 | 0.001 | 0.01 | 8 | ft | 5.54 | 13 | single | 29 | — | — | 1 | 1 | 1 | 345,834 | 17.8 | 0.00 |
| τ=1 | 0.001 | 0.001 | 8 | star | 7.59 | 8 | richardson | 33 | raw | 8 | 1 | 0 | — | 762,300 | 0.4 | 4.81 |
| τ=1 | 0.001 | 0.001 | 8 | hybrid | 5.85 | 8 | richardson | 29 | raw | 16 | 4 | 8 | — | 701,176 | 0.0688 | 2.34 |
| τ=1 | 0.001 | 0.001 | 8 | ft | 5.86 | 8 | richardson | 29 | — | — | 8 | 32 | 32 | 729,116 | 0.135 | 0.90 |
| τ=1 | 0.0001 | 0.01 | 8 | star | 4.85 | 8 | single | 13 | raw | 1 | 1 | 0 | — | 70,980 | 0.542 | 0.88 |
| τ=1 | 0.0001 | 0.01 | 8 | hybrid | 4.86 | 12 | single | 13 | raw | 1 | 1 | 1 | — | 72,342 | 0.447 | 0.38 |
| τ=1 | 0.0001 | 0.01 | 8 | ft | 4.87 | 8 | single | 13 | — | — | 1 | 2 | 4 | 73,366 | 3.5 | 0.00 |
| τ=1 | 0.0001 | 0.001 | 8 | star | 5.50 | 8 | richardson | 15 | raw | 8 | 1 | 0 | — | 157,500 | 0.17 | 2.49 |
| τ=1 | 0.0001 | 0.001 | 8 | hybrid | 5.07 | 8 | richardson | 13 | raw | 4 | 2 | 4 | — | 116,312 | 0.0918 | 1.36 |
| τ=1 | 0.0001 | 0.001 | 8 | ft | 5.23 | 8 | richardson | 13 | — | — | 4 | 32 | 32 | 169,996 | 0.109 | 0.90 |
| τ=L/4 | 0.001 | 0.01 | 8 | star | 6.07 | 16 | single | 29 | raw | 8 | 1 | 0 | — | 588,700 | 0.362 | 5.11 |
| τ=L/4 | 0.001 | 0.01 | 8 | hybrid | 5.50 | 18 | single | 27 | raw | 2 | 1 | 1 | — | 317,410 | 1.12 | 2.18 |
| τ=L/4 | 0.001 | 0.01 | 8 | ft | 5.55 | 19 | single | 29 | — | — | 1 | 1 | 4 | 354,244 | 20.9 | 0.00 |
| τ=L/4 | 0.001 | 0.001 | 8 | star | 9.90 | 16 | richardson | 37 | raw | 8 | 1 | 0 | — | 958,300 | 0.896 | 8.79 |
| τ=L/4 | 0.001 | 0.001 | 8 | hybrid | 6.88 | 16 | richardson | 31 | raw | 64 | 32 | 32 | — | 1,531,928 | 0.0316 | 3.78 |
| τ=L/4 | 0.001 | 0.001 | 8 | ft | 6.00 | 16 | richardson | 29 | — | — | 32 | 64 | 64 | 1,009,138 | 0.115 | 0.90 |
| τ=L/4 | 0.0001 | 0.01 | 8 | star | 4.85 | 17 | single | 13 | raw | 1 | 1 | 0 | — | 70,980 | 1.15 | 1.87 |
| τ=L/4 | 0.0001 | 0.01 | 8 | hybrid | 4.86 | 23 | single | 13 | raw | 1 | 1 | 1 | — | 72,342 | 0.854 | 0.72 |
| τ=L/4 | 0.0001 | 0.01 | 8 | ft | 4.94 | 16 | single | 13 | — | — | 1 | 1 | 32 | 86,200 | 2.98 | 0.00 |
| τ=L/4 | 0.0001 | 0.001 | 8 | star | 6.50 | 16 | richardson | 15 | raw | 8 | 1 | 0 | — | 157,500 | 0.34 | 4.12 |
| τ=L/4 | 0.0001 | 0.001 | 8 | hybrid | 5.20 | 16 | richardson | 13 | raw | 16 | 4 | 16 | — | 158,344 | 0.0543 | 1.82 |
| τ=L/4 | 0.0001 | 0.001 | 8 | ft | 5.40 | 16 | richardson | 15 | — | — | 8 | 32 | 64 | 249,218 | 0.143 | 0.90 |
| τ=L/4 | 0.001 | 0.01 | 16 | star | 24.62 | 32 | single | 55 | raw | 16 | 1 | 0 | — | 7,247,900 | 2.92 | 43.34 |
| τ=L/4 | 0.001 | 0.01 | 16 | hybrid | 12.61 | 32 | single | 39 | raw | 128 | 32 | 512 | — | 26,648,984 | 0.132 | 16.47 |
| τ=L/4 | 0.001 | 0.01 | 16 | ft | 6.28 | 41 | single | 27 | — | — | 16 | 64 | 32 | 1,902,616 | 0.566 | 0.00 |
| τ=L/4 | 0.001 | 0.001 | 16 | star | 46.75 | 32 | richardson | 85 | noisy_pcec | 16 | 1 | 0 | — | 17,311,100 | 8.76 | 88.82 |
| τ=L/4 | 0.001 | 0.001 | 16 | hybrid | 22.57 | 32 | richardson | 53 | noisy_pcec | 128 | 32 | 512 | — | 48,771,672 | 0.313 | 35.49 |
| τ=L/4 | 0.001 | 0.001 | 16 | ft | 7.58 | 32 | richardson | 35 | — | — | 32 | 1024 | 256 | 37,798,576 | 0.491 | 0.90 |
| τ=L/4 | 0.0001 | 0.01 | 16 | star | 7.43 | 41 | single | 25 | raw | 16 | 1 | 0 | — | 1,497,500 | 16.8 | 4.43 |
| τ=L/4 | 0.0001 | 0.01 | 16 | hybrid | 5.84 | 34 | single | 15 | raw | 32 | 64 | 32 | — | 686,168 | 0.103 | 4.32 |
| τ=L/4 | 0.0001 | 0.01 | 16 | ft | 5.64 | 40 | single | 15 | — | — | 1 | 1 | 128 | 438,424 | 17.1 | 0.00 |
| τ=L/4 | 0.0001 | 0.001 | 16 | star | 9.80 | 186 | single | 25 | raw | 16 | 1 | 0 | — | 1,497,500 | 34.8 | 4.69 |
| τ=L/4 | 0.0001 | 0.001 | 16 | hybrid | 8.31 | 41 | richardson | 25 | raw | 128 | 64 | 128 | — | 4,706,072 | 1.2 | 3.67 |
| τ=L/4 | 0.0001 | 0.001 | 16 | ft | 6.97 | 32 | richardson | 17 | — | — | 4 | 256 | 256 | 3,129,024 | 1.48 | 0.90 |

## Scope

Analytical capacity model with the integrated_v1 component instruments. Concurrency widths, per-lane and per-port tile charges, and the source-to-data transfer are declared assumptions, not compiled schedules. Richardson residuals are the measured finite-patch envelope of |(4E(2r)-E(r))/3| r^4 at r >= 8 tau. The time-growing family uses the running maximum of all measured patches through each time; no patch larger than 4x4 is measured. Interacting bulk convergence is not established. v2 uses the Stim-calibrated surface-code law for both encoded arms; its fit is memory-only at d=3-9 and is extrapolated to the distances used here.

Raw-RUS gains rest on first-order coherent cancellation in the restricted instrument; Toshio et al. App. E finds raw RUS slightly worse than PCEC. See [RAW_RUS_NOTE.md](../../docs/RAW_RUS_NOTE.md).

Both encoded arms use the Stim-calibrated surface-code law of [SURFACE_LAW.md](../../docs/SURFACE_LAW.md), fitted at d = 3–9 and extrapolated to the distances used here.

Validation: PASS: round-v1 NISQ replay, calibrated surface law, budgets, fault ledgers, switch domains, PCEC >= raw. Reproduce with `make optimization-round` (compute node recommended; `pbs/optimization_round.pbs`), redraw with `make optimization-round-figures`, audit with `make optimization-round-check`.
