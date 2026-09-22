"""Analytical reference curves, not a compiled or calibrated device estimate.

Angles use exp(-i theta P). See analytical_forms.tex for derivations and scope.
Only the standard library is required for the model itself.
"""
from dataclasses import dataclass
from math import ceil, exp, isfinite, log, log1p, log2


@dataclass(frozen=True)
class Workload:
    L: int
    tau: float
    steps: int
    algorithm: str
    # (magnitude in exp(-i theta P), multiplicity), before Clifford extraction.
    dynamic_angles: tuple
    fixed_T: int

    @property
    def sites(self):
        return self.L**2

    @property
    def rotations(self):
        return sum(count for _, count in self.dynamic_angles)

    @property
    def theta(self):
        return sum(angle * count for angle, count in self.dynamic_angles)


def workload(L, tau, steps, algorithm="plaquette", u_over_j=4.0):
    if not isinstance(L, int) or L < 2 or L % 2:
        raise ValueError("This first model requires even L >= 2 and open boundaries")
    if not isfinite(tau) or tau <= 0 or not isinstance(steps, int) or steps < 1:
        raise ValueError("Require positive finite tau and integer steps >= 1")
    if not isfinite(u_over_j) or u_over_j <= 0:
        raise ValueError("This count model requires finite U/J > 0")
    N, E, dt = L * L, 2 * L * (L - 1), tau / steps
    u = u_over_j / 4
    if algorithm == "termwise":
        # K=4E+N ordered Pauli terms, onsite term first and last. Merge
        # neighboring equal terms at each Strang center and step boundary.
        # Last term is an onsite ZZ. No commuting/Clifford optimizations.
        hist = ((dt / 4, 8 * E * steps),
                (u * dt / 2, 2 * steps * (N - 2) + 2),
                (u * dt, 2 * steps - 1))
        fixed = 0
    elif algorithm == "plaquette":
        # HI/2, Hp/2, Hg, Hp/2, HI/2. Open boundary gold tiles
        # include B isolated bonds. Four-mode tile: 8 T + 2 rotations/spin.
        P, G, B = (L // 2)**2, (L // 2 - 1)**2, 2 * (L - 2)
        pink_calls = steps if L == 2 else 2 * steps  # Hg empty for L=2
        hist = ((u * dt / 2, 2 * N), (u * dt, (steps - 1) * N),
                (dt if L == 2 else dt / 2, 4 * P * pink_calls),
                (dt, 4 * G * steps), (dt / 2, 4 * B * steps))
        fixed = 16 * (P * pink_calls + G * steps)
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")
    return Workload(L, tau, steps, algorithm,
                    tuple((a, n) for a, n in hist if n), fixed)


def accuracy(c):
    eps = c["relative_error"] * c["doublon_reference"]
    if not 0 < c["doublon_reference"] <= 1 or not 0 < eps < 1:
        raise ValueError("Specify a nonzero intensive reference for relative error")
    if not 0 < c["delta"] < 1:
        raise ValueError("Require 0 < delta < 1")
    fractions = [c[key + "_fraction"] for key in
                 ("stat", "trotter", "logical", "synthesis", "magic")]
    if min(fractions) <= 0 or abs(sum(fractions) - 1) > 1e-12:
        raise ValueError("Positive error-budget fractions must sum to one")
    # Center D in [-1/2,1/2]; ideal PEC samples then lie in [-gamma/2,gamma/2].
    n0 = ceil(log(2 / c["delta"]) / (2 * (eps * c["stat_fraction"])**2))
    return eps, n0


def log_gamma2_dephasing(q):
    if not 0 <= q < 0.5:
        raise ValueError("Invertible dephasing reference requires 0 <= q < 1/2")
    return -2 * log1p(-2 * q)


def log_gamma2_depolarizing(p, qubits):
    """p = total probability of a uniformly random nonidentity Pauli."""
    k = 4**qubits
    if not 0 <= p < (k - 1) / k:
        raise ValueError("Depolarizing channel is outside the reference branch")
    # Inverse coefficients: identity 1+p/s, others -p/((k-1)s).
    s = 1 - k * p / (k - 1)
    return 2 * log1p(2 * p / s)


def star_exposure(w, c, arm):
    q0 = 4 * c["p"] / 15
    fixed = w.fixed_T * log_gamma2_dephasing(q0)
    if arm == "star_original":
        return (w.rotations + w.fixed_T) * log_gamma2_dephasing(q0)
    floor = c["star_floor"] if arm == "star_floor" else 0.0
    # Fixed T gates use the original injection reference; small-angle fit
    # is applied only to variable-angle rotations, with RUS already folded in.
    return fixed + sum(n * log_gamma2_dephasing(c["star_alpha"] * c["p"] * a + floor)
                       for a, n in w.dynamic_angles)


def physical_proxy(w, c):
    R, T, q = w.rotations, w.fixed_T, 2 * w.sites
    g2 = c["g2_per_dynamic_rotation"] * R + c["g2_per_fixed_T"] * T
    g1 = c["g1_per_dynamic_rotation"] * R + c["g1_per_fixed_T"] * T
    busy = 2 * g2 + g1
    layers = ceil(c["scheduling_overhead"] * busy / q)
    idle = max(0, q * layers - busy)
    chi = (g2 * log_gamma2_depolarizing(c["p"], 2)
           + (g1 + w.sites) * log_gamma2_depolarizing(c["p"] / 10, 1)
           + idle * log_gamma2_depolarizing(c["p"], 1)
           + (q + 2) * log_gamma2_dephasing(c["p"]))
    # q preparation flips, N Neel X gates, two final observable-bit flips.
    seconds = c["reset_seconds"] + (layers + 1) * c["gate_seconds"] + c["measurement_seconds"]
    return chi, seconds, g2, idle


def minimum_distance(tiles, seconds_at_d, log_budget, c):
    if not 0 < c["p"] < c["surface_threshold"]:
        raise ValueError("Reference surface law requires 0 < p < p_th")
    for d in range(3, c["max_distance"] + 1, 2):
        seconds = seconds_at_d(d)
        cycles = ceil(seconds / c["cycle_seconds"])
        log_fault_bound = (log(tiles * cycles * c["surface_A"])
                           + (d + 1) / 2 * log(c["p"] / c["surface_threshold"]))
        if log_fault_bound <= log_budget:
            return d, seconds, log_fault_bound
    return None


def factory(level, p):
    """Deliberately serial 15-to-1 reference; shared d for all stages.

    Level two uses two 11-tile blocks and 15 storage tiles. Produce its
    15 inputs sequentially; no unmodeled pipelined supply is assumed.
    Acceptance and output infidelity use leading-order 15-to-1 formulas.
    """
    a1, q1 = 1 - 15 * p, 35 * p**3
    if a1 <= 0:
        raise ValueError("15-to-1 leading-order reference is outside its regime")
    if level == 1:
        return 12, 11 / a1, q1
    if level == 2:
        return 37, (15 * 11 / a1 + 11) / (1 - 15 * q1), 35 * q1**3
    raise ValueError("Only one or two levels are implemented")


def finish(w, arm, chi, qcopy, seconds, c, **details):
    _, n0 = accuracy(c)
    log_shots = log(n0) + chi
    # Throughput approximation; ceil replicas where it matters. Geometric
    # placement is deferred. Q denotes an area budget in this first model.
    if seconds > c["week_seconds"]:
        log_qmin = float("inf")
    else:
        lm = max(0, log_shots + log(seconds / c["week_seconds"]))
        if lm < 30:
            lm = log(ceil(exp(lm)))
        log_qmin = (log(qcopy) + lm) / log(10)
    replicas = c["device_qubits"] // qcopy
    runtime = (max(log(seconds), log_shots + log(seconds) - log(replicas)) / log(10)
               if replicas else float("inf"))
    return dict(L=w.L, N=w.sites, tau=w.tau, r=w.steps, algorithm=w.algorithm,
                arm=arm, rotations=w.rotations, fixed_T=w.fixed_T, theta=w.theta,
                log10_gamma2=chi / log(10), seconds_per_attempt=seconds,
                qubits_per_copy=qcopy, log10_Q_min=log_qmin,
                log10_seconds_at_Q=runtime, **details)


def estimate(w, c, arm):
    eps, _ = accuracy(c)
    if arm == "nisq":
        chi, seconds, g2, idle = physical_proxy(w, c)
        return finish(w, arm, chi, 2 * w.sites, seconds, c,
                      distance=0, factory_copies=0, factory_level=0,
                      T_count=0, two_qubit_gates=g2, idle_locations=idle)
    data_tiles = ceil(1.5 * (2 * w.sites + 6))
    if arm.startswith("star_"):
        chi = star_exposure(w, c, arm)
        tiles = data_tiles + c["star_injection_tiles"]
        clocks = 2 + c["star_clocks_per_rotation"] * (w.rotations + w.fixed_T)
        def duration(d):
            return (clocks * d * c["cycle_seconds"] + c["reset_seconds"]
                    + c["measurement_seconds"]
                    + (w.rotations + w.fixed_T) * c["feedforward_seconds"])
        solved = minimum_distance(tiles, duration,
                                  log(eps * c["logical_fraction"]) - chi / 2, c)
        if solved is None:
            raise ValueError("Increase max_distance; no admissible STAR distance")
        d, seconds, _ = solved
        return finish(w, arm, chi, 2 * tiles * d * d, seconds, c,
                      distance=d, factory_copies=0, factory_level=0, T_count=0)
    if arm != "ft":
        raise ValueError(f"Unknown arm: {arm}")
    # ||U-U_tilde|| <= sum eta_j; observable bias <= 2 sum eta_j for ||D||<=1.
    eta = eps * c["synthesis_fraction"] / (2 * w.rotations)
    t_per_rotation = ceil(3.02 * log2(1 / eta) + 0.77)
    nt = w.fixed_T + w.rotations * t_per_rotation
    candidates = []
    for level in (1, 2):
        ftiles, fclocks, qmagic = factory(level, c["p"])
        if nt * qmagic > eps * c["magic_fraction"]:
            continue
        for copies in c["factory_copies"]:
            tiles = data_tiles + copies * ftiles
            # Nine clocks per T in the compact data block; add 2 per
            # variable rotation as a tunable Clifford-work placeholder.
            clocks = 2 + max(9 * nt + 2 * w.rotations, ceil(nt / copies) * fclocks)
            def duration(d):
                return (clocks * d * c["cycle_seconds"] + nt * c["feedforward_seconds"]
                        + c["reset_seconds"] + c["measurement_seconds"])
            solved = minimum_distance(tiles, duration, log(eps * c["logical_fraction"]), c)
            if solved is None:
                continue
            d, seconds, _ = solved
            candidates.append(finish(w, arm, 0, 2 * tiles * d * d, seconds, c,
                              distance=d, factory_copies=copies, factory_level=level,
                              T_count=nt))
    if not candidates:
        raise ValueError("No admissible FT factory/distance configuration")
    # Optimize area for the week cutoff, not runtime at the separate fixed Q.
    best = min(candidates, key=lambda x: x["log10_Q_min"])
    fastest = min(candidates, key=lambda x: x["log10_seconds_at_Q"])
    best["log10_seconds_at_Q"] = fastest["log10_seconds_at_Q"]
    for key in ("distance", "factory_copies", "factory_level", "qubits_per_copy", "seconds_per_attempt"):
        best["fixed_Q_" + key] = fastest[key]
    return best
