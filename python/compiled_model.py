"""Replace analytical circuit proxies with CUDA-Q counts; retain named QEC models."""
from math import ceil, expm1, log, log1p

from model import (accuracy, estimate, factory, finish, log_gamma2_dephasing,
                   log_gamma2_depolarizing, minimum_distance, star_exposure, workload)


def compiled_physical(record, c, idle_convention="per_100ns", routing="jw"):
    if routing not in ('jw','fswap'):
        raise ValueError(routing)
    w = workload(record["L"],record["tau"],record["steps"],record["algorithm"],c["u_over_j"])
    schedule = record["native_schedule" if routing=="jw" else "fswap_schedule"]
    p_idle = c["p"]
    if idle_convention == "per_100ns":
        # Compose the one-qubit depolarizing eigenvalue, so a 100ns idle
        # has total nonidentity-Pauli probability exactly p.
        p_idle = -.75*expm1((schedule["gate_ns"]/100)*log1p(-4*c["p"]/3))
    elif idle_convention != "per_layer":
        raise ValueError(idle_convention)
    pieces = {"two_qubit": schedule["two_qubit_gates"]*log_gamma2_depolarizing(c["p"],2),
              "one_qubit": schedule["one_qubit_gates"]*log_gamma2_depolarizing(c["p"]/10,1),
              "idle": schedule["idle_qubit_layers"]*log_gamma2_depolarizing(p_idle,1),
              "spam": (2*w.sites+2)*log_gamma2_dephasing(c["p"])}
    row = finish(w,"nisq",sum(pieces.values()),schedule["physical_copy_qubits"],
                 schedule["seconds_per_attempt"],c,distance=0,
                 idle_convention=idle_convention,**{"chi_"+key:value for key,value in pieces.items()})
    row["implementation"] = "compiled_"+routing
    return row


def encoded(record, c, arm="ft", logical_bias="worst_case", t_clocks=9,
            injection_error=None):
    """Actual synthesis counts plus the declared tile/factory reference.

    hybrid uses analog dynamic rotations and distilled fixed pi/8 rotations.
    commuting_pauli is a sensitivity requiring branch-independent Pauli faults
    outside complete PEC blocks, or commuting with the inverted channel.
    """
    if logical_bias not in ("worst_case","commuting_pauli"):
        raise ValueError(logical_bias)
    if arm not in ("ft","hybrid"):
        raise ValueError(arm)
    # Keep the v1 call convention reproducible. New scenarios must specify
    # injected-state error separately from physical circuit error.
    pin=c['p'] if injection_error is None else injection_error
    if not 0<pin<1/15:
        raise ValueError('Require 0 < injected-state error < 1/15')
    w = workload(record["L"],record["tau"],record["steps"],record["algorithm"],c["u_over_j"])
    eps,_ = accuracy(c)
    if arm == "ft":
        nt, chi, analog = record["synthesis"]["T_total"], 0, 0
    else:
        nt, analog = w.fixed_T, w.rotations
        # Remove the original-injection channel for the now-distilled T gates.
        chi = star_exposure(w,c,"star_floor") - nt*log_gamma2_dephasing(4*c["p"]/15)
    log_amp = chi/2 if logical_bias == "worst_case" else 0
    candidates = []
    for level in (0,1,2):
        if level == 0:
            if nt:
                continue
            ftiles,fclocks,qmagic = 0,0,0
        else:
            if not nt:
                continue
            ftiles,fclocks,qmagic = factory(level,pin)
            if log(nt*qmagic)+log_amp > log(eps*c["magic_fraction"]):
                continue
        for copies in ([0] if level == 0 else c["factory_copies"]):
            tiles = ceil(1.5*(2*w.sites+6)) + (c["star_injection_tiles"] if analog else 0) + copies*ftiles
            consume = t_clocks*nt + c["star_clocks_per_rotation"]*analog + (2*w.rotations if not analog else 0)
            clocks = 2 + max(consume, ceil(nt/copies)*fclocks if copies else 0)
            def duration(d):
                return (clocks*d*c["cycle_seconds"] + (nt+analog)*c["feedforward_seconds"]
                        + c["reset_seconds"]+c["measurement_seconds"])
            solved = minimum_distance(tiles,duration,log(eps*c["logical_fraction"])-log_amp,c)
            if solved is None:
                continue
            d,seconds,_ = solved
            candidates.append(finish(w,arm,chi,2*tiles*d*d,seconds,c,distance=d,
                                      T_count=nt,factory_copies=copies,factory_level=level,
                                      logical_bias=logical_bias,t_clocks=t_clocks,
                                      injection_error=pin,magic_output_error=qmagic))
    if not candidates:
        raise ValueError("No encoded reference configuration meets the budget")
    result = min(candidates,key=lambda row:row["log10_Q_min"])
    fastest = min(candidates,key=lambda row:row["log10_seconds_at_Q"])
    result['log10_seconds_at_Q']=fastest['log10_seconds_at_Q']
    for key in ('distance','factory_copies','factory_level','qubits_per_copy','seconds_per_attempt'):
        result['fixed_Q_'+key]=fastest[key]
    result["implementation"] = "CUDA_Q_counts_with_reference_QEC_and_supply"
    return result


def compiled_star(record,c,logical_bias="worst_case"):
    if logical_bias not in ("worst_case","commuting_pauli"):
        raise ValueError(logical_bias)
    w = workload(record["L"],record["tau"],record["steps"],record["algorithm"],c["u_over_j"])
    if logical_bias == "worst_case":
        row = estimate(w,c,"star_floor")
    else:
        # Use the same channel/supply as v1; change only the logical-bias bound.
        eps,_ = accuracy(c)
        chi = star_exposure(w,c,"star_floor")
        tiles=ceil(1.5*(2*w.sites+6))+c["star_injection_tiles"]
        clocks=2+c["star_clocks_per_rotation"]*(w.rotations+w.fixed_T)
        def duration(d):
            return (clocks*d*c["cycle_seconds"]+(w.rotations+w.fixed_T)*c["feedforward_seconds"]
                    +c["reset_seconds"]+c["measurement_seconds"])
        d,seconds,_=minimum_distance(tiles,duration,log(eps*c["logical_fraction"]),c)
        row=finish(w,"star_floor",chi,2*tiles*d*d,seconds,c,distance=d)
    row["logical_bias"]=logical_bias
    row["implementation"]="compiled_angle_counts_with_reference_STAR_supply"
    return row
