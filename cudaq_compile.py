"""Run CUDA-Q Logical on the common Hubbard workload and its NN realization.

Compile the analog program, native parity circuit, and Clifford+T path.
Synthesis is cached by signed angle class; an optional whole-program pass
checks the aggregated T count. No dynamics simulation is involved.
"""
import argparse
import hashlib
import importlib.metadata
import json
from math import isclose
from pathlib import Path
import time

import cudaq.logical as cql
from cudaq.logical.algebra import X, Y, Z

from hubbard_circuits import (angle_summary, build_workload, native_path,
                              native_fermionic_route, schedule_native)
from model import accuracy, workload

ROOT = Path(__file__).resolve().parent


def program_from_operations(operations, width, name):
    @cql.program
    def program() -> None:
        qubits = cql.allocate(width, state=cql.types.zero)
        for op in operations:
            if op.name == "rotate":
                product = None
                for index, letter in zip(op.qubits, op.paulis):
                    factor = {"X": X, "Y": Y, "Z": Z}[letter](qubits[index])
                    product = factor if product is None else product @ factor
                old = tuple(qubits[index] for index in op.qubits)
                values = cql.ops.rotate(product, angle=2*op.angle)
                successors = {id(before): after for before,after in zip(product.operands, values)}
                for index, before in zip(op.qubits, old):
                    qubits[index] = successors[id(before)]
            elif op.name == "rz":
                q = op.qubits[0]
                qubits[q] = cql.rz(qubits[q], op.angle)
            elif op.name == "cx":
                a,b = op.qubits
                qubits[a], qubits[b] = cql.cx(qubits[a], qubits[b])
            else:
                q = op.qubits[0]
                qubits[q] = getattr(cql, op.name)(qubits[q])
        # A full bit string also permits symmetry-selection studies later.
        for i in range(width):
            cql.measure_z(qubits[i])
    return program


def logical_counts(build):
    return cql.estimate(build, tier=cql.estimate.Tier.LOGICAL).to_dict()


def t_count(counts):
    return sum(n for name,n in counts["actions"].items()
               if name in ("qlx_standard_t", "qlx_standard_tdg"))


def synthesize_angles(ops, precision, directory):
    from collections import Counter
    # Use exact float hex keys, so no rounding silently coalesces angles.
    histogram = Counter(op.angle.hex() for op in ops if op.name == "rotate" and op.role == "dynamic")
    classes = []
    for i,(hex_angle,count) in enumerate(sorted(histogram.items())):
        theta = float.fromhex(hex_angle)
        from hubbard_circuits import Operation
        program = program_from_operations([Operation("rz", (0,), 2*theta)], 1, f"angle_{i}")
        built = cql.compiler.synthesize(program, gate_set=cql.compiler.clifford_t, precision=precision)
        counts = logical_counts(built)
        (directory / f"angle_{i}_ft.mlir").write_text(built.to_mlir())
        classes.append({"theta": theta, "theta_hex": hex_angle, "multiplicity": count,
                        "T_per_rotation": t_count(counts), "logical_estimate": counts})
    fixed = sum(op.name == "rotate" and op.role == "fixed_T" for op in ops)
    return {"operator_precision_per_rotation": precision, "angle_classes": classes,
            "T_fixed": fixed,
            "T_total": fixed + sum(row["multiplicity"]*row["T_per_rotation"] for row in classes),
            "method": "CUDA-Q GridSynth per signed angle class; additive counts before inter-gate optimization"}


def compile_case(L, tau, steps, algorithm, config, out, full_synthesis=False):
    directory = out / f"{algorithm}_L{L}_tau{tau:g}_r{steps}"
    directory.mkdir(parents=True, exist_ok=True)
    ops, central = build_workload(L,tau,steps,algorithm,config["u_over_j"])
    summary = angle_summary(ops)
    expected = workload(L,tau,steps,algorithm,config["u_over_j"])
    assert summary["dynamic_rotations"] == expected.rotations
    assert summary["fixed_T"] == expected.fixed_T
    assert isclose(summary["theta_dynamic"], expected.theta, rel_tol=1e-11)
    print(f"Compiling {directory.name}: analog {len(ops)} actions", flush=True)
    start = time.monotonic()
    analog = cql.compile(program_from_operations(ops,2*L*L,"hubbard_analog"))
    analog_counts = logical_counts(analog)
    (directory / "analog_p0.mlir").write_text(analog.to_mlir())
    native_ops = native_path(ops)
    print(f"Compiling {directory.name}: native {len(native_ops)} gates", flush=True)
    native = cql.compile(program_from_operations(native_ops,2*L*L,"hubbard_nn"))
    native_counts = logical_counts(native)
    scheduled = schedule_native(native_ops, 2*L*L, gate_ns=config["gate_seconds"]*1e9,
                                reset_ns=config["reset_seconds"]*1e9,
                                measurement_ns=config["measurement_seconds"]*1e9)
    assert native_counts["actions"].get("qlx_standard_cx",0) == scheduled["two_qubit_gates"]
    (directory / "native_p0.mlir").write_text(native.to_mlir())
    swapped,permutation,swaps=native_fermionic_route(ops,2*L*L)
    print(f"Compiling {directory.name}: fSWAP comparator {len(swapped)} gates",flush=True)
    swapped_build=cql.compile(program_from_operations(swapped,2*L*L,"hubbard_fswap"))
    swapped_counts=logical_counts(swapped_build)
    swapped_schedule=schedule_native(swapped,2*L*L,gate_ns=config["gate_seconds"]*1e9,
                                     reset_ns=config["reset_seconds"]*1e9,
                                     measurement_ns=config["measurement_seconds"]*1e9)
    assert swapped_counts["actions"].get("qlx_standard_cx",0)==swapped_schedule["two_qubit_gates"]
    (directory/"fswap_p0.mlir").write_text(swapped_build.to_mlir())
    eps,_ = accuracy(config)
    precision = eps*config["synthesis_fraction"]/(2*summary["dynamic_rotations"])
    synthesized = synthesize_angles(ops,precision,directory)
    if full_synthesis:
        print(f"Compiling {directory.name}: whole-program Clifford+T", flush=True)
        full = cql.compiler.synthesize(analog, gate_set=cql.compiler.clifford_t, precision=precision)
        counts = logical_counts(full)
        synthesized["whole_program"] = counts
        synthesized["whole_program_T"] = t_count(counts)
        synthesized["whole_program_matches_angle_classes"] = t_count(counts) == synthesized["T_total"]
        (directory / "full_ft_p0.mlir").write_text(full.to_mlir())
        pbc = cql.compiler.to_pbc(full)
        synthesized["pauli_based_computation"] = logical_counts(pbc)
        (directory / "full_ft_pbc.mlir").write_text(pbc.to_mlir())
    result = {"L": L, "tau": tau, "steps": steps, "algorithm": algorithm,
              "configuration": config,
              "source_sha256": {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                for name in ('hubbard_circuits.py','cudaq_compile.py','model.py')},
              "central_modes": central, "angles": summary, "analog_logical": analog_counts,
              "native_logical": native_counts, "native_schedule": scheduled,
              "fswap_logical":swapped_counts,"fswap_schedule":swapped_schedule,
              "fswap_mode_to_wire":permutation,"fswap_count":swaps,
              "fswap_central_wires":[permutation[q] for q in central],
              "synthesis": synthesized, "compile_wall_seconds": time.monotonic()-start,
              "status": {"physical_routing": "explicit NN Hamiltonian path, all-gate ASAP schedule",
                         "STAR": "compiled analog P0; encoded analog gadget timing remains external",
                         "FT": "actual GridSynth counts; factory and surface timing remain external"}}
    (directory / "result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(f"Finished {directory.name}: G2={scheduled['two_qubit_gates']}, "
          f"depth={scheduled['depth_layers']}, T={synthesized['T_total']}", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=[2,4])
    parser.add_argument("--algorithms", nargs="+", default=["termwise","plaquette"])
    parser.add_argument("--tau", type=float, default=1.)
    parser.add_argument("--steps", type=int, default=4)
    parser.add_argument("--full-synthesis", action="store_true")
    parser.add_argument("--manifest-only", action="store_true",
                        help="Refresh artifact/source hashes after generating plots and diagnostics")
    parser.add_argument("--out", type=Path, default=ROOT/"results/cudaq_v1")
    args = parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    config = json.loads((ROOT/"assumptions.json").read_text())
    for L in ([] if args.manifest_only else args.sizes):
        for algorithm in args.algorithms:
            compile_case(L,args.tau,args.steps,algorithm,config,args.out,args.full_synthesis)
    files = [ROOT/name for name in ("hubbard_circuits.py","cudaq_compile.py","assumptions.json",
             "model.py","compiled_model.py","plot_compiled.py","cudaq_surface_probe.py",
             "cudaq_encoded_reference.py","test_hubbard_circuits.py","test_compiled_model.py",
             "README.md","WORKING_DECISIONS.md","CUDA_Q_RESULTS.md","requirements-cudaq.lock")]
    manifest = {"packages": {name: importlib.metadata.version(name)
                            for name in ("cudaq-logical","cuda-quantum-cu12","numpy","scipy")},
                "configuration": config,
                "sources": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in files if path.exists()},
                "artifacts": {str(path.relative_to(args.out)): hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in sorted(args.out.rglob("*")) if path.is_file() and path.name != "manifest.json"}}
    (args.out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")


if __name__ == "__main__":
    main()
