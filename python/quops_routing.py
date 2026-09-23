#!/usr/bin/env python3
"""Reproducible QUOPS routing census; writes only its new result directory/figure.

One layer is a uniform CNOT matching followed by iid Pauli rotations (paper II B,
full final layer). All physical sites, including unused ones, are mapped because
Cirq 1.4 RouteCQC otherwise restricts routing to the initially occupied subgraph.
Run with -B, OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 on the login node.
"""
from __future__ import annotations

import argparse
from collections import defaultdict, deque
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import socket
import sys
import time

import cirq
import networkx as nx
import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/quops_routing_v1"
WIDTHS = [4, 6, 9, 12, 16, 20, 25, 30, 36, 49, 64, 81, 100, 121, 144]
METRICS = ["swaps_per_layer", "added_cnot_depth_per_layer",
           "fused_cnots_per_layer", "fused_cnot_depth_per_layer"]


def near_square(n):
    cols = math.ceil(math.sqrt(n))
    return math.ceil(n / cols), cols


def snake(rows, cols, offset=(0, 0)):
    return [cirq.GridQubit(r + offset[0], c + offset[1])
            for r in range(rows)
            for c in (range(cols) if r % 2 == 0 else reversed(range(cols)))]


def device(width, area):
    rows, cols = near_square(math.ceil(area * width))
    base_rows, base_cols = near_square(width)
    physical = list(cirq.GridQubit.rect(rows, cols))
    offset = ((rows - base_rows) // 2, (cols - base_cols) // 2)
    # The same compact initial patch is embedded centrally in the larger device.
    occupied = snake(base_rows, base_cols, offset)[:width]
    occupied_set = set(occupied)
    ordered = occupied + [q for q in physical if q not in occupied_set]
    mapping = dict(zip(cirq.LineQubit.range(rows * cols), ordered))
    graph = nx.Graph()
    graph.add_nodes_from(physical)
    graph.add_edges_from((q, z) for q in physical for z in
                         (cirq.GridQubit(q.row + 1, q.col),
                          cirq.GridQubit(q.row, q.col + 1)) if z in graph)
    return rows, cols, graph, mapping, occupied


def ensemble(width, layers, seed, capacity):
    rng = np.random.default_rng(seed)
    logical = cirq.LineQubit.range(capacity)
    moments, stages = [], []
    if capacity > width:
        moments.append(cirq.Moment(cirq.I(q) for q in logical[width:]))
    index = 0
    for _ in range(layers):
        perm = rng.permutation(width)
        pairs = [(int(a), int(b)) for a, b in zip(perm[::2], perm[1::2])]
        cnots, rotations = [], []
        for a, b in pairs:
            cnots.append(cirq.CNOT(logical[a], logical[b]).with_tags(index))
            index += 1
        for q in logical[:width]:
            axis = int(rng.integers(3))
            angle = float(rng.uniform(0, 2 * np.pi))
            rotations.append((cirq.rx, cirq.ry, cirq.rz)[axis](angle)(q)
                             .with_tags(index))
            index += 1
        moments.extend([cirq.Moment(cnots), cirq.Moment(rotations)])
        stages.append((pairs, cnots, rotations))
    return cirq.Circuit(moments), stages


def constructive(stages, rows, cols, initial, occupied):
    """Match pairs to nearby dominoes, then row-column-row permutation routing.

    A C-regular bipartite multigraph between source/destination rows decomposes
    into C perfect matchings, assigning an intermediate column to every token.
    Parallel odd-even sorts in rows, columns, rows need <= 2*C+R swap rounds.
    Transposing the board gives the smaller 2*min(R,C)+max(R,C) bound.
    """
    qlogical = list(sorted(initial))
    position = {q.x: p for q, p in initial.items()}
    inverse = {p: q.x for q, p in initial.items()}
    output = []
    swap_rounds = 0
    # Rows are the longer dimension, so the row sorts use the shorter dimension.
    transpose = cols > rows
    nrow, ncol = (cols, rows) if transpose else (rows, cols)

    def cell(r, c):
        return cirq.GridQubit(c, r) if transpose else cirq.GridQubit(r, c)

    def coordinates(q):
        return (q.col, q.row) if transpose else (q.row, q.col)

    def dist(q, z):
        return abs(q.row - z.row) + abs(q.col - z.col)

    for pairs, cnots, rotations in stages:
        slots = list(zip(occupied[::2], occupied[1::2]))
        costs, reverse = np.zeros((len(pairs), len(slots))), {}
        for i, (a, b) in enumerate(pairs):
            for j, (x, y) in enumerate(slots):
                straight = dist(position[a], x) + dist(position[b], y)
                flipped = dist(position[a], y) + dist(position[b], x)
                costs[i, j] = min(straight, flipped)
                reverse[i, j] = flipped < straight
        target = {}
        for i, j in zip(*linear_sum_assignment(costs)):
            a, b = pairs[i]
            x, y = slots[j]
            target[a], target[b] = (y, x) if reverse[i, j] else (x, y)
        if len(occupied) % 2:
            unpaired = next(i for i in range(len(occupied)) if i not in target)
            target[unpaired] = occupied[-1]
        free = set(inverse) - set(target.values())
        # Keep dummy tokens still whenever possible, then minimize their travel.
        for token in range(len(occupied), rows * cols):
            if position[token] in free:
                target[token] = position[token]
                free.remove(position[token])
        remaining = [i for i in range(rows * cols) if i not in target]
        destinations = sorted(free)
        if remaining:
            cost = [[dist(position[i], q) for q in destinations] for i in remaining]
            for i, j in zip(*linear_sum_assignment(cost)):
                target[remaining[i]] = destinations[j]
        target_rc = {i: coordinates(q) for i, q in target.items()}
        parallel_edges = defaultdict(list)
        for token, q in position.items():
            parallel_edges[coordinates(q)[0], target_rc[token][0]].append(token)
        intermediate = {}
        for col in range(ncol):
            graph = nx.Graph()
            graph.add_nodes_from(range(2 * nrow))
            graph.add_edges_from((a, nrow + b) for (a, b), tokens in
                                 sorted(parallel_edges.items()) if tokens)
            match = nx.bipartite.maximum_matching(graph, top_nodes=range(nrow))
            assert len(match) == 2 * nrow
            for a in range(nrow):
                token = parallel_edges[a, match[a] - nrow].pop()
                intermediate[token] = col

        def sort_lines(columnwise, key):
            nonlocal swap_rounds
            length, count = (nrow, ncol) if columnwise else (ncol, nrow)
            for phase in range(length):
                swaps = []
                for line in range(count):
                    for start in range(phase % 2, length - 1, 2):
                        x = cell(start, line) if columnwise else cell(line, start)
                        y = cell(start + 1, line) if columnwise else cell(line, start + 1)
                        a, b = inverse[x], inverse[y]
                        if key(a) > key(b):
                            swaps.append(cirq.SWAP(x, y))
                            inverse[x], inverse[y] = b, a
                            position[a], position[b] = y, x
                if swaps:
                    output.append(cirq.Moment(swaps))
                    swap_rounds += 1

        sort_lines(False, intermediate.__getitem__)
        sort_lines(True, lambda token: target_rc[token][0])
        sort_lines(False, lambda token: target_rc[token][1])
        assert position == target
        mapping = {q: position[q.x] for q in qlogical}
        output.append(cirq.Moment(op.transform_qubits(mapping) for op in cnots))
        output.append(cirq.Moment(op.transform_qubits(mapping) for op in rotations))
    return cirq.Circuit(output), {initial[q]: position[q.x] for q in qlogical}, swap_rounds


# Minimal CNOT words for the six reversible linear maps on two bits. Fusion is
# exact for maximal same-edge blocks with no intervening operation on either end.
IDENTITY = (0, 1, 2, 3)
CX01 = (0, 1, 3, 2)
CX10 = (0, 3, 2, 1)
SWAP = (0, 2, 1, 3)


def compose(new, old):
    return tuple(new[x] for x in old)


def minimal_words():
    words, queue = {IDENTITY: ()}, deque([IDENTITY])
    while queue:
        old = queue.popleft()
        for direction, transform in enumerate((CX01, CX10)):
            new = compose(transform, old)
            if new not in words:
                words[new] = words[old] + (direction,)
                queue.append(new)
    assert len(words) == 6 and max(map(len, words.values())) == 3
    return words


WORDS = minimal_words()


def fused_operations(operations):
    blocks, last = [], {}
    for op in operations:
        if cirq.num_qubits(op) == 1:
            if op.gate == cirq.I:
                continue
            q = op.qubits[0]
            i = last.get(q, -1)
            if i >= 0 and isinstance(blocks[i], list):
                blocks[i][2].append(op)
            else:
                last[q] = len(blocks)
                blocks.append(op)
            continue
        x, y = sorted(op.qubits)
        i = last.get(x, -1)
        if i >= 0 and i == last.get(y, -2) and isinstance(blocks[i], list):
            blocks[i][2].append(op)
        else:
            i = len(blocks)
            blocks.append([x, y, [op]])
        last[x] = last[y] = i
    fused = []
    for block in blocks:
        if not isinstance(block, list):
            fused.append(block)
        else:
            x, y, local_ops = block
            entanglers = [op for op in local_ops if cirq.num_qubits(op) == 2]
            if len(entanglers) == 1:
                for op in local_ops:
                    if op.gate == cirq.SWAP:
                        fused.extend([cirq.CNOT(x, y), cirq.CNOT(y, x), cirq.CNOT(x, y)])
                    else:
                        fused.append(op)
            elif len(entanglers) == len(local_ops):
                transform = IDENTITY
                for op in local_ops:
                    g = SWAP if op.gate == cirq.SWAP else (
                        CX01 if op.qubits[0] == x else CX10)
                    transform = compose(g, transform)
                fused.extend(cirq.CNOT(x, y) if d == 0 else cirq.CNOT(y, x)
                             for d in WORDS[transform])
            else:
                # KAK includes arbitrary intervening 1q rotations; full CZs are
                # CNOT-equivalent up to 1q gates. No partial-CZ cost discount.
                unitary = cirq.Circuit(local_ops).unitary(qubit_order=[x, y])
                compiled = cirq.two_qubit_matrix_to_cz_operations(
                    x, y, unitary, allow_partial_czs=False, atol=1e-10)
                assert all(cirq.num_qubits(op) == 1 or op.gate == cirq.CZ for op in compiled)
                fused.extend(compiled)
    return fused


def two_qubit_depth(operations, swap_weight=1):
    clocks = defaultdict(int)
    for op in operations:
        if cirq.num_qubits(op) == 2:
            finish = max(clocks[q] for q in op.qubits) + (
                swap_weight if op.gate == cirq.SWAP else 1)
            for q in op.qubits:
                clocks[q] = finish
    return max(clocks.values(), default=0)


def validate_route(original, routed, initial, final):
    """Check every original gate and its causal order, adjacency, final mapping."""
    expected, order = {}, defaultdict(list)
    for op in original.all_operations():
        if op.tags:
            index = op.tags[0]
            expected[index] = op
            for q in op.qubits:
                order[q].append(index)
    at = {p: q for q, p in initial.items()}
    seen, actual_order = set(), defaultdict(list)
    for op in routed.all_operations():
        if cirq.num_qubits(op) == 2:
            a, b = op.qubits
            assert abs(a.row - b.row) + abs(a.col - b.col) == 1
        if op.gate == cirq.SWAP:
            a, b = op.qubits
            at[a], at[b] = at[b], at[a]
        elif op.gate != cirq.I:
            assert len(op.tags) == 1
            index = op.tags[0]
            logical = expected[index]
            assert index not in seen and op.gate == logical.gate
            assert tuple(at[q] for q in op.qubits) == logical.qubits
            seen.add(index)
            for q in logical.qubits:
                actual_order[q].append(index)
    assert seen == set(expected) and dict(actual_order) == dict(order)
    assert all(at[final[p]] == q for q, p in initial.items())


def measure(routed, width, layers):
    ops = list(routed.all_operations())
    swaps = sum(op.gate == cirq.SWAP for op in ops)
    cnots = sum(op.gate == cirq.CNOT for op in ops)
    assert cnots == layers * (width // 2)
    fused = fused_operations(ops)
    fused_count = sum(cirq.num_qubits(op) == 2 for op in fused)
    primitive_depth = two_qubit_depth(ops)
    cnot_depth = two_qubit_depth(ops, swap_weight=3)
    fused_depth = two_qubit_depth(fused)
    return dict(swaps=swaps, original_cnots=cnots,
                unfused_cnots=cnots + 3 * swaps, fused_cnots=fused_count,
                primitive_2q_depth=primitive_depth, cnot_depth=cnot_depth,
                fused_cnot_depth=fused_depth,
                swaps_per_layer=swaps / layers,
                swaps_per_cnot=swaps / cnots,
                added_swap_unit_depth_per_layer=primitive_depth / layers - 1,
                added_cnot_depth_per_layer=cnot_depth / layers - 1,
                fused_cnots_per_layer=fused_count / layers,
                fused_cnot_depth_per_layer=fused_depth / layers,
                added_fused_cnot_depth_per_layer=fused_depth / layers - 1,
                fusion_saving_fraction=1 - fused_count / (cnots + 3 * swaps))


def small_checks():
    """Statevector checks complement the all-circuit symbolic routing audits."""
    for transform, word in WORDS.items():
        result = IDENTITY
        for d in word:
            result = compose((CX01, CX10)[d], result)
        assert result == transform
    # Explicitly exercise fusion across disjoint gates and its 1q barriers.
    qa, qb, qc = cirq.LineQubit.range(3)
    for middle in ([], [cirq.rx(.31)(qa)], [cirq.ry(.17)(qc)]):
        circuit = cirq.Circuit([cirq.CNOT(qa, qb), *middle, cirq.SWAP(qa, qb)])
        fused = cirq.Circuit(fused_operations(circuit.all_operations()))
        cirq.testing.assert_allclose_up_to_global_phase(
            circuit.unitary(qubit_order=[qa, qb, qc]),
            fused.unitary(qubit_order=[qa, qb, qc]), atol=1e-12)
    checked = 0
    for width in (4, 5):
        for area in (1, 2):
            rows, cols, graph, initial, occupied = device(width, area)
            original, stages = ensemble(width, 3, 901 + width, rows * cols)
            cqc, _, final_cqc = cirq.RouteCQC(graph).route_circuit(
                original, initial_mapper=cirq.HardCodedInitialMapper(initial))
            construct, final_construct, _ = constructive(
                stages, rows, cols, initial, occupied)
            physical = sorted(initial.values())
            reference = original.transform_qubits(initial)
            rng = np.random.default_rng(871)
            state = rng.normal(size=2**len(physical)) + 1j * rng.normal(size=2**len(physical))
            state /= np.linalg.norm(state)
            sim = cirq.Simulator(dtype=np.complex128)
            target = sim.simulate(reference, qubit_order=physical, initial_state=state).final_state_vector
            for routed, final in ((cqc, final_cqc), (construct, final_construct)):
                validate_route(original, routed, initial, final)
                for candidate in (routed, cirq.Circuit(fused_operations(routed.all_operations()))):
                    result = sim.simulate(candidate, qubit_order=[final[q] for q in physical],
                                          initial_state=state.reshape([2] * len(physical)).transpose(
                                              [physical.index(final[q]) for q in physical]).reshape(-1)
                                          ).final_state_vector
                    cirq.testing.assert_allclose_up_to_global_phase(result, target, atol=1e-10)
                    checked += 1
    return {"statevector_comparisons": checked, "fusion_group_elements": len(WORDS),
            "all_circuits": "gate identities, causal order, adjacency, final token permutation"}


def run_case(task):
    width, layers, area, sample, seed, output, deadline, old = task
    if time.monotonic() > deadline:
        return None
    path = Path(output) / "cases" / f"w{width}_D{layers}_a{area}_s{sample}_seed{seed}.json"
    if path.exists():
        saved = json.loads(path.read_text())
        assert saved["status"] == "complete"
        return saved["records"]
    rows, cols, graph, initial, occupied = device(width, area)
    original, stages = ensemble(width, layers, seed, rows * cols)
    circuit_hash = hashlib.sha256(str(original).encode()).hexdigest()
    result = []
    for name in ("cqc", "constructive"):
        if name in old:
            assert int(old[name]["seed"]) == seed
            result.append(old[name])
            continue
        t0 = time.monotonic()
        if name == "cqc":
            routed, returned_initial, final = cirq.RouteCQC(graph).route_circuit(
                original, initial_mapper=cirq.HardCodedInitialMapper(initial),
                lookahead_radius=8)
            assert returned_initial == initial
            swap_rounds = -1
        else:
            routed, final, swap_rounds = constructive(stages, rows, cols, initial, occupied)
            if area == 2:
                # Spare area cannot worsen the achievable cost: include the
                # compact network embedded in the larger array as a candidate.
                br, bc, _, compact_initial, compact_occupied = device(width, 1)
                compact, compact_final, compact_rounds = constructive(
                    stages, br, bc, compact_initial, compact_occupied)
                dr, dc = occupied[0].row, occupied[0].col
                shift = lambda q: cirq.GridQubit(q.row + dr, q.col + dc)
                compact = compact.transform_qubits(shift)
                compact_final = {shift(q): shift(z) for q, z in compact_final.items()}
                compact_final.update({q: q for q in initial.values() if q not in compact_final})
                full_metrics, compact_metrics = measure(routed, width, layers), measure(compact, width, layers)
                if (compact_metrics["fused_cnots"], compact_metrics["fused_cnot_depth"]) <= (
                        full_metrics["fused_cnots"], full_metrics["fused_cnot_depth"]):
                    routed, final, swap_rounds = compact, compact_final, compact_rounds
        routing_seconds = time.monotonic() - t0
        validate_route(original, routed, initial, final)
        result.append(dict(width=width, layers=layers, sample=sample, seed=seed,
                           area_multiplier=area, rows=rows, cols=cols,
                           physical_qubits=rows * cols, router=name,
                           circuit_sha256=circuit_hash, routing_seconds=routing_seconds,
                           constructive_swap_rounds=swap_rounds,
                           **measure(routed, width, layers)))
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(dict(status="complete", records=result,
                                         host=socket.gethostname(), pid=os.getpid()), indent=2) + "\n")
    temporary.replace(path)
    return result


def census(args):
    output = Path(args.output)
    (output / "cases").mkdir(parents=True, exist_ok=True)
    raw_path = output / "raw.csv"
    existing = defaultdict(dict)
    if raw_path.exists():
        with raw_path.open() as f:
            for r in csv.DictReader(f):
                key = int(r["width"]), int(r["layers"]), int(r["area_multiplier"]), int(r["sample"])
                existing[key][r["router"]] = r
    start = time.monotonic()
    validations = small_checks()
    tasks = []
    for width in args.widths:
        for layers in args.layers:
            for area in (1, 2):
                for sample in range(args.samples):
                    seed = args.seed + 100000 * width + 1000 * layers + sample
                    tasks.append((width, layers, area, sample, seed, str(output),
                                  start + args.max_seconds, existing[width, layers, area, sample]))
    # Expensive cases first reduce the final straggler; identity and seeds are fixed.
    tasks.sort(key=lambda t: (-t[0]**3 * t[1] * t[2], t[3]))
    collected, completed = [], 0
    if args.workers == 1:
        iterator = map(run_case, tasks)
        pool = None
    else:
        pool = ProcessPoolExecutor(max_workers=args.workers)
        iterator = (f.result() for f in as_completed([pool.submit(run_case, t) for t in tasks]))
    try:
        for result in iterator:
            if result is not None:
                collected.extend(result)
                completed += 1
                if completed % 20 == 0 or completed == len(tasks):
                    print(f"{completed}/{len(tasks)} paired cases; "
                          f"elapsed {time.monotonic()-start:.1f}s", flush=True)
    finally:
        if pool is not None:
            pool.shutdown()
    collected.sort(key=lambda r: (int(r["width"]), int(r["layers"]),
                                  int(r["area_multiplier"]), int(r["sample"]), r["router"]))
    if collected:
        temporary = raw_path.with_suffix(".tmp")
        with temporary.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(collected[0]))
            writer.writeheader()
            writer.writerows(collected)
        temporary.replace(raw_path)
    validations.update(elapsed_seconds_this_run=time.monotonic() - start,
                       host=socket.gethostname(), workers=args.workers,
                       expected_pairs=len(tasks), completed_pairs=completed,
                       complete=completed == len(tasks), command=sys.argv,
                       versions={"cirq": cirq.__version__, "networkx": nx.__version__},
                       source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (output / "validation.json").write_text(json.dumps(validations, indent=2) + "\n")
    if completed != len(tasks):
        raise SystemExit("Time budget reached; atomic cases are resumable.")


def load_raw(output):
    with (output / "raw.csv").open() as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        for key in row:
            if key not in ("router", "circuit_sha256"):
                row[key] = float(row[key])
    return rows


def fit_power(groups, metric, rng, bootstraps, minimum_width):
    widths = sorted(w for w in groups if w >= minimum_width)
    values = [np.array([r[metric] for r in groups[w]]) for w in widths]
    means = np.array([np.mean(v) for v in values])
    x = np.log(widths)

    def solve(y):
        b, loga = np.polyfit(x, np.log(y), 1)
        return np.array([np.exp(loga), b])

    estimate = solve(means)
    # Resample independent circuits within each width AND depth stratum. Equal
    # numbers at each D make the overall arithmetic mean equally depth-weighted.
    strata = [[np.array([r[metric] for r in groups[w] if r["layers"] == d])
               for d in sorted({r["layers"] for r in groups[w]})] for w in widths]
    draws = np.array([solve([np.mean([np.mean(rng.choice(v, len(v), replace=True))
                                     for v in group]) for group in strata])
                      for _ in range(bootstraps)])
    prediction = estimate[0] * np.array(widths)**estimate[1]
    residuals = means / prediction
    return dict(form="a * width ** b", a=estimate[0], b=estimate[1],
                a_ci95=np.quantile(draws[:, 0], [.025, .975]).tolist(),
                b_ci95=np.quantile(draws[:, 1], [.025, .975]).tolist(),
                min_width=minimum_width, widths=widths,
                max_relative_fit_residual=float(np.max(np.abs(residuals - 1))),
                rms_log_residual=float(np.sqrt(np.mean(np.log(residuals)**2))))


def summarize(args):
    output = Path(args.output)
    rows = load_raw(output)
    rng = np.random.default_rng(918271)
    fits, grouped, comparisons, means = {}, {}, {}, []
    for area in (1, 2):
        for router in ("cqc", "constructive"):
            name = f"{router}_area{area}"
            groups = defaultdict(list)
            for row in rows:
                if row["router"] == router and row["area_multiplier"] == area:
                    groups[int(row["width"])].append(row)
            grouped[name] = groups
            fits[name] = {m: fit_power(groups, m, rng, args.bootstraps, 16) for m in METRICS}
            fits[name]["all_width_sensitivity"] = {
                m: fit_power(groups, m, rng, args.bootstraps, 4) for m in METRICS}
            fits[name]["depth_sensitivity"] = {}
            for d in sorted({r["layers"] for r in rows}):
                sub = {w: [r for r in group if r["layers"] == d]
                       for w, group in groups.items()}
                fits[name]["depth_sensitivity"][str(int(d))] = {
                    m: fit_power(sub, m, rng, args.bootstraps, 16) for m in METRICS}
            for w, group in sorted(groups.items()):
                entry = dict(router=router, area_multiplier=area, width=w,
                             n=len(group), physical_qubits=int(group[0]["physical_qubits"]))
                for m in METRICS + ["swaps_per_cnot", "fusion_saving_fraction",
                                     "added_swap_unit_depth_per_layer"]:
                    v = np.array([r[m] for r in group])
                    entry[m] = float(v.mean())
                    entry[m + "_se"] = float(v.std(ddof=1) / np.sqrt(len(v)))
                means.append(entry)
        indexed = {(r["width"], r["layers"], r["sample"], r["router"]): r
                   for r in rows if r["area_multiplier"] == area}
        comparisons[f"area{area}"] = {}
        for m in METRICS:
            ratios = [r[m] / indexed[(w, d, s, "cqc")][m]
                      for (w, d, s, router), r in indexed.items()
                      if router == "constructive" and w >= 16]
            comparisons[f"area{area}"][m] = dict(
                median_constructive_over_cqc=float(np.median(ratios)),
                p05_p95=np.quantile(ratios, [.05, .95]).tolist())
    result = dict(
        schema=1, paper="https://arxiv.org/abs/2609.12146", paper_section="II B",
        ensemble="uniform disjoint CNOT matching then iid R_P(theta), P=X/Y/Z, theta uniform [0,2pi); zeta=1",
        depth_convention="D benchmark layers = paper depth d/2; no final mapping restoration",
        geometry="cols=ceil(sqrt(area*w)), rows=ceil(area*w/cols); same centered compact initial patch",
        routers={"cqc": "Cirq RouteCQC 1.4.0, fixed initial map, lookahead_radius=8, all sites mapped",
                 "constructive": "minimum-distance pair-to-domino assignment; bipartite edge coloring; row-column-row odd-even sorting; 2x area also admits embedded compact routing and selects fewer fused CNOTs, then lower depth"},
        counting="1 SWAP=3 CNOT; 2q critical-path depth ignores 1q duration, respects gate dependencies; maximal same-edge blocks including intervening 1q gates are fused via exact linear synthesis or KAK into full CZs (one CNOT-equivalent each), no partial-CZ discount",
        uncertainty="95% stratified circuit bootstrap; excludes compiler/model/extrapolation uncertainty",
        fit_method="unweighted log-linear fit of arithmetic means per width, equally mixing sampled D; principal width>=16",
        asymptotic_constructive_bounds="<=2*min(rows,cols)+max(rows,cols) SWAP rounds and floor(Q/2) times that many SWAPs per benchmark layer",
        versions=dict(python=platform.python_version(), cirq=cirq.__version__, networkx=nx.__version__, numpy=np.__version__),
        records=len(rows), fits=fits, router_comparison=comparisons, means=means,
        sha256={"raw.csv": hashlib.sha256((output / "raw.csv").read_bytes()).hexdigest(),
                "python/quops_routing.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "refs/2609.12146v1.pdf": hashlib.sha256((ROOT / "refs/2609.12146v1.pdf").read_bytes()).hexdigest()})
    (output / "fit.json").write_text(json.dumps(result, indent=2) + "\n")
    os.environ.setdefault("MPLCONFIGDIR", str(output / ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7), constrained_layout=True)
    plotted = ["swaps_per_layer", "added_cnot_depth_per_layer", "fused_cnots_per_layer"]
    titles = ["Added SWAPs / QUOPS layer", "Added CNOT depth / QUOPS layer",
              "Fused CNOTs / QUOPS layer"]
    for ax, metric, title in zip(axes, plotted, titles):
        for area, style in ((1, "-"), (2, "--")):
            for router, color in (("cqc", "#2066ac"), ("constructive", "#c56721")):
                group = [r for r in means if r["router"] == router and r["area_multiplier"] == area]
                x = np.array([r["width"] for r in group])
                y = np.array([r[metric] for r in group])
                err = np.array([1.96 * r[metric + "_se"] for r in group])
                label = f"{'RouteCQC' if router == 'cqc' else 'Constructive'}, {area}× area"
                ax.errorbar(x, y, yerr=err, color=color, linestyle=style, marker="o",
                            markersize=3, linewidth=1.3, label=label, capsize=2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Logical width w")
        ax.grid(alpha=.18, which="both")
        ax.set_xticks([4, 16, 36, 64, 144], [4, 16, 36, 64, 144])
    axes[0].legend(fontsize=8, loc="upper left")
    fig.suptitle("Nearest-neighbor QUOPS routing • full random matchings • D = 5, 10, 20", fontsize=11)
    figure = Path(args.figure)
    figure.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure, dpi=200)
    plt.close(fig)
    print(json.dumps({k: {m: [round(v[m]["a"], 4), round(v[m]["b"], 4)]
                               for m in METRICS} for k, v in fits.items()}, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=["run", "summarize", "check"])
    p.add_argument("--output", default=str(DEFAULT_OUTPUT))
    p.add_argument("--figure", default=str(ROOT / "figures/quops_routing.png"))
    p.add_argument("--widths", nargs="+", type=int, default=WIDTHS)
    p.add_argument("--layers", nargs="+", type=int, default=[5, 10, 20])
    p.add_argument("--samples", type=int, default=20)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--seed", type=int, default=260912146)
    p.add_argument("--max-seconds", type=float, default=2100)
    p.add_argument("--bootstraps", type=int, default=1000)
    args = p.parse_args()
    if args.command == "run":
        census(args)
    elif args.command == "summarize":
        summarize(args)
    else:
        print(json.dumps(small_checks(), indent=2))


if __name__ == "__main__":
    main()
