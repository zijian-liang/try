"""Circuit-distance bounds, with certified and heuristic results kept separate.

All weights count noisy physical locations, not propagated data-qubit weight.
The complete, *undecomposed* Stim DEM is used.  Its unit-weight column problem
is a relaxation of the physical fault problem, so every lower bound is sound.
An upper bound is accepted only after every selected column is mapped back to
an actual circuit error.  The supplied circuits use Pauli-closed depolarizing
channels: multiple selected Paulis at one location can be combined into one
allowed Pauli, never increasing the number of faulty locations.

Usage: python distance_audit.py --codes w7 bb72 surface7 --seconds 120
Longer optional SAT certification: --sat-seconds 3600 --max-prove 6
This does not estimate circuit distance from Monte Carlo slopes.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import threading
import time
from typing import Any

import stim


def circuit_hash(circuit: stim.Circuit) -> str:
    return hashlib.sha256(str(circuit).encode()).hexdigest()


def _indices(mask: int):
    while mask:
        b = mask & -mask
        yield b.bit_length() - 1
        mask ^= b


def dem_columns(dem: stim.DetectorErrorModel) -> list[tuple[int, int]]:
    """Distinct (detector mask, observable mask) for supported single faults."""
    result = set()
    for inst in dem.flattened():
        if inst.type != "error" or inst.args_copy()[0] == 0:
            continue
        det = obs = 0
        for target in inst.targets_copy():
            if target.is_separator():
                raise ValueError("Distance audit requires decompose_errors=False")
            if target.is_relative_detector_id():
                det ^= 1 << target.val
            elif target.is_logical_observable_id():
                obs ^= 1 << target.val
        if det or obs:
            result.add((det, obs))
    return sorted(result)


def columns_dem(columns: list[tuple[int, int]]) -> stim.DetectorErrorModel:
    result = stim.DetectorErrorModel()
    for det, obs in columns:
        targets = [stim.target_relative_detector_id(i) for i in _indices(det)]
        targets += [stim.target_logical_observable_id(i) for i in _indices(obs)]
        result.append("error", 1, targets)
    return result


def physical_witness(circuit: stim.Circuit, columns: list[tuple[int, int]]) -> dict:
    det = obs = 0
    for d, o in columns:
        det ^= d
        obs ^= o
    if det or not obs:
        raise ValueError("Invalid witness: must have zero detectors and nonzero logical action")
    explanation = circuit.explain_detector_error_model_errors(
        dem_filter=columns_dem(columns), reduce_to_one_representative_error=True
    )
    if len(explanation) != len(columns) or any(not e.circuit_error_locations for e in explanation):
        raise ValueError("A DEM witness column has no matching physical single-fault location")
    locations = []
    for error in explanation:
        loc = error.circuit_error_locations[0]
        locations.append({"tick_offset": loc.tick_offset,
                          "target_range_start": loc.instruction_targets.target_range_start,
                          "target_range_end": loc.instruction_targets.target_range_end,
                          "stack": [[f.instruction_offset, f.iteration_index] for f in loc.stack_frames]})
    return {
        "weight_upper_bound": len(columns),
        "logical_observables_flipped": list(_indices(obs)),
        "detector_parity_zero": True,
        "dem_errors": str(columns_dem(columns)),
        "physical_errors": [str(e) for e in explanation],
        "physical_location_keys": locations,
        "distinct_physical_locations": len({json.dumps(v, sort_keys=True) for v in locations}),
    }


def static_logical_physical_witness(circuit, qubits, sector):
    """Realize a static logical by data-only faults after each last noisy CX.

    The constructed single-fault responses must exist in the original DEM;
    their combined response is checked again by physical_witness.  Thus this
    is a physical witness, not an assumed d_circuit <= static_distance rule.
    """
    operations = list(circuit.flattened())
    last = {}
    for i, op in enumerate(operations):
        if op.name == "DEPOLARIZE2" and op.gate_args_copy()[0] > 0:
            for t in op.targets_copy():
                if t.value in qubits:
                    last[t.value] = i
    if len(last) != len(qubits):
        raise ValueError("Static witness includes data without a supported final CNOT fault")
    clean_parts = []
    for op in operations:
        part = stim.Circuit()
        part.append(op)
        clean_parts.append(part.without_noise())
    columns = []
    for q in qubits:
        single = stim.Circuit()
        for i, part in enumerate(clean_parts):
            single += part
            if i == last[q]:
                single.append(sector + "_ERROR", [q], .125)
        response = dem_columns(single.detector_error_model(decompose_errors=False))
        if len(response) != 1:
            raise ValueError("Static witness qubit has no unique nontrivial fault response")
        columns += response
    return physical_witness(circuit, columns)


def small_weight_audit(columns: list[tuple[int, int]], *, triple_seconds: float = 30) -> dict:
    """Exactly exclude weights 1,2 and, if completed, 3.

    For a zero-syndrome triple, sort detector degrees a>=b>=c.  Its first two
    columns overlap on (a+b-c)/2 >= ceil(a/2) detectors.  The inverted index
    therefore enumerates every possible logical triple without checking all
    unrelated column pairs.  Stopping early NEVER certifies weight-3 absence.
    """
    seen: dict[int, tuple[int, int]] = {}
    for idx, (det, obs) in enumerate(columns):
        if det == 0 and obs:
            return {"lower_bound": 1, "witness_indices": [idx], "method": "exact_weight_1"}
    for idx, (det, obs) in enumerate(columns):
        if det in seen and seen[det][1] != obs:
            return {"lower_bound": 2, "witness_indices": [seen[det][0], idx], "method": "exact_weights_1_2"}
        seen[det] = (idx, obs)
    result = {"lower_bound": 3, "method": "exact_weights_1_2", "triple_complete": False}
    if triple_seconds <= 0:
        return result
    start = time.monotonic()
    degrees = [d.bit_count() for d, _ in columns]
    incident: dict[int, list[int]] = defaultdict(list)
    for idx, (det, _) in enumerate(columns):
        for i in _indices(det):
            incident[i].append(idx)
    for ai in sorted(range(len(columns)), key=lambda i: degrees[i], reverse=True):
        if time.monotonic() - start >= triple_seconds:
            result["triple_seconds"] = time.monotonic() - start
            return result
        ad, ao = columns[ai]
        overlap: dict[int, int] = defaultdict(int)
        for i in _indices(ad):
            for bi in incident[i]:
                if bi != ai and degrees[bi] <= degrees[ai]:
                    overlap[bi] += 1
        threshold = (degrees[ai] + 1) // 2
        for bi, count in overlap.items():
            if count < threshold:
                continue
            bd, bo = columns[bi]
            match = seen.get(ad ^ bd)
            if match is not None:
                ci, co = match
                if ci not in (ai, bi) and degrees[ci] <= degrees[bi] and (ao ^ bo ^ co):
                    result.update(witness_indices=[ai, bi, ci], method="exact_weights_1_2_3")
                    return result
    result.update(lower_bound=4, method="exact_weights_1_2_3", triple_complete=True,
                  triple_seconds=time.monotonic() - start)
    return result


def _project_columns(columns, detector_ids, observable_ids):
    detmask = sum(1 << i for i in detector_ids)
    obsmask = sum(1 << i for i in observable_ids)
    return sorted({(d & detmask, o & obsmask) for d, o in columns if (d & detmask) or (o & obsmask)})


def exact_four_audit(columns, *, seconds=60):
    """Meet-in-the-middle exact exclusion of all weight-4 logical histories.

    Run small_weight_audit first: this routine assumes no weight-1/2/3
    logical.  A dictionary retains complete integer syndrome masks, so hash
    collisions cannot invalidate the proof.  About four million entries per
    CSS sector are typical for this campaign; memory can approach 1 GB.
    """
    started = time.monotonic()
    pairs = {}
    for i, (di, oi) in enumerate(columns):
        if time.monotonic() - started >= seconds:
            return {"lower_bound": 4, "complete": False, "method": "exact_weight4_pair_enumeration_timeout"}
        for j in range(i):
            dj, oj = columns[j]
            d, o = di ^ dj, oi ^ oj
            previous = pairs.get(d)
            if previous is not None:
                if previous[0] != o:
                    return {"lower_bound": 4, "complete": True,
                            "witness_indices": [previous[1], previous[2], i, j],
                            "method": "exact_weight4_pair_enumeration"}
            else:
                pairs[d] = (o, i, j)
    return {"lower_bound": 5, "complete": True, "method": "exact_weight4_pair_enumeration",
            "pair_syndromes": len(pairs), "seconds": time.monotonic() - started}


def projected_four_audit(columns, metadata, *, seconds=60):
    dg, og = metadata.get("detector_groups", {}), metadata.get("observable_groups", {})
    if set(dg) != {"X", "Z"} or set(og) != {"X", "Z"}:
        return None
    out, physical_columns = {}, None
    for sector in ("X", "Z"):
        check = {"X": "Z", "Z": "X"}[sector]
        projected = _project_columns(columns, dg[check], og[sector])
        small = small_weight_audit(projected, triple_seconds=seconds)
        indices = small.pop("witness_indices", None)
        result = small
        if result["lower_bound"] >= 4:
            result = exact_four_audit(projected, seconds=seconds)
            indices = result.pop("witness_indices", None)
        out[sector] = result
        if indices is not None:
            # A projected upper bound alone is NOT a full physical upper
            # bound.  Require pure-sector columns present verbatim in the
            # undecomposed full circuit DEM before accepting this lift.
            full = set(columns)
            selected = [projected[i] for i in indices]
            if all(col in full for col in selected):
                physical_columns = selected
    return {"lower_bound": min(v["lower_bound"] for v in out.values()),
            "method": "exact_complete_CSS_projections_through_weight4",
            "sectors": out, "physical_columns": physical_columns}


def exact_five_audit(columns, *, seconds=120):
    """Exclude weight 5, ASSUMING weights <=4 are already excluded.

    A minimum logical fault set is connected in the graph joining fault
    columns with overlapping detector support: otherwise one connected
    component is a smaller zero-syndrome logical.  Every connected graph
    with five vertices contains a connected induced triple.  Enumerate all
    connected triples and match their syndromes against ALL fault pairs.
    This is exhaustive under the stated lower-bound assumption, not a
    locality cutoff.  Integer syndrome keys are exact.
    """
    started = time.monotonic()
    pairs = {}
    for i, (di, oi) in enumerate(columns):
        for j in range(i):
            dj, oj = columns[j]
            d, o = di ^ dj, oi ^ oj
            previous = pairs.get(d)
            if previous is not None and previous[0] != o:
                raise ValueError("Weight-4 logical found: exact_five_audit precondition is false")
            if previous is None:
                pairs[d] = (o, i, j)
    incident = defaultdict(list)
    for i, (d, _) in enumerate(columns):
        if not d:
            raise ValueError("Zero-detector column violates weight>=5 precondition")
        for bit in _indices(d):
            incident[bit].append(i)
    triples = 0
    for i, (di, oi) in enumerate(columns):
        if time.monotonic() - started > seconds:
            return {"lower_bound": 5, "complete": False, "method": "exact_connected_weight5_timeout",
                    "triples_checked": triples, "seconds": time.monotonic() - started}
        adjacent = set()
        for bit in _indices(di):
            adjacent.update(incident[bit])
        adjacent.discard(i)
        adjacent = sorted(adjacent)
        for jj, j in enumerate(adjacent):
            dj, oj = columns[j]
            first_d, first_o = di ^ dj, oi ^ oj
            for k in adjacent[:jj]:
                # A triangle has three possible centers; choose the smallest
                # center.  A path has one center and is never removed.
                dk, ok = columns[k]
                if (dj & dk) and i > min(j, k):
                    continue
                triples += 1
                previous = pairs.get(first_d ^ dk)
                if previous is not None and previous[0] != first_o ^ ok:
                    witness = [i, j, k, previous[1], previous[2]]
                    if len(set(witness)) != 5:
                        raise ValueError("Lower-weight logical violates weight>=5 precondition")
                    return {"lower_bound": 5, "complete": True, "method": "exact_connected_weight5",
                            "witness_indices": witness, "triples_checked": triples,
                            "seconds": time.monotonic() - started}
    return {"lower_bound": 6, "complete": True, "method": "exact_connected_weight5",
            "assumption": "complete exclusion of all logical fault sets with weight<=4",
            "triples_checked": triples, "pair_syndromes": len(pairs),
            "seconds": time.monotonic() - started}


def projected_five_audit(columns, metadata, weight4_certificate, *, seconds=120):
    """Apply the connected-triple proof independently to covering CSS projections."""
    if any(weight4_certificate["sectors"][s]["lower_bound"] < 5 for s in ("X", "Z")):
        return None
    results, physical_columns = {}, None
    for sector, check in (("X", "Z"), ("Z", "X")):
        p = _project_columns(columns, metadata["detector_groups"][check], metadata["observable_groups"][sector])
        result = exact_five_audit(p, seconds=seconds)
        ids = result.pop("witness_indices", None)
        if ids is not None:
            selected = [p[i] for i in ids]
            if all(col in set(columns) for col in selected):
                physical_columns = selected
        results[sector] = result
    return {"lower_bound": min(r["lower_bound"] for r in results.values()),
            "method": "exact_complete_CSS_projections_through_weight5",
            "sectors": results, "physical_columns": physical_columns}


def sector_graph_lower_bound(columns, metadata) -> dict | None:
    """Lower bound via a covering set of detector/observable projections.

    Each nontrivial full fault history remains nontrivial in at least one
    observable projection.  Discarding the other detectors only relaxes the
    zero-record condition.  Hence min(projected exact distances) is a lower
    bound; no independence of X/Z faults is assumed here.
    """
    dg = metadata.get("detector_groups", {})
    og = metadata.get("observable_groups", {})
    if not dg or not og:
        return None
    all_obs = set(i for _, o in columns for i in _indices(o))
    covered = set(i for indices in og.values() for i in indices)
    if not all_obs <= covered:
        return None
    out = {}
    for name, obs_ids in og.items():
        # circuits.py names observable groups by residual ERROR component:
        # X observables are logical-Z/ref-Z parities and pair with Z checks.
        check_name = {"X": "Z", "Z": "X"}.get(name, name)
        if check_name not in dg:
            return None
        projected = _project_columns(columns, dg[check_name], obs_ids)
        if any(d.bit_count() > 2 for d, _ in projected):
            return None
        try:
            witness = columns_dem(projected).shortest_graphlike_error(ignore_ungraphlike_errors=False)
        except ValueError:
            return None
        out[name] = len(witness)
    return {"lower_bound": min(out.values()), "sector_distances": out,
            "method": "exact_graphlike_observable_projections"}


def sat_certify(columns, *, lower_bound: int, max_prove: int, seconds: float) -> dict:
    """UNSAT at weight <=w certifies physical distance >=w+1.

    Cardinality variables select complete single-fault DEM columns.  XOR
    constraints force every detector to zero and at least one logical parity
    to one.  A wall timeout returns UNKNOWN, never UNSAT.  This mode is
    optional because proving distance 5--7 may require substantial resources.
    """
    from pysat.card import ITotalizer
    from pysat.solvers import Solver
    started = time.monotonic()
    n = len(columns)
    result: dict[str, Any] = {"lower_bound": lower_bound, "attempts": [], "method": "full_dem_xor_sat"}
    if seconds <= 0 or lower_bound > max_prove:
        return result
    clauses: list[list[int]] = []
    last_var = n

    def parity(lits):
        nonlocal last_var
        if not lits:
            return None
        out = lits[0]
        for lit in lits[1:]:
            last_var += 1
            z = last_var
            clauses.extend([[-out, -lit, -z], [out, lit, -z], [out, -lit, z], [-out, lit, z]])
            out = z
        return out

    detector_lits = defaultdict(list)
    observable_lits = defaultdict(list)
    for idx, (d, o) in enumerate(columns, 1):
        for j in _indices(d):
            detector_lits[j].append(idx)
        for j in _indices(o):
            observable_lits[j].append(idx)
    for lits in detector_lits.values():
        out = parity(lits)
        if out is not None:
            clauses.append([-out])
    logical_parities = [parity(lits) for lits in observable_lits.values()]
    clauses.append([v for v in logical_parities if v is not None])
    with ITotalizer(lits=list(range(1, n + 1)), ubound=max_prove, top_id=last_var) as totalizer:
        clauses.extend(totalizer.cnf.clauses)
        with Solver(name="minisat22", bootstrap_with=clauses) as solver:
            del clauses
            result["variables"] = solver.nof_vars()
            result["clauses"] = solver.nof_clauses()
            for weight in range(lower_bound, max_prove + 1):
                remaining = seconds - (time.monotonic() - started)
                if remaining <= 0:
                    result["status"] = "timeout"
                    break
                timer = threading.Timer(remaining, solver.interrupt)
                timer.daemon = True
                timer.start()
                at = time.monotonic()
                try:
                    answer = solver.solve_limited(assumptions=[-totalizer.rhs[weight]], expect_interrupt=True)
                finally:
                    timer.cancel()
                attempt = {"weight_at_most": weight, "seconds": time.monotonic() - at,
                           "status": "SAT" if answer is True else "UNSAT" if answer is False else "UNKNOWN"}
                result["attempts"].append(attempt)
                if answer is None:
                    result["status"] = "timeout"
                    break
                if answer:
                    model = set(v for v in solver.get_model() if v > 0)
                    result["witness_indices"] = [i for i in range(n) if i + 1 in model]
                    result["status"] = "witness"
                    break
                result["lower_bound"] = weight + 1
            else:
                result["status"] = "completed"
    result["seconds"] = time.monotonic() - started
    return result


def _heuristic_child(circuit_text, size, degree, output):
    try:
        c = stim.Circuit(circuit_text)
        errors = c.search_for_undetectable_logical_errors(
            dont_explore_detection_event_sets_with_size_above=size,
            dont_explore_edges_with_degree_above=degree,
            dont_explore_edges_increasing_symptom_degree=True,
            canonicalize_circuit_errors=True,
        )
        cols = []
        for error in errors:
            d = o = 0
            if not error.circuit_error_locations:
                raise ValueError("Search returned an error with no circuit explanation")
            for term in error.dem_error_terms:
                t = term.dem_target
                if t.is_relative_detector_id():
                    d ^= 1 << t.val
                elif t.is_logical_observable_id():
                    o ^= 1 << t.val
            cols.append((d, o))
        output.put({"columns": cols})
    except Exception as ex:
        output.put({"error": repr(ex)})


def heuristic_upper(circuit, *, seconds=60, size=8, degree=8):
    """Run Stim's heuristic with a hard subprocess wall-clock deadline."""
    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    queue = ctx.Queue()
    proc = ctx.Process(target=_heuristic_child, args=(str(circuit), size, degree, queue))
    proc.start()
    try:
        import queue as queue_module
        try:
            result = queue.get(timeout=seconds)
        except queue_module.Empty:
            result = {"error": "heuristic wall-clock timeout"}
    finally:
        if proc.is_alive():
            proc.terminate()
        proc.join()
        queue.close()
    if "columns" in result:
        return physical_witness(circuit, result["columns"])
    return result


def audit_circuit(circuit, metadata, *, seconds=60, sat_seconds=0, max_prove=6,
                  triple_seconds=30, pair_seconds=0, projected_five_seconds=0) -> dict:
    started = time.monotonic()
    dem = circuit.detector_error_model(decompose_errors=False, flatten_loops=True,
                                       approximate_disjoint_errors=False)
    columns = dem_columns(dem)
    result: dict[str, Any] = {
        "circuit_sha256": circuit_hash(circuit), "rounds": metadata.get("rounds"),
        "idle_scale": metadata.get("idle_scale"),
        "schedule_sha256": hashlib.sha256(json.dumps(metadata.get("schedule"), sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
        "schedule": metadata.get("schedule"),
        "noise": {"model": metadata.get("noise_model"), "p": metadata.get("p"),
                  "idle_scale": metadata.get("idle_scale")}, "boundary": metadata.get("boundary"),
        "fault_model": "complete_undecomposed_dem_single_physical_Pauli_faults",
        "num_detectors": circuit.num_detectors, "num_observables": circuit.num_observables,
        "num_distinct_fault_responses": len(columns), "lower_bound": 1,
        "upper_bound": None, "exact": False, "certificates": [],
    }

    def accept_witness(witness, method):
        if "weight_upper_bound" not in witness:
            result.setdefault("search_notes", []).append({"method": method, **witness})
            return
        weight = witness["weight_upper_bound"]
        if result["upper_bound"] is None or weight < result["upper_bound"]:
            result["upper_bound"] = weight
            result["witness"] = witness
            result["upper_bound_method"] = method

    saved_path = Path(__file__).with_name("distance_bounds.json")
    if saved_path.exists():
        saved = json.loads(saved_path.read_text()).get("codes", {}).get(metadata.get("code"), {})
        if saved.get("circuit_sha256") == result["circuit_sha256"] and "witness" in saved:
            prior_cols = dem_columns(stim.DetectorErrorModel(saved["witness"]["dem_errors"]))
            accept_witness(physical_witness(circuit, prior_cols), "saved_physical_witness_independently_revalidated")

    exact_small = small_weight_audit(columns, triple_seconds=triple_seconds)
    result["lower_bound"] = exact_small["lower_bound"]
    idx = exact_small.pop("witness_indices", None)
    result["certificates"].append(exact_small)
    if idx is not None:
        accept_witness(physical_witness(circuit, [columns[i] for i in idx]), exact_small["method"])
    projected = sector_graph_lower_bound(columns, metadata)
    if projected:
        result["lower_bound"] = max(result["lower_bound"], projected["lower_bound"])
        result["certificates"].append(projected)
    static_path = Path(__file__).with_name("static_distance_results.json")
    if static_path.exists():
        static_results = json.loads(static_path.read_text()).get(metadata.get("code"), {})
        for sector, entry in static_results.items():
            if sector in ("X", "Z") and "witness" in entry:
                try:
                    accept_witness(static_logical_physical_witness(circuit, entry["witness"], sector),
                                   "static_logical_realized_at_last_noisy_CNOT_locations")
                except ValueError as ex:
                    result.setdefault("search_notes", []).append({"static_witness_rejected": str(ex)})
    if pair_seconds > 0 and result["lower_bound"] < 5:
        four = projected_four_audit(columns, metadata, seconds=pair_seconds)
        if four:
            physical_cols = four.pop("physical_columns")
            result["lower_bound"] = max(result["lower_bound"], four["lower_bound"])
            result["certificates"].append(four)
            if physical_cols:
                accept_witness(physical_witness(circuit, physical_cols), "exact_projected_weight4_physical_witness")
            if projected_five_seconds > 0 and result["lower_bound"] < 6 and four["lower_bound"] >= 5:
                five = projected_five_audit(columns, metadata, four, seconds=projected_five_seconds)
                if five:
                    physical_cols = five.pop("physical_columns")
                    result["lower_bound"] = max(result["lower_bound"], five["lower_bound"])
                    result["certificates"].append(five)
                    if physical_cols:
                        accept_witness(physical_witness(circuit, physical_cols), "exact_projected_weight5_physical_witness")
    if result["upper_bound"] != result["lower_bound"]:
        try:
            graph = dem.shortest_graphlike_error(ignore_ungraphlike_errors=True)
            accept_witness(physical_witness(circuit, dem_columns(graph)), "graphlike_subset_physical_witness")
        except ValueError as ex:
            result.setdefault("search_notes", []).append({"graphlike_subset": str(ex)})
    if seconds > 0 and result["upper_bound"] != result["lower_bound"]:
        accept_witness(heuristic_upper(circuit, seconds=seconds), "stim_truncated_search_upper_bound_only")
    if sat_seconds > 0 and result["upper_bound"] != result["lower_bound"]:
        target = max_prove
        if result["upper_bound"] is not None:
            target = min(target, result["upper_bound"] - 1)
        sat = sat_certify(columns, lower_bound=result["lower_bound"], max_prove=target, seconds=sat_seconds)
        result["lower_bound"] = max(result["lower_bound"], sat["lower_bound"])
        idx = sat.pop("witness_indices", None)
        result["certificates"].append(sat)
        if idx is not None:
            accept_witness(physical_witness(circuit, [columns[i] for i in idx]), "full_dem_sat_physical_witness")
    if result["upper_bound"] is not None and result["lower_bound"] > result["upper_bound"]:
        raise AssertionError("Inconsistent distance bounds; do not publish this audit")
    result["exact"] = result["upper_bound"] == result["lower_bound"]
    result["method"] = "certified_lower_bound_and_explicit_physical_fault_witness"
    result["elapsed_seconds"] = time.monotonic() - started
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--codes", nargs="+", default=["w7", "bb72", "surface7"])
    parser.add_argument("--rounds", type=int, default=None, help="Default: static distance of each code")
    parser.add_argument("--idle-scale", type=float, default=0)
    parser.add_argument("--seconds", type=float, default=60, help="maximum Stim heuristic seconds per code")
    parser.add_argument("--triple-seconds", type=float, default=30)
    parser.add_argument("--sat-seconds", type=float, default=0, help="optional SAT proof budget per code")
    parser.add_argument("--pair-seconds", type=float, default=60, help="exact weight4 pair enumeration budget per CSS sector; about 1 GB peak")
    parser.add_argument("--projected-five-seconds", "--five-seconds", type=float, default=120,
                        help="exact weight5 exclusion through connected triples, seconds per CSS sector; 0 disables")
    parser.add_argument("--max-prove", type=int, default=6, help="largest weight to exclude by SAT")
    parser.add_argument("--output", default="distance_bounds.json")
    args = parser.parse_args()
    from circuits import build_circuit, get_code
    destination = Path(args.output)
    report = {"schema_version": 1, "codes": {}}
    if destination.exists():
        report = json.loads(destination.read_text())
    for name in args.codes:
        circuit, meta = build_circuit(name, p=0.001, rounds=args.rounds if args.rounds is not None else get_code(name)["static_distance"], idle_scale=args.idle_scale)
        result = audit_circuit(circuit, meta, seconds=args.seconds, sat_seconds=args.sat_seconds,
                               max_prove=args.max_prove, triple_seconds=args.triple_seconds,
                               pair_seconds=args.pair_seconds, projected_five_seconds=args.projected_five_seconds)
        report["codes"][name] = result
        temp = destination.with_suffix(destination.suffix + ".tmp")
        temp.write_text(json.dumps(report, indent=2), encoding="utf-8")
        temp.replace(destination)
        symbol = "=" if result["exact"] else " in "
        print(f"{name}: d_circuit{symbol}[{result['lower_bound']}, {result['upper_bound']}], "
              f"{result['elapsed_seconds']:.1f} s", flush=True)
    print(f"Saved {destination.resolve()}", flush=True)


if __name__ == "__main__":
    main()
