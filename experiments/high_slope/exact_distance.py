"""Rigorous, resumable W7 circuit-distance decision at weight <= 6.

A completed UNSAT answer in every CSS/logical-bit/minimum-index shard excludes
ALL undetectable logical histories with at most six faulty locations.  A SAT
answer is an upper bound only after its response is lifted to actual faults in
the complete undecomposed circuit DEM.  Timeouts are UNKNOWN, never UNSAT.

This is a bounded exact search: its mathematics is exact, but finishing an
UNSAT proof can be expensive.  Increase seconds_per_task to retry UNKNOWNs.
"""
from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed, CancelledError
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import threading
import time
from typing import Any

from screening import (_atomic_json, _canonical_hash, _read_json,
                       _available_memory, _cpu_budget, _thread_limits,
                       _worker_initialize, _source_identity, _identity_for,
                       _complete_dem_hash, screen_candidates)

ALGORITHM = "w7-full-css-weight-at-most-six-sat-minimum-index-v1"
BYTES_PER_WORKER = 256 * 1024**2
ROOT = Path(__file__).resolve().parent


def _bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def _source_hash():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def _valid_selected(columns, indices, observable, start, stop, max_weight):
    """Check a solver/cached model independently of the CNF encoding."""
    if not indices or len(indices) != len(set(indices)) or len(indices) > max_weight:
        return False
    if min(indices) < start or min(indices) >= stop:
        return False
    d = o = 0
    for i in indices:
        if not isinstance(i, int) or not 0 <= i < len(columns):
            return False
        di, oi = columns[i]
        d ^= di
        o ^= oi
    return d == 0 and bool(o & (1 << observable))


def solve_projection_task(columns, *, observable, start, stop, max_weight=6,
                          seconds=300, solver_name="minisat22"):
    """Exact decision for one partition; public for small independent checks.

    ``start <= min(selected indices) < stop`` partitions all nonempty fault
    sets.  Each supported projected response is a Boolean variable.  Repeated
    responses cancel in GF(2), so considering each distinct response once is
    an exhaustive relaxation of physical histories of at most max_weight.
    """
    begin = time.monotonic()
    if not 0 <= start < stop <= len(columns):
        raise ValueError("Invalid first-selected-index interval")
    if seconds <= 0:
        return {"status": "UNKNOWN", "reason": "zero_time_budget", "seconds": 0.0}
    if max_weight < 1:
        raise ValueError("max_weight must be positive")
    from pysat.card import ITotalizer
    from pysat.solvers import Solver
    n = len(columns)
    clauses = []
    last = n

    def parity(lits):
        nonlocal last
        if not lits:
            return None
        out = lits[0]
        for lit in lits[1:]:
            last += 1
            z = last
            clauses.extend([[-out, -lit, -z], [out, lit, -z],
                            [out, -lit, z], [-out, lit, z]])
            out = z
        return out

    detectors = defaultdict(list)
    observable_lits = []
    for i, (d, o) in enumerate(columns, 1):
        for bit in _bits(d):
            detectors[bit].append(i)
        if o & (1 << observable):
            observable_lits.append(i)
    for lits in detectors.values():
        out = parity(lits)
        clauses.append([-out])
    out = parity(observable_lits)
    clauses.append([out] if out is not None else [])
    # All responses before start are absent; at least one in this interval
    # must occur.  Later responses remain unrestricted apart from weight.
    clauses.extend([[-i] for i in range(1, start + 1)])
    clauses.append(list(range(start + 1, stop + 1)))
    totalizer = None
    assumptions = []
    try:
        if n > max_weight:
            totalizer = ITotalizer(lits=list(range(1, n + 1)),
                                   ubound=max_weight, top_id=last)
            clauses.extend(totalizer.cnf.clauses)
            assumptions = [-totalizer.rhs[max_weight]]
        remaining = seconds - (time.monotonic() - begin)
        if remaining <= 0:
            return {"status": "UNKNOWN", "reason": "encoding_time_budget",
                    "seconds": time.monotonic() - begin}
        with Solver(name=solver_name, bootstrap_with=clauses) as solver:
            del clauses
            remaining = seconds - (time.monotonic() - begin)
            if remaining <= 0:
                return {"status": "UNKNOWN", "reason": "solver_setup_time_budget",
                        "seconds": time.monotonic() - begin}
            timer = threading.Timer(remaining, solver.interrupt)
            timer.daemon = True
            timer.start()
            try:
                answer = solver.solve_limited(assumptions=assumptions,
                                              expect_interrupt=True)
            finally:
                timer.cancel()
                timer.join()
            result = {"status": "SAT" if answer is True else
                      "UNSAT" if answer is False else "UNKNOWN",
                      "seconds": time.monotonic() - begin,
                      "variables": solver.nof_vars(), "clauses": solver.nof_clauses(),
                      "solver_statistics": solver.accum_stats()}
            if answer is True:
                model = {v for v in solver.get_model() if v > 0}
                indices = [i for i in range(n) if i + 1 in model]
                if not _valid_selected(columns, indices, observable, start, stop, max_weight):
                    raise RuntimeError("SAT model failed independent parity/cardinality verification")
                result["witness_indices"] = indices
                result["weight"] = len(indices)
            if answer is None:
                result["reason"] = "solver_interrupted_at_time_budget"
            return result
    finally:
        if totalizer is not None:
            totalizer.delete()


def _task_worker(model_path, task, seconds, destination):
    began = time.monotonic()
    record = {"task": task, "seconds_budget": seconds}
    try:
        model = _read_json(Path(model_path))
        if not model:
            raise ValueError("Missing prepared model")
        payload = model["payload"]
        if _canonical_hash(payload) != task["model_sha256"]:
            raise ValueError("Prepared model hash mismatch")
        if payload["engine_sha256"] != _source_hash():
            raise ValueError("Engine changed after task preparation")
        columns = [tuple(pair) for pair in payload["projected_columns"][task["sector"]]]
        result = solve_projection_task(columns, observable=task["observable"],
                                       start=task["start"], stop=task["stop"],
                                       max_weight=task["max_weight"], seconds=seconds,
                                       solver_name=task["solver"])
        record.update(result)
    except Exception as ex:
        record.update(status="ERROR", error=f"{type(ex).__name__}: {ex}")
    record["elapsed_seconds"] = time.monotonic() - began
    _atomic_json(Path(destination), record)
    return record


def _verified_prior_lower(record, expected_identity):
    """Import only current pipeline results, not bounds typed into old logs."""
    if not record or record.get("screening_schema_version") != 1:
        return 1
    if record.get("identity") != expected_identity or record.get("running"):
        return 1
    audit = record.get("audit", {})
    if not audit.get("lower_bound_recomputed_for_this_pipeline"):
        return 1
    certificates = audit.get("certificates", [])
    sector_lower = {}
    for sector in ("X", "Z"):
        stages = {c.get("stage"): c for c in certificates if c.get("sector") == sector}
        a, b, c = (stages.get(stage, {}) for stage in ("weights1to3", "weight4", "weight5"))
        value = 1
        if a.get("lower_bound", 1) >= 4 and a.get("triple_complete"):
            value = 4
            if b.get("lower_bound", 1) >= 5 and b.get("complete"):
                value = 5
                if c.get("lower_bound", 1) >= 6 and c.get("complete"):
                    value = 6
        sector_lower[sector] = value
    return min(min(sector_lower.values()), audit.get("lower_bound", 1))


def run_distance_campaign(candidates, output, *, workers=376, seconds_per_task=300,
                          shards=16, screening_records=None, screening_workers=8,
                          proof_seconds=300, heuristic_seconds=10,
                          retry_unknown=False, solver_name="minisat22"):
    """Return candidate records with audit.lower_bound/upper_bound/exact.

    By default independently recompute the <=5 exclusion first.  A caller may
    supply records from screening.screen_candidates for these same source,
    schedule, circuit and complete DEM hashes.  Results and individual tasks
    are atomic JSON checkpoints below ``output/exact_distance``.
    """
    if workers < 1 or shards < 1 or seconds_per_task < 0:
        raise ValueError("workers/shards must be positive and time budget nonnegative")
    if solver_name not in ("minisat22", "glucose4", "glucose3"):
        raise ValueError("Use an interruptible solver: minisat22, glucose4 or glucose3")
    _thread_limits()
    import stim
    import pysat
    from circuits import build_circuit
    from distance_audit import (dem_columns, _project_columns, physical_witness,
                                static_logical_physical_witness)
    candidates = list(candidates)
    root = Path(output).resolve()
    folder = root / "exact_distance"
    folder.mkdir(parents=True, exist_ok=True)
    if screening_records is None:
        screening_records = screen_candidates({"candidates": candidates}, root / "preproof",
            workers=screening_workers, min_distance=6, heuristic_seconds=heuristic_seconds,
            proof_seconds=proof_seconds)
    prior_by_id = {r["id"]: r for r in screening_records}
    sources = _source_identity()
    states = {}
    pending = []
    static = _read_json(ROOT / "static_distance_results.json") or {}
    known_witnesses = _read_json(ROOT / "data" / "known_physical_witnesses.json") or {}

    for candidate in candidates:
        cid = candidate["id"]
        # The preparatory screener validates ids; keep this API safe when
        # caller supplies its own records too.
        import re
        if not isinstance(cid, str) or not re.fullmatch(r"[A-Za-z0-9_-]{1,128}", cid):
            raise ValueError("Unsafe candidate id")
        if cid in states:
            raise ValueError("Duplicate candidate id")
        circuit, meta = build_circuit("w7", .003, rounds=7, idle_scale=0,
                                      schedule=candidate["schedule"])
        dem = circuit.detector_error_model(decompose_errors=False, flatten_loops=True,
                                           approximate_disjoint_errors=False)
        expected = _identity_for(candidate, circuit, dem, meta, sources, 6)
        full_columns = dem_columns(dem)
        supported = set(full_columns)
        projected = {sector: _project_columns(full_columns, meta["detector_groups"][check],
                                             meta["observable_groups"][sector])
                     for sector, check in (("X", "Z"), ("Z", "X"))}
        all_obs = set(range(circuit.num_observables))
        if set(meta["observable_groups"]["X"]) | set(meta["observable_groups"]["Z"]) != all_obs:
            raise ValueError("CSS projections do not cover every logical observable")
        prior = prior_by_id.get(cid)
        lower = _verified_prior_lower(prior, expected)
        audit = {"lower_bound": lower, "upper_bound": None, "exact": False,
                 "method": ALGORITHM, "rounds": 7, "idle_scale": 0,
                 "circuit_sha256": expected["circuit_sha256"],
                 "dem_sha256": expected["dem_sha256"],
                 "lower_bound_source": "current_hash_validated_complete_through5_screening" if lower >= 6 else
                                       "no_complete_current_through5_certificate",
                 "search_notes": []}

        def accept(witness, source):
            response = dem_columns(stim.DetectorErrorModel(witness["dem_errors"]))
            if not all(pair in supported for pair in response):
                raise ValueError("Witness response not supported by current complete DEM")
            checked = physical_witness(circuit, response)
            if audit["upper_bound"] is None or checked["weight_upper_bound"] < audit["upper_bound"]:
                audit.update(upper_bound=checked["weight_upper_bound"], witness=checked,
                             upper_bound_method=source)

        if prior and prior.get("audit", {}).get("witness"):
            try:
                accept(prior["audit"]["witness"], "preproof_physical_witness_independently_revalidated")
            except ValueError as ex:
                audit["search_notes"].append({"rejected_prior_witness": str(ex)})
        bundled = known_witnesses.get(cid)
        if bundled:
            if bundled.get("circuit_sha256") != expected["circuit_sha256"]:
                audit["search_notes"].append({"bundled_witness_skipped": "circuit_sha256_mismatch"})
            else:
                try:
                    accept(bundled["witness"], "bundled_physical_witness_independently_revalidated")
                except (KeyError, TypeError, ValueError) as ex:
                    audit["search_notes"].append({"rejected_bundled_witness": str(ex)})
        if audit["upper_bound"] is None or audit["upper_bound"] > 7:
            for sector, item in static.get("w7", {}).items():
                if sector in ("X", "Z") and item.get("witness"):
                    accept(static_logical_physical_witness(circuit, item["witness"], sector),
                           "static_logical_lifted_to_current_physical_CNOT_faults")
                    break
        payload = {"algorithm": ALGORITHM, "engine_sha256": _source_hash(),
                   "identity": expected, "pysat_version": pysat.__version__,
                   "projected_columns": projected,
                   "observable_groups": meta["observable_groups"]}
        digest = _canonical_hash(payload)
        candidate_dir = folder / cid
        model_path = candidate_dir / "model.json"
        _atomic_json(model_path, {"model_sha256": digest, "payload": payload})
        state = {"candidate": candidate, "metadata": meta, "audit": audit,
                 "model_sha256": digest, "model_path": str(model_path),
                 "projected_columns": projected, "circuit": circuit, "supported": supported,
                 "tasks": {}, "results": {}, "directory": candidate_dir,
                 "pure_sector_projection_has_supported_full_lift":
                     {s: all(pair in supported for pair in cols) for s, cols in projected.items()}}
        states[cid] = state
        if audit["upper_bound"] is not None and audit["upper_bound"] <= 6:
            continue
        for sector in ("X", "Z"):
            n = len(projected[sector])
            count = min(shards, n)
            for observable in meta["observable_groups"][sector]:
                for shard in range(count):
                    start, stop = n * shard // count, n * (shard + 1) // count
                    task = {"model_sha256": digest, "sector": sector,
                            "observable": observable, "start": start, "stop": stop,
                            "max_weight": 6, "solver": solver_name}
                    task_key = f"{sector}_L{observable:02d}_part{shard:03d}"
                    task["task_sha256"] = _canonical_hash(task)
                    state["tasks"][task_key] = task
                    destination = candidate_dir / "tasks" / f"{task_key}.json"
                    cached = _read_json(destination)
                    matches = cached and cached.get("task") == task
                    if matches and cached.get("status") == "SAT":
                        matches = _valid_selected(projected[sector], cached.get("witness_indices", []),
                                                  observable, start, stop, 6)
                    reusable = matches and (cached.get("status") in ("SAT", "UNSAT") or
                        (cached.get("status") == "UNKNOWN" and not retry_unknown and
                         seconds_per_task <= cached.get("seconds_budget", -1)))
                    if reusable:
                        state["results"][task_key] = cached
                    else:
                        pending.append((cid, task_key, task, str(destination)))

    def apply_sat(cid, key, result):
        state = states[cid]
        if result.get("status") != "SAT":
            return
        task = state["tasks"][key]
        cols = state["projected_columns"][task["sector"]]
        ids = result.get("witness_indices", [])
        if not _valid_selected(cols, ids, task["observable"], task["start"], task["stop"], 6):
            raise RuntimeError("Worker SAT model failed parent verification")
        selected = [cols[i] for i in ids]
        audit = state["audit"]
        if not all(pair in state["supported"] for pair in selected):
            audit["search_notes"].append({"task": key, "projected_SAT_without_pure_physical_lift": True})
            return
        witness = physical_witness(state["circuit"], selected)
        if audit["upper_bound"] is None or witness["weight_upper_bound"] < audit["upper_bound"]:
            audit.update(upper_bound=witness["weight_upper_bound"], witness=witness,
                         upper_bound_method="exact_projected_SAT_revalidated_full_physical_witness")

    def snapshot(cid):
        state = states[cid]
        audit = state["audit"]
        values = list(state["results"].values())
        counts = {name: sum(v.get("status") == name for v in values)
                  for name in ("SAT", "UNSAT", "UNKNOWN", "ERROR")}
        complete_unsat = bool(state["tasks"]) and counts["UNSAT"] == len(state["tasks"])
        if complete_unsat:
            audit["lower_bound"] = max(audit["lower_bound"], 7)
            audit["lower_bound_source"] = "all_covering_CSS_logical_bit_index_shards_UNSAT_at_weight_le6"
        lower, upper = audit["lower_bound"], audit["upper_bound"]
        if upper is not None and lower > upper:
            raise RuntimeError("Inconsistent exact distance bounds: stop and investigate")
        audit["exact"] = lower == upper
        status = f"exact_{lower}" if audit["exact"] else "inconclusive"
        record = {**state["candidate"], "metadata": state["metadata"], "audit": audit,
                  "status": status, "certified": lower >= 6,
                  "model_sha256": state["model_sha256"],
                  "pure_sector_projection_has_supported_full_lift":
                      state["pure_sector_projection_has_supported_full_lift"],
                  "tasks_summary": {"expected": len(state["tasks"]), "completed": len(values),
                                    **counts, "all_unsat": complete_unsat},
                  "distance_report_path": str(state["directory"] / "result.json")}
        _atomic_json(state["directory"] / "result.json", record)
        return record

    for cid, state in states.items():
        for key, result in state["results"].items():
            apply_sat(cid, key, result)
        snapshot(cid)
    # Do not spend any more SAT effort after a physical <=6 witness was found.
    pending = [job for job in pending if states[job[0]]["audit"]["upper_bound"] is None or
               states[job[0]]["audit"]["upper_bound"] > 6]
    cpus, affinity = _cpu_budget()
    memory = _available_memory()
    reserve = max(1024**3, memory // 5)
    memory_workers = max(1, (memory - reserve) // BYTES_PER_WORKER)
    effective = max(1, min(workers, cpus, memory_workers, max(1, len(pending))))
    resources = {"requested_workers": workers, "effective_workers": effective,
                 "bytes_per_worker_estimate": BYTES_PER_WORKER,
                 "available_memory_bytes": memory, "reserved_memory_bytes": reserve,
                 "seconds_per_task": seconds_per_task, "shards_per_logical_bit": shards,
                 "pending_tasks": len(pending)}
    _atomic_json(folder / "resources.json", resources)
    print(f"Exact weight<=6 search: {len(pending)} tasks; {effective}/{workers} workers; "
          f"{seconds_per_task:g}s/task; UNKNOWN is inconclusive", flush=True)
    if pending:
        executor = ProcessPoolExecutor(max_workers=effective, mp_context=mp.get_context("spawn"),
                                       initializer=_worker_initialize, initargs=(affinity,))
        futures = {}
        interrupted = False
        try:
            for cid, key, task, destination in pending:
                future = executor.submit(_task_worker, states[cid]["model_path"], task,
                                         seconds_per_task, destination)
                futures[future] = (cid, key)
            completed = 0
            for future in as_completed(futures):
                cid, key = futures[future]
                try:
                    result = future.result()
                except CancelledError:
                    continue
                except Exception as ex:
                    result = {"task": states[cid]["tasks"][key], "status": "ERROR",
                              "error": f"Worker failure: {type(ex).__name__}: {ex}"}
                states[cid]["results"][key] = result
                apply_sat(cid, key, result)
                record = snapshot(cid)
                completed += 1
                if result.get("status") == "SAT" or completed % max(1, effective // 4) == 0:
                    a = record["audit"]
                    print(f"{cid}: distance [{a['lower_bound']},{a['upper_bound']}]; "
                          f"tasks {record['tasks_summary']}", flush=True)
                upper = record["audit"]["upper_bound"]
                if upper is not None and upper <= 6:
                    for other, (other_id, _) in futures.items():
                        if other_id == cid and not other.done():
                            other.cancel()
        except BaseException:
            interrupted = True
            for future in futures:
                future.cancel()
            # Do not block Ctrl-C for the full SAT budget.  Every finished
            # shard was already written atomically by its worker.
            for process in list(getattr(executor, "_processes", {}).values()):
                if process.is_alive():
                    process.terminate()
            raise
        finally:
            executor.shutdown(wait=not interrupted, cancel_futures=True)
    records = [snapshot(candidate["id"]) for candidate in candidates]
    _atomic_json(folder / "summary.json", {"algorithm": ALGORITHM, "records": records,
                                           "resources": resources})
    return records
