"""Resume-safe physical circuit-distance screening for the W7 schedule scan.

Only complete R=7, p=0.003, idle=0 circuits are screened.  A heuristic supplies
UPPER bounds only.  By default a circuit qualifies only after an independently
computed lower bound reaches the requested distance.  The lower-bound proof
projects every supported physical fault into each CSS component, then exactly
excludes all logical sets of up to five faults.  Projection is a relaxation,
so its lower bounds remain sound for correlated X/Z physical faults.

``proof_seconds`` is a total lower-bound proof budget per candidate; a timeout
is inconclusive.  Completed per-candidate files are atomic checkpoints.  The
complete circuit, DEM, source files, and algorithm version identify a cached
proof.  Increasing a budget automatically retries an inconclusive result.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import re
import signal
import tempfile
import time
import traceback
from typing import Any


ALGORITHM_VERSION = "w7-screen-full-r7-connected-five-v1"
ROUNDS = 7
PHYSICAL_P = 0.003
IDLE_SCALE = 0.0
BYTES_PER_JOB = 1536 * 1024**2
_DIRECTORY = Path(__file__).resolve().parent


def _canonical_hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def _atomic_json(path: Path, value: Any):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _read_json(path: Path) -> dict | None:
    try:
        result = json.loads(path.read_text(encoding="utf-8"))
        return result if isinstance(result, dict) else None
    except (OSError, ValueError):
        return None


def _available_memory() -> int:
    limits = []
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                limits.append(int(line.split()[1]) * 1024)
                break
    except (OSError, ValueError):
        pass
    try:
        ceiling = Path("/sys/fs/cgroup/memory.max").read_text().strip()
        used = int(Path("/sys/fs/cgroup/memory.current").read_text())
        if ceiling != "max":
            limits.append(max(0, int(ceiling) - used))
    except (OSError, ValueError):
        pass
    if not limits:
        try:
            limits.append(os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE"))
        except (ValueError, OSError, AttributeError):
            # Unknown memory: permit one conservative proof worker only.
            return 2 * BYTES_PER_JOB
    return min(limits)


def _cpu_budget() -> tuple[int, list[int]]:
    try:
        affinity = sorted(os.sched_getaffinity(0))
    except AttributeError:
        affinity = list(range(os.cpu_count() or 1))
    budget = max(1, len(affinity))
    try:
        quota, period = Path("/sys/fs/cgroup/cpu.max").read_text().split()
        if quota != "max":
            budget = min(budget, max(1, math.floor(int(quota) / int(period))))
    except (OSError, ValueError):
        pass
    return budget, affinity


def _thread_limits():
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        os.environ[key] = "1"


def _worker_initialize(affinity):
    _thread_limits()
    # Put proof workers and their optional Stim-search children into separate
    # process groups, so cancelling the scan can stop the entire worker tree.
    if hasattr(os, "setsid"):
        try:
            os.setsid()
        except OSError:
            pass
    identity = mp.current_process()._identity
    if affinity and hasattr(os, "sched_setaffinity"):
        try:
            os.sched_setaffinity(0, {affinity[((identity[-1] if identity else 1) - 1) % len(affinity)]})
        except OSError:
            pass


def _source_identity() -> dict:
    import stim
    result = {"algorithm": ALGORITHM_VERSION, "stim_version": stim.__version__}
    for filename in ("screening.py", "distance_audit.py", "circuits.py"):
        result[filename] = hashlib.sha256((_DIRECTORY / filename).read_bytes()).hexdigest()
    return result


def _complete_dem_hash(dem) -> str:
    """Hash full labeled DEM distribution, retaining all hyperedge correlations.

    Identical detector/observable columns represent independent Bernoulli XORs
    in the exact DEM.  Merge their probabilities by multiplying (1-2p), in a
    sorted order.  Hexadecimal floating representations avoid tolerance-based
    false deduplication.  Detectors and logical-observable labels are retained.
    """
    groups: dict[tuple[int, int], list[float]] = {}
    for instruction in dem.flattened():
        if instruction.type != "error":
            continue
        probability = instruction.args_copy()[0]
        if probability == 0:
            continue
        detectors = observables = 0
        for target in instruction.targets_copy():
            if target.is_separator():
                raise ValueError("Screening requires an undecomposed complete DEM")
            if target.is_relative_detector_id():
                detectors ^= 1 << target.val
            elif target.is_logical_observable_id():
                observables ^= 1 << target.val
        if detectors or observables:
            groups.setdefault((detectors, observables), []).append(probability)
    rows = []
    for (detectors, observables), probabilities in sorted(groups.items()):
        complement = math.prod(1 - 2*p for p in sorted(probabilities))
        rows.append((hex(detectors), hex(observables), complement.hex()))
    return _canonical_hash({"num_detectors": dem.num_detectors,
                            "num_observables": dem.num_observables, "errors": rows})


def _classification(audit, minimum, allow_uncertified):
    lower, upper = audit.get("lower_bound", 1), audit.get("upper_bound")
    if upper is not None and upper < minimum:
        return "rejected", False
    if lower >= minimum:
        return "qualified", True
    if allow_uncertified and upper is not None and upper >= minimum:
        return "qualified", False
    return "inconclusive", False


def _identity_for(candidate, circuit, dem, metadata, sources, minimum):
    from distance_audit import circuit_hash
    identity = {
        "sources": sources, "circuit_sha256": circuit_hash(circuit),
        "dem_sha256": _complete_dem_hash(dem), "rounds": ROUNDS,
        "p": PHYSICAL_P, "idle_scale": IDLE_SCALE,
        "schedule_sha256": _canonical_hash(candidate["schedule"]),
        "boundary": metadata["boundary"], "min_distance": minimum,
    }
    return {**identity, "pipeline_sha256": _canonical_hash(identity)}


def _new_record(candidate, identity, metadata, destination, minimum, heuristic_seconds,
                proof_seconds, allow_uncertified):
    return {
        **candidate, "screening_schema_version": 1, "identity": identity,
        "metadata": metadata, "status": "inconclusive", "certified": False,
        "distance_report_path": str(destination.resolve()), "running": True,
        "phase": "pending", "min_distance": minimum,
        "allow_uncertified": allow_uncertified,
        "budgets": {"heuristic_seconds": heuristic_seconds, "proof_seconds": proof_seconds},
        "audit": {"lower_bound": 1, "upper_bound": None, "exact": False,
                  "rounds": ROUNDS, "p": PHYSICAL_P, "idle_scale": IDLE_SCALE,
                  "circuit_sha256": identity["circuit_sha256"],
                  "dem_sha256": identity["dem_sha256"], "schedule": candidate["schedule"],
                  "certificates": [], "search_notes": [],
                  "method": "physical_upper_witness_and_exact_complete_CSS_projection_lower_bound"},
    }


def _screen_one(candidate, identity, destination_text, minimum, heuristic_seconds,
                proof_seconds, allow_uncertified):
    import stim
    from circuits import build_circuit
    from distance_audit import (circuit_hash, dem_columns, _project_columns,
                                physical_witness, static_logical_physical_witness,
                                heuristic_upper, small_weight_audit,
                                exact_four_audit, exact_five_audit)
    started, cpu_started = time.monotonic(), time.process_time()
    destination = Path(destination_text)
    previous = _read_json(destination)
    record = None
    try:
        circuit, metadata = build_circuit("w7", PHYSICAL_P, rounds=ROUNDS,
                                          idle_scale=IDLE_SCALE, schedule=candidate["schedule"])
        if circuit_hash(circuit) != identity["circuit_sha256"]:
            raise ValueError("Circuit changed between preflight and worker execution")
        dem = circuit.detector_error_model(decompose_errors=False, flatten_loops=True,
                                           approximate_disjoint_errors=False)
        if _complete_dem_hash(dem) != identity["dem_sha256"]:
            raise ValueError("Complete DEM changed between preflight and worker execution")
        columns = dem_columns(dem)
        supported = set(columns)
        record = _new_record(candidate, identity, metadata, destination, minimum,
                             heuristic_seconds, proof_seconds, allow_uncertified)
        audit = record["audit"]
        audit.update(num_detectors=circuit.num_detectors, num_observables=circuit.num_observables,
                     num_distinct_fault_responses=len(columns))

        def checkpoint(phase):
            record["phase"] = phase
            audit["exact"] = audit["upper_bound"] == audit["lower_bound"]
            record["elapsed_seconds"] = time.monotonic() - started
            record["cpu_seconds"] = time.process_time() - cpu_started
            _atomic_json(destination, record)

        def accept_witness(witness, source):
            if "dem_errors" not in witness:
                audit["search_notes"].append({"source": source, **witness})
                return False
            selected = dem_columns(stim.DetectorErrorModel(witness["dem_errors"]))
            if not all(column in supported for column in selected):
                audit["search_notes"].append({"source": source, "note": "witness does not belong to this complete R7 DEM"})
                return False
            checked = physical_witness(circuit, selected)
            upper = checked["weight_upper_bound"]
            if audit["upper_bound"] is None or upper < audit["upper_bound"]:
                audit.update(upper_bound=upper, witness=checked, upper_bound_method=source)
                checkpoint("physical_upper_witness_verified")
            return True

        def below_threshold():
            return audit["upper_bound"] is not None and audit["upper_bound"] < minimum

        def finish():
            if audit["upper_bound"] is not None and audit["lower_bound"] > audit["upper_bound"]:
                raise AssertionError("Inconsistent lower and physical upper bounds")
            record["status"], record["certified"] = _classification(audit, minimum, allow_uncertified)
            record["running"] = False
            audit["lower_bound_recomputed_for_this_pipeline"] = True
            if record["status"] == "qualified" and not record["certified"]:
                record["qualification_note"] = "Explicit allow_uncertified: upper bound is not a distance guarantee"
            checkpoint("complete")
            return record

        checkpoint("upper_bound_screen")
        seeds = []
        if previous and previous.get("audit", {}).get("witness"):
            seeds.append((previous["audit"]["witness"], "previous_result_witness_revalidated_on_current_circuit"))
        baseline = _read_json(_DIRECTORY / "baseline_distance.json") or {}
        baseline_entry = baseline.get("codes", {}).get("w7", {})
        if baseline_entry.get("witness"):
            seeds.append((baseline_entry["witness"], "baseline_witness_only_revalidated_no_lower_bound_imported"))
        for witness, source in seeds:
            try:
                accept_witness(witness, source)
            except ValueError as ex:
                audit["search_notes"].append({"source": source, "rejected_seed": str(ex)})
        if below_threshold():
            return finish()

        # Static logicals are converted to actual supported last-CNOT faults.
        if audit["upper_bound"] is None or audit["upper_bound"] > 7:
            static = _read_json(_DIRECTORY / "static_distance_results.json") or {}
            for sector, item in static.get("w7", {}).items():
                if sector in ("X", "Z") and item.get("witness"):
                    accept_witness(static_logical_physical_witness(circuit, item["witness"], sector),
                                   "static_logical_realized_as_physical_CNOT_faults")
                    break

        search_started = time.monotonic()
        if heuristic_seconds > 0:
            quick = min(10.0, heuristic_seconds)
            accept_witness(heuristic_upper(circuit, seconds=quick, size=6, degree=6),
                           "full_R7_stim_search_size6_degree6_upper_only")
            if below_threshold():
                return finish()
            remaining = heuristic_seconds - (time.monotonic() - search_started)
            if remaining > .25 and (audit["upper_bound"] is None or audit["upper_bound"] > minimum):
                accept_witness(heuristic_upper(circuit, seconds=remaining, size=8, degree=8),
                               "full_R7_stim_search_size8_degree8_upper_only")
                if below_threshold():
                    return finish()

        proof_started = time.monotonic()
        sector_state = {"X": {"lower_bound": 1}, "Z": {"lower_bound": 1}}
        projected = {
            sector: _project_columns(columns, metadata["detector_groups"][check], metadata["observable_groups"][sector])
            for sector, check in (("X", "Z"), ("Z", "X"))
        }
        # The two observable groups cover every residual logical Pauli bit.
        all_obs = set(range(circuit.num_observables))
        if set(metadata["observable_groups"]["X"]) | set(metadata["observable_groups"]["Z"]) != all_obs:
            raise ValueError("CSS projections do not cover all joint logical observables")

        for stage, required_lower in (("weights1to3", 1), ("weight4", 4), ("weight5", 5)):
            if audit["lower_bound"] >= minimum:
                break
            for sector in ("X", "Z"):
                remaining = proof_seconds - (time.monotonic() - proof_started)
                if remaining <= 0:
                    audit["search_notes"].append({"phase": "lower_bound_proof", "status": "time_budget_exhausted"})
                    return finish()
                if sector_state[sector]["lower_bound"] < required_lower:
                    continue
                checkpoint(f"prove_{stage}_{sector}")
                if stage == "weights1to3":
                    result = small_weight_audit(projected[sector], triple_seconds=remaining)
                elif stage == "weight4":
                    result = exact_four_audit(projected[sector], seconds=remaining)
                else:
                    result = exact_five_audit(projected[sector], seconds=remaining)
                indices = result.pop("witness_indices", None)
                sector_state[sector] = result
                audit["certificates"].append({"sector": sector, "stage": stage, **result})
                audit["lower_bound"] = max(audit["lower_bound"], min(v["lower_bound"] for v in sector_state.values()))
                if indices is not None:
                    selected = [projected[sector][i] for i in indices]
                    if all(column in supported for column in selected):
                        accept_witness(physical_witness(circuit, selected), f"exact_{stage}_{sector}_physical_witness")
                        if below_threshold():
                            return finish()
                    else:
                        audit["search_notes"].append({"sector": sector, "stage": stage,
                                                      "note": "projected witness has no pure-sector physical lift; no full upper bound claimed"})
                checkpoint(f"proved_{stage}_{sector}")
        audit["proof_elapsed_seconds"] = time.monotonic() - proof_started
        return finish()
    except Exception as ex:
        if record is None:
            record = {**candidate, "identity": identity, "audit": {"lower_bound": 1, "upper_bound": None, "exact": False},
                      "distance_report_path": str(destination.resolve())}
        record.update(status="error", certified=False, running=False,
                      error=f"{type(ex).__name__}: {ex}", traceback=traceback.format_exc(),
                      elapsed_seconds=time.monotonic() - started)
        _atomic_json(destination, record)
        return record


def screen_candidates(manifest: dict, output: Path, workers: int = 32, min_distance: int = 6,
                      heuristic_seconds: float = 30, proof_seconds: float = 120,
                      allow_uncertified: bool = False, retry_inconclusive: bool = False) -> list[dict]:
    """Return every candidate's checkpoint record, in manifest order.

    Status is qualified/rejected/inconclusive/duplicate/error.  Qualified with
    ``certified=True`` means the recomputed lower bound reaches min_distance;
    exact distance need not be known (e.g. [6,7] qualifies).  With the explicit
    allow_uncertified override, [5,7] can qualify with certified=False.
    """
    if workers < 1 or not 1 <= min_distance <= 7:
        raise ValueError("workers must be positive; min_distance must lie in 1..7")
    if heuristic_seconds < 0 or proof_seconds < 0:
        raise ValueError("Screening budgets cannot be negative")
    _thread_limits()
    from circuits import build_circuit
    candidates = manifest.get("candidates")
    if not isinstance(candidates, list):
        raise ValueError("Manifest requires a candidates list")
    ids = [candidate.get("id") for candidate in candidates]
    if any(not isinstance(i, str) or not re.fullmatch(r"[A-Za-z0-9_-]{1,128}", i) for i in ids):
        raise ValueError("Candidate ids must contain only letters, digits, underscore, hyphen")
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate candidate ids in manifest")
    output = Path(output).resolve()
    distance_dir = output / "distance"
    distance_dir.mkdir(parents=True, exist_ok=True)
    sources = _source_identity()
    cpus, affinity = _cpu_budget()
    available = _available_memory()
    reserve = max(1024**3, available // 5)
    memory_workers = max(1, (max(0, available-reserve)) // BYTES_PER_JOB)
    effective = max(1, min(workers, cpus, memory_workers, max(1, len(candidates))))
    resource_info = {"requested_workers": workers, "effective_workers": effective,
                     "available_memory_bytes": available, "reserved_memory_bytes": reserve,
                     "estimated_bytes_per_proof_worker": BYTES_PER_JOB, "available_cpus": cpus}
    print(f"Distance screening: {effective} workers (requested {workers}); full R={ROUNDS}, p={PHYSICAL_P:g}, idle=0", flush=True)
    records: dict[str, dict] = {}
    identities, metadata_by_id, representatives, duplicates, pending = {}, {}, {}, {}, []
    for position, candidate in enumerate(candidates, 1):
        cid = candidate["id"]
        destination = distance_dir / f"{cid}.json"
        try:
            if candidate.get("code", "w7") != "w7":
                raise ValueError("This screener only accepts W7 candidates")
            circuit, metadata = build_circuit("w7", PHYSICAL_P, rounds=ROUNDS,
                                              idle_scale=IDLE_SCALE, schedule=candidate["schedule"])
            dem = circuit.detector_error_model(decompose_errors=False, flatten_loops=True,
                                               approximate_disjoint_errors=False)
            identity = _identity_for(candidate, circuit, dem, metadata, sources, min_distance)
            identities[cid], metadata_by_id[cid] = identity, metadata
            digest = identity["dem_sha256"]
            if digest in representatives:
                duplicates[cid] = representatives[digest]
                continue
            representatives[digest] = cid
            existing = _read_json(destination)
            matches = existing and existing.get("identity", {}).get("pipeline_sha256") == identity["pipeline_sha256"]
            retry = False
            if matches and not existing.get("running", False):
                status, certified = _classification(existing.get("audit", {}), min_distance, allow_uncertified)
                old_budgets = existing.get("budgets", {})
                increased = heuristic_seconds > old_budgets.get("heuristic_seconds", 0) or proof_seconds > old_budgets.get("proof_seconds", 0)
                needs_more_proof = (not certified and status != "rejected") or existing.get("status") == "error"
                retry = needs_more_proof and (retry_inconclusive or increased)
                if existing.get("status") not in ("error", "duplicate") and not retry:
                    existing.update(status=status, certified=certified, allow_uncertified=allow_uncertified,
                                    cache_reused=True, cache_source="identical_pipeline_and_complete_circuit_hash",
                                    distance_report_path=str(destination))
                    records[cid] = existing
                    _atomic_json(destination, existing)
                    continue
                if existing.get("status") == "error" and not retry:
                    records[cid] = existing
                    continue
            pending.append(candidate)
        except Exception as ex:
            record = {**candidate, "status": "error", "certified": False, "running": False,
                      "audit": {"lower_bound": 1, "upper_bound": None, "exact": False},
                      "error": f"{type(ex).__name__}: {ex}", "traceback": traceback.format_exc(),
                      "distance_report_path": str(destination)}
            records[cid] = record
            _atomic_json(destination, record)
        if position % 10 == 0 or position == len(candidates):
            print(f"DEM preflight {position}/{len(candidates)}; {len(representatives)} distinct complete models", flush=True)

    def show(record):
        audit = record["audit"]
        print(f"{record['id']}: {record['status']} [{audit.get('lower_bound')}, {audit.get('upper_bound')}] "
              f"certified={record.get('certified', False)}", flush=True)

    if pending and effective == 1:
        for candidate in pending:
            cid = candidate["id"]
            record = _screen_one(candidate, identities[cid], str(distance_dir / f"{cid}.json"),
                                 min_distance, heuristic_seconds, proof_seconds, allow_uncertified)
            records[cid] = record
            show(record)
    elif pending:
        executor = ProcessPoolExecutor(max_workers=min(effective, len(pending)), mp_context=mp.get_context("spawn"),
                                       initializer=_worker_initialize, initargs=(affinity,))
        cancelled = False
        jobs = {}
        try:
            jobs = {
                executor.submit(_screen_one, candidate, identities[candidate["id"]],
                                str(distance_dir / f"{candidate['id']}.json"), min_distance,
                                heuristic_seconds, proof_seconds, allow_uncertified): candidate
                for candidate in pending
            }
            for future in as_completed(jobs):
                candidate = jobs[future]
                cid = candidate["id"]
                try:
                    record = future.result()
                except Exception as ex:
                    record = {**candidate, "identity": identities[cid], "status": "error", "certified": False,
                              "running": False, "error": f"worker failed: {type(ex).__name__}: {ex}",
                              "audit": {"lower_bound": 1, "upper_bound": None, "exact": False},
                              "distance_report_path": str(distance_dir / f"{cid}.json")}
                    _atomic_json(Path(record["distance_report_path"]), record)
                records[cid] = record
                show(record)
        except KeyboardInterrupt:
            cancelled = True
            for future in jobs:
                future.cancel()
            processes = list((getattr(executor, "_processes", None) or {}).values())
            for process in processes:
                try:
                    if hasattr(os, "killpg") and os.getpgid(process.pid) == process.pid:
                        os.killpg(process.pid, signal.SIGTERM)
                    else:
                        process.terminate()
                except (OSError, ProcessLookupError):
                    pass
            executor.shutdown(wait=False, cancel_futures=True)
            for process in processes:
                process.join(timeout=1)
                if process.is_alive():
                    process.kill()
            print("Screening interrupted; completed checkpoints saved. Unfinished candidates will be retried.", flush=True)
            raise
        finally:
            if not cancelled:
                executor.shutdown(wait=True)

    by_id = {candidate["id"]: candidate for candidate in candidates}
    for cid, representative in duplicates.items():
        record = _new_record(by_id[cid], identities[cid], metadata_by_id[cid], distance_dir / f"{cid}.json",
                             min_distance, heuristic_seconds, proof_seconds, allow_uncertified)
        record.update(status="duplicate", certified=False, running=False, phase="complete",
                      duplicate_of=representative,
                      duplicate_report_path=records[representative]["distance_report_path"],
                      duplicate_reason="identical full labeled DEM including merged error probabilities",
                      qualification_note="Not simulated twice; no independent lower-bound proof claimed for this circuit")
        records[cid] = record
        _atomic_json(Path(record["distance_report_path"]), record)
    ordered = [records[cid] for cid in ids]
    counts = {status: sum(record["status"] == status for record in ordered)
              for status in ("qualified", "rejected", "inconclusive", "duplicate", "error")}
    _atomic_json(output / "screening_summary.json", {
        "schema_version": 1, "algorithm": ALGORITHM_VERSION,
        "manifest_sha256": _canonical_hash(manifest), "resources": resource_info,
        "rounds": ROUNDS, "p": PHYSICAL_P, "idle_scale": IDLE_SCALE,
        "min_distance": min_distance, "allow_uncertified": allow_uncertified,
        "counts": counts, "candidate_ids": ids,
        "qualified_ids": [r["id"] for r in ordered if r["status"] == "qualified"],
        "certified_qualified_ids": [r["id"] for r in ordered if r["status"] == "qualified" and r["certified"]],
    })
    print("Screening finished: " + ", ".join(f"{key}={value}" for key, value in counts.items()), flush=True)
    return ordered
