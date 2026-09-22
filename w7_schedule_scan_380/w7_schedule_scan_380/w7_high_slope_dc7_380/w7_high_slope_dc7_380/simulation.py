"""Reproducible, full-wave Monte Carlo simulation for the W7 schedule scan.

Each wave is committed before sampling. All its shots are counted before the
stopping rule is evaluated; provisional completion-order statistics never
choose which shots enter a completed wave. Sequential stopping is still an
approximate statistical procedure, not an anytime-valid confidence guarantee.
"""
from __future__ import annotations

import os
for _thread_variable in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import math
import multiprocessing as mp
from pathlib import Path
import re
import signal
import sys
import time

DECODER_SETTINGS = {
    "det_beam": 2,
    "beam_climbing": False,
    "no_revisit_dets": True,
    "merge_errors": True,
    "pqlimit": 200_000,
    "sparsify_errors": False,
    "num_det_orders": 20,
    "seed": 2384753,
}
ROUNDS = 7
IDLE_SCALE = 0.0
_WORK_CIRCUIT = None
_WORK_DECODER = None
_WORK_K = None


def _utc():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _json_default(value):
    if hasattr(value, "tolist"):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value)}")


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=_json_default)


def _hash(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, default=_json_default)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _statistics(counts):
    n, f = counts["N"], counts["Fjoint"]
    p = f / n if n else None
    sigma = math.sqrt(p * (1 - p) / n) if n else None
    lo = hi = None
    if n:
        denominator = 1 + 1 / n
        center = (p + 0.5 / n) / denominator
        half = math.sqrt(p * (1 - p) / n + 0.25 / n**2) / denominator
        lo, hi = max(0.0, center - half), min(1.0, center + half)
    return {
        "joint_probability": p,
        "joint_sigma": sigma,
        "relative_sigma": sigma / p if p else None,
        "wilson_z": 1.0,
        "joint_wilson_lower": lo,
        "joint_wilson_upper": hi,
        "sector_X_known_failures": counts["FX"],
        "sector_Z_known_failures": counts["FZ"],
        "sector_X_operational_upper_failures": counts["FX"] + counts["decodefail"],
        "sector_Z_operational_upper_failures": counts["FZ"] + counts["decodefail"],
    }


def _stop(counts, settings):
    n, f = counts["N"], counts["Fjoint"]
    precision = (n >= settings["min_shots"] and f > 0 and
                 _statistics(counts)["relative_sigma"] <= settings["relative_error"])
    if precision:
        return "precision_reached"
    if n >= settings["max_shots"]:
        return "max_shots"
    return None


def _seed(run_seed, experiment_id, batch_id):
    text = f"w7-scan-joint-wave-v1|{run_seed}|{experiment_id}|{batch_id}"
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "little")


def _worker_init(circuit_text, dem_text, settings, k):
    global _WORK_CIRCUIT, _WORK_DECODER, _WORK_K
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    import stim
    from decoder import CheckedTesseract
    _WORK_CIRCUIT = stim.Circuit(circuit_text)
    _WORK_DECODER = CheckedTesseract(stim.DetectorErrorModel(dem_text), settings)
    _WORK_K = k


def _worker_batch(batch_id, n, seed):
    started = time.monotonic()
    detectors, observables = _WORK_CIRCUIT.compile_detector_sampler(seed=seed).sample(
        n, separate_observables=True)
    counts, exceptions = _WORK_DECODER.count_batch(detectors, observables, _WORK_K)
    return {"batch_id": batch_id, "counts": counts,
            "seconds": time.monotonic() - started, "exception_examples": exceptions}


def _probe(circuit_text, dem_text, settings, k):
    _worker_init(circuit_text, dem_text, settings, k)
    # Diagnostic shots do not enter any campaign count.
    result = _worker_batch(-1, 4, 8124717)
    peak = 0
    try:
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        peak *= 1 if sys.platform == "darwin" else 1024
    except (ImportError, AttributeError):
        pass
    return {"peak_rss_bytes": peak, "seconds_for_four_shots": result["seconds"],
            "decodefail": result["counts"]["decodefail"]}


def _available_memory():
    available = None
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                available = int(line.split()[1]) * 1024
                break
    except (OSError, ValueError):
        pass
    try:
        maximum = Path("/sys/fs/cgroup/memory.max").read_text().strip()
        if maximum != "max":
            remaining = int(maximum) - int(Path("/sys/fs/cgroup/memory.current").read_text())
            available = remaining if available is None else min(available, remaining)
    except (OSError, ValueError):
        pass
    return available


def _worker_limit(requested, probe):
    try:
        cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        cpus = os.cpu_count() or 1
    try:
        quota, period = Path("/sys/fs/cgroup/cpu.max").read_text().split()
        if quota != "max":
            cpus = min(cpus, max(1, math.floor(int(quota) / int(period))))
    except (OSError, ValueError):
        pass
    available = _available_memory()
    estimate = max(128 * 1024**2, int(1.5 * probe["peak_rss_bytes"]))
    actual = max(1, min(requested, cpus))
    if available is not None:
        actual = min(actual, max(1, int(0.75 * available // estimate)))
    return actual, {"requested_workers": requested, "workers": actual,
                    "available_cpus": cpus, "available_memory_bytes": available,
                    "estimated_bytes_per_worker": estimate,
                    "note": "Memory estimate, not a hard per-shot memory bound."}


def _merge(total, delta):
    for key, value in delta.items():
        if isinstance(value, list):
            if len(total[key]) != len(value):
                raise ValueError(f"Count vector length changed: {key}")
            total[key] = [x + y for x, y in zip(total[key], value)]
        else:
            total[key] += value
    if total["Fjoint"] != total["FX"] + total["FZ"] - total["Fboth"] + total["decodefail"]:
        raise AssertionError("Same-shot union count consistency failure")


def _next_wave_size(counts, settings, workers):
    n, f = counts["N"], counts["Fjoint"]
    if n < settings["min_shots"]:
        desired = settings["min_shots"] - n
    elif f:
        estimate = f / n
        desired = max(1, math.ceil((1 - estimate) / (estimate * settings["relative_error"]**2)) - n)
    else:
        desired = max(n, settings["min_shots"])
    return max(0, min(desired, settings["max_shots"] - n,
                      workers * settings["batch_size"]))


def _plan_wave(result, settings, workers):
    ledger = result["ledger"]
    if ledger["pending"] or ledger["active_wave"] is not None:
        raise AssertionError("Cannot create a wave before finishing the previous one")
    n = _next_wave_size(result["counts"], settings, workers)
    if n <= 0:
        raise AssertionError("No next wave despite unmet stop rule")
    chunk = min(settings["batch_size"], max(1, math.ceil(n / workers)))
    batch_ids = []
    remaining = n
    while remaining:
        batch_id = ledger["next_batch_id"]
        ledger["next_batch_id"] += 1
        batch_n = min(chunk, remaining)
        ledger["pending"][str(batch_id)] = batch_n
        batch_ids.append(batch_id)
        remaining -= batch_n
    ledger["active_wave"] = {
        "wave_id": ledger["next_wave_id"], "planned_shots": n,
        "start_N": result["counts"]["N"], "batch_ids": batch_ids,
        "created_utc": _utc(),
    }
    ledger["next_wave_id"] += 1


def _progress(result, active=0, label="provisional"):
    counts = result["counts"]
    relative = _statistics(counts)["relative_sigma"]
    rel = "n/a" if relative is None else f"{100*relative:.1f}%"
    wave = result["ledger"]["active_wave"]
    wave_id = wave["wave_id"] if wave else "-"
    print(f'[{_utc()}] {result["candidate_id"]} p={result["p"]:g} '
          f'N={counts["N"]:,} Fjoint={counts["Fjoint"]:,} FX={counts["FX"]:,} '
          f'FZ={counts["FZ"]:,} Fboth={counts["Fboth"]:,} '
          f'decodefail={counts["decodefail"]:,} relative1sigma={rel} '
          f'active={active} pending_batches={len(result["ledger"]["pending"])} '
          f'wave={wave_id} {label}', flush=True)


def _terminate(pool):
    # Python 3.12 lacks ProcessPoolExecutor.terminate_workers().
    for process in list((getattr(pool, "_processes", None) or {}).values()):
        process.terminate()


def simulate_point(candidate: dict, audit: dict, p: float, output: Path,
                   workers: int = 380, min_shots: int = 1000,
                   max_shots: int = 100_000_000, relative_error: float = .20,
                   batch_size: int = 64, run_seed: int = 20260921) -> dict:
    """Simulate one independently seeded candidate at a physical error rate p.

    The caller serializes points and owns the campaign writer lock. On resume,
    already committed waves finish even if a newly lowered max_shots cap has
    been passed. Changing workers/precision/chunk sizes never changes identity.
    First Ctrl+C completes and saves the entire current wave, then returns a
    result with status ``interrupted``; the caller must stop its outer scan.
    A second Ctrl+C terminates children and raises KeyboardInterrupt.
    """
    if not math.isfinite(p) or not 0 < p < .5:
        raise ValueError("Require finite 0 < p < .5")
    if min(workers, min_shots, max_shots, batch_size) < 1:
        raise ValueError("workers/min_shots/max_shots/batch_size must be positive")
    if not 0 < relative_error < 1 or run_seed < 0:
        raise ValueError("Require 0<relative_error<1 and run_seed>=0")
    if candidate.get("code", "w7") != "w7":
        raise ValueError("This schedule scan accepts W7 candidates only")
    from circuits import build_circuit
    from decoder import empty_counts

    circuit, metadata = build_circuit("w7", p, rounds=ROUNDS,
                                     idle_scale=IDLE_SCALE, schedule=candidate["schedule"])
    k = int(metadata["k"])
    if circuit.num_observables != 2 * k:
        raise ValueError("Circuit must contain simultaneous 2*k logical observables")
    dem = circuit.detector_error_model(decompose_errors=False, flatten_loops=True,
                                       approximate_disjoint_errors=False)
    circuit_text, dem_text = str(circuit), str(dem)
    candidate_id = str(candidate["id"])
    versions = {name: version(name) for name in ("stim", "tesseract-decoder", "numpy")}
    identity = {
        "protocol": "w7_joint_same_shot_scan_full_waves_v1",
        "candidate_id": candidate_id, "code": "w7", "p": p,
        "rounds": ROUNDS, "idle_scale": IDLE_SCALE, "run_seed": run_seed,
        "schedule": metadata["schedule"], "circuit_sha256": _hash(circuit_text),
        "dem_sha256": _hash(dem_text), "versions": versions,
        "decoder": "tesseract", "decoder_settings": dict(DECODER_SETTINGS),
        "decoder_adapter_sha256": _hash(Path(__file__).with_name("decoder.py").read_text()),
        "count_rule": "Fjoint=FX+FZ-Fboth+decodefail; no postselection",
    }
    identity = json.loads(_canonical(identity))
    experiment_id = _hash(_canonical(identity))
    safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", candidate_id).strip(".") or "candidate"
    point_dir = Path(output).expanduser().resolve() / "points" / f"{safe_id}_p{p:.8g}_{experiment_id[:12]}"
    point_dir.mkdir(parents=True, exist_ok=True)
    path = point_dir / "result.json"
    settings = {"workers": workers, "min_shots": min_shots, "max_shots": max_shots,
                "relative_error": relative_error, "batch_size": batch_size}
    if path.exists():
        result = json.loads(path.read_text(encoding="utf-8"))
        if result["identity"] != identity or result["experiment_id"] != experiment_id:
            raise RuntimeError(f"Incompatible simulation identity: {path}")
    else:
        result = {
            "schema_version": 1, "experiment_id": experiment_id,
            "candidate_id": candidate_id, "code_name": "w7", "p": p,
            "rounds": ROUNDS, "k": k, "identity": identity,
            "counts": empty_counts(k), "status": "prepared", "point_done": False,
            "created_utc": _utc(), "elapsed_seconds": 0.0, "worker_seconds": 0.0,
            "ledger": {"next_batch_id": 0, "next_wave_id": 0,
                       "completed_batches": 0, "pending": {}, "active_wave": None},
            "wave_history": [], "decoder_exception_examples": [],
        }
    result.update(candidate=candidate, audit=audit, metadata=metadata,
                  stop_settings=settings, result_path=str(path))
    (point_dir / "circuit.stim").write_text(circuit_text, encoding="utf-8")
    (point_dir / "model.dem").write_text(dem_text, encoding="utf-8")
    _atomic_json(point_dir / "config.json", {"identity": identity, "candidate": candidate,
                                            "metadata": metadata, "audit": audit,
                                            "stop_settings": settings})
    elapsed_before = result["elapsed_seconds"]
    started = time.monotonic()

    def save():
        result["updated_utc"] = _utc()
        result["elapsed_seconds"] = elapsed_before + time.monotonic() - started
        result["statistics"] = _statistics(result["counts"])
        _atomic_json(path, result)

    def close_wave():
        ledger = result["ledger"]
        if ledger["pending"]:
            raise AssertionError("Cannot close an incomplete wave")
        wave = ledger["active_wave"]
        if wave is None:
            return
        if result["counts"]["N"] - wave["start_N"] != wave["planned_shots"]:
            raise AssertionError("Completed wave shot count does not match its committed plan")
        result["wave_history"].append({
            "wave_id": wave["wave_id"], "planned_shots": wave["planned_shots"],
            "N": result["counts"]["N"], "Fjoint": result["counts"]["Fjoint"],
            "relative_sigma": _statistics(result["counts"])["relative_sigma"],
            "completed_utc": _utc(),
        })
        ledger["active_wave"] = None
        result["status"] = "wave_end"
        save()
        _progress(result, label="wave_end")

    ledger = result["ledger"]
    if ledger["pending"] and ledger["active_wave"] is None:
        raise RuntimeError("Pending batches have no committed wave; refusing an unsafe resume")
    if not ledger["pending"]:
        close_wave()  # Recover a crash between the last batch and wave commit.
        reason = _stop(result["counts"], settings)
        if reason:
            result.update(status=reason, point_done=True)
            save()
            _progress(result, label=f"point_done {reason} resumed_without_sampling")
            return result

    result.update(status="preflight", point_done=False)
    save()
    print(f"Preparing {candidate_id} p={p:g}; R=7, idle=0, full correlated DEM, fixed beam=2.", flush=True)
    context = mp.get_context("spawn")
    probe_pool = ProcessPoolExecutor(max_workers=1, mp_context=context)
    try:
        future = probe_pool.submit(_probe, circuit_text, dem_text, DECODER_SETTINGS, k)
        while not future.done():
            wait([future], timeout=30)
            if not future.done():
                print(f"[{_utc()}] {candidate_id} preflight still running; diagnostic shots only.", flush=True)
        preflight = future.result()
    except BaseException:
        _terminate(probe_pool)
        probe_pool.shutdown(wait=False, cancel_futures=True)
        result["status"] = "interrupted" if sys.exc_info()[0] is KeyboardInterrupt else "error"
        save()
        raise
    else:
        probe_pool.shutdown(wait=True)
    actual_workers, runtime = _worker_limit(workers, preflight)
    result["runtime"] = {**runtime, "preflight": preflight}
    print(f"Workers {actual_workers}/{workers}; single-thread processes; "
          f"stop checked only after complete waves; memory estimate "
          f"{runtime['estimated_bytes_per_worker']/1024**2:.0f} MiB/worker.", flush=True)

    pool = ProcessPoolExecutor(max_workers=actual_workers, mp_context=context,
                               initializer=_worker_init,
                               initargs=(circuit_text, dem_text, DECODER_SETTINGS, k))
    interrupted = False
    interrupt_count = 0
    interruption_announced = False
    terminated = False
    last_progress = 0.0
    futures = {}
    to_submit = []
    previous_sigint_handler = signal.getsignal(signal.SIGINT)

    def request_interrupt(signum, frame):
        nonlocal interrupted, interrupt_count
        interrupt_count += 1
        interrupted = True
        if interrupt_count >= 2:
            raise KeyboardInterrupt

    # The first signal is a flag, not an exception at an arbitrary point in
    # the count/pending transaction. The second may interrupt anywhere; its
    # error path below preserves the last fully committed disk checkpoint.
    signal.signal(signal.SIGINT, request_interrupt)
    try:
        while True:
            try:
                if interrupted and not interruption_announced:
                    interruption_announced = True
                    result.update(status="draining_after_interrupt", point_done=False)
                    save()
                    print("Ctrl+C: finish/save the entire committed wave, then stop. "
                          "Press Ctrl+C again to terminate; uncounted batches will resume.", flush=True)
                if not ledger["pending"]:
                    close_wave()
                    if interrupted:
                        result.update(status="interrupted", point_done=False)
                        save()
                        _progress(result, label="interrupted wave_complete")
                        return result
                    reason = _stop(result["counts"], settings)
                    if reason:
                        result.update(status=reason, point_done=True)
                        save()
                        _progress(result, label=f"point_done {reason}")
                        return result
                    _plan_wave(result, settings, actual_workers)
                    result["status"] = "provisional"
                    save()  # Commit every batch size/id before any submission.
                    to_submit = sorted(int(i) for i in ledger["pending"])
                elif not futures and not to_submit:
                    # Resume the existing wave before examining current stop settings.
                    to_submit = sorted(int(i) for i in ledger["pending"])
                    result["status"] = "provisional"
                    save()

                while to_submit and len(futures) < actual_workers:
                    batch_id = to_submit.pop(0)
                    n = ledger["pending"][str(batch_id)]
                    seed = _seed(run_seed, experiment_id, batch_id)
                    futures[pool.submit(_worker_batch, batch_id, n, seed)] = batch_id
                finished, _ = wait(futures, timeout=2, return_when=FIRST_COMPLETED)
                for future in finished:
                    batch_id = futures.pop(future)
                    batch = future.result()
                    if batch["batch_id"] != batch_id:
                        raise AssertionError("Returned batch has the wrong id")
                    expected = ledger["pending"][str(batch_id)]
                    if batch["counts"]["N"] != expected:
                        raise AssertionError("Returned batch has the wrong number of shots")
                    _merge(result["counts"], batch["counts"])
                    del ledger["pending"][str(batch_id)]
                    ledger["completed_batches"] += 1
                    result["worker_seconds"] += batch["seconds"]
                    for example in batch["exception_examples"]:
                        if example not in result["decoder_exception_examples"] and len(result["decoder_exception_examples"]) < 10:
                            result["decoder_exception_examples"].append(example)
                    save()  # Counts and deletion from pending are one transaction.
                if time.monotonic() - last_progress >= 30:
                    _progress(result, len(futures), "provisional draining_wave" if interrupted else "provisional")
                    last_progress = time.monotonic()
            except KeyboardInterrupt:
                terminated = True
                _terminate(pool)
                raise
    except BaseException as exc:
        # An exception may occur between in-memory count and pending edits.
        # Change only the status of the last atomic checkpoint, never persist
        # the possibly interrupted transaction held in `result`.
        checkpoint = json.loads(path.read_text(encoding="utf-8"))
        checkpoint.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "error",
                          point_done=False, last_error=type(exc).__name__ + ": " + str(exc),
                          updated_utc=_utc())
        _atomic_json(path, checkpoint)
        if not terminated:
            terminated = True
            _terminate(pool)
        raise
    finally:
        signal.signal(signal.SIGINT, previous_sigint_handler)
        pool.shutdown(wait=not terminated, cancel_futures=terminated)
