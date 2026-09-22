#!/usr/bin/env python3
"""Run the joint X/Z circuit-level campaign. See README.md for definitions."""
from __future__ import annotations

import os
for _thread_variable in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import math
import multiprocessing as mp
from pathlib import Path
import signal
import sys
import time

import campaign_config as defaults

BASE_DIR = Path(__file__).resolve().parent
SCHEMA_VERSION = 1
_WORK_CIRCUIT = None
_WORK_DECODER = None
_WORK_K = None


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def json_default(value):
    if hasattr(value, "tolist"):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Not JSON serializable: {type(value)}")


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=json_default)


def sha256_text(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as file:
        json.dump(value, file, indent=2, ensure_ascii=False, default=json_default)
        file.write("\n")
        file.flush()
        os.fsync(file.fileno())
    os.replace(temporary, path)


def wilson_interval(failures, shots, z=1.0):
    if shots == 0:
        return None, None
    p = failures / shots
    denominator = 1 + z * z / shots
    center = (p + z * z / (2 * shots)) / denominator
    half = z * math.sqrt(p * (1 - p) / shots + z * z / (4 * shots * shots)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def statistics(counts):
    n, f = counts["N"], counts["Fjoint"]
    lo, hi = wilson_interval(f, n)
    p = f / n if n else None
    # This estimated binomial standard error defines the approved stopping
    # rule. Wilson z=1 intervals, not a Wald interval, are used for the plot.
    sigma = math.sqrt(p * (1 - p) / n) if n else None
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


def stopping_reason(counts, args):
    if counts["N"] >= args.max_shots:
        return "max_shots"
    if counts["N"] < args.min_shots or counts["Fjoint"] == 0:
        return None
    relative = statistics(counts)["relative_sigma"]
    if relative is not None and relative <= args.relative_error:
        return "precision_reached"
    return None


def target_shots_for_next_wave(counts, args):
    """Avoid launching 376 * 64 shots when 1000 high-p shots will suffice."""
    n, f = counts["N"], counts["Fjoint"]
    if n < args.min_shots:
        target = args.min_shots
    elif f:
        estimate = f / n
        target = max(args.min_shots, math.ceil((1 - estimate) / (estimate * args.relative_error**2)))
    else:
        target = max(2 * args.min_shots, 2 * n)
    return min(args.max_shots, max(n + 1, target))


def _worker_initialize(circuit_text, dem_text, decoder_settings, k):
    global _WORK_CIRCUIT, _WORK_DECODER, _WORK_K
    # Only the parent handles Ctrl+C and saves the shared aggregate.
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    import stim
    from decoder import CheckedTesseract
    _WORK_CIRCUIT = stim.Circuit(circuit_text)
    _WORK_DECODER = CheckedTesseract(stim.DetectorErrorModel(dem_text), decoder_settings)
    _WORK_K = k


def _worker_batch(batch_id, num_shots, seed):
    started = time.monotonic()
    sampler = _WORK_CIRCUIT.compile_detector_sampler(seed=seed)
    detectors, observables = sampler.sample(num_shots, separate_observables=True)
    counts, exceptions = _WORK_DECODER.count_batch(detectors, observables, _WORK_K)
    return {"batch_id": batch_id, "counts": counts,
            "seconds": time.monotonic() - started,
            "exception_examples": exceptions}


def _preflight_worker(circuit_text, dem_text, decoder_settings, k):
    """One real child measures loaded-decoder memory before launching 376."""
    _worker_initialize(circuit_text, dem_text, decoder_settings, k)
    result = _worker_batch(0, 4, 8124717)
    rss = 0
    try:
        import resource
        # Linux reports KiB; macOS reports bytes.
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        rss *= 1 if sys.platform == "darwin" else 1024
    except (ImportError, AttributeError):
        pass
    return {"peak_rss_bytes": rss, "seconds_for_four_shots": result["seconds"],
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
    # Respect a container's cgroup v2 limit when applicable.
    try:
        maximum = Path("/sys/fs/cgroup/memory.max").read_text().strip()
        if maximum != "max":
            remaining = int(maximum) - int(Path("/sys/fs/cgroup/memory.current").read_text())
            available = remaining if available is None else min(remaining, available)
    except (OSError, ValueError):
        pass
    return available


def effective_workers(requested, preflight, no_memory_cap=False):
    try:
        cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        cpus = os.cpu_count() or 1
    workers = max(1, min(requested, cpus))
    available = _available_memory()
    memory_per_worker = max(128 * 1024**2, int(preflight["peak_rss_bytes"] * 1.5))
    if available is not None and not no_memory_cap:
        workers = min(workers, max(1, int(0.75 * available // memory_per_worker)))
    return workers, {
        "requested_workers": requested, "available_cpus": cpus,
        "workers": workers, "available_memory_bytes": available,
        "estimated_bytes_per_worker": memory_per_worker,
        "memory_cap_enabled": not no_memory_cap,
        "note": "Memory estimate includes 50% headroom; a harder shot can still use more memory.",
    }


def add_counts(total, delta):
    for name, value in delta.items():
        if isinstance(value, list):
            total[name] = [a + b for a, b in zip(total[name], value)]
        else:
            total[name] += value
    if total["Fjoint"] != total["FX"] + total["FZ"] - total["Fboth"] + total["decodefail"]:
        raise AssertionError("Same-shot union count consistency failure")


def point_seed(run_seed, experiment_id, batch_id):
    payload = f"joint-memory-v1|{run_seed}|{experiment_id}|{batch_id}"
    return int.from_bytes(hashlib.sha256(payload.encode()).digest()[:8], "little")


def _progress(result, workers, active, draining=False):
    c = result["counts"]
    s = statistics(c)
    rel = "n/a" if s["relative_sigma"] is None else f'{100*s["relative_sigma"]:.1f}%'
    p_joint = "n/a" if s["joint_probability"] is None else f'{s["joint_probability"]:.5g}'
    state = "draining" if draining else "running"
    print(f'[{utc_now()}] {result["code_name"]} p={result["p"]:g} '
          f'N={c["N"]:,} Fjoint={c["Fjoint"]:,} FX={c["FX"]:,} FZ={c["FZ"]:,} '
          f'Fboth={c["Fboth"]:,} decodefail={c["decodefail"]:,} '
          f'Pjoint={p_joint} relative1sigma={rel} workers={workers} active={active} {state}', flush=True)


def prepare_point(code_name, p, args, versions):
    from circuits import build_circuit
    from decoder import empty_counts
    rounds = defaults.ROUNDS_BY_CODE[code_name] if args.rounds is None else args.rounds
    circuit, metadata = build_circuit(code_name, p, rounds=rounds, idle_scale=args.idle_scale)
    k = int(metadata["k"])
    if circuit.num_observables != 2 * k:
        raise ValueError(f"{code_name}: circuit must contain exactly 2*k simultaneous observables")
    dem = circuit.detector_error_model(decompose_errors=False, flatten_loops=True)
    circuit_text, dem_text = str(circuit), str(dem)
    identity = {
        "protocol": "joint_bell_reference_full_circuit_v1",
        "code_name": code_name, "p": p, "rounds": rounds, "idle_scale": args.idle_scale,
        "circuit_sha256": sha256_text(circuit_text), "dem_sha256": sha256_text(dem_text),
        "decoder_adapter_sha256": sha256_text((BASE_DIR / "decoder.py").read_text()),
        "decoder": "tesseract", "decoder_settings": defaults.DECODER_SETTINGS,
        "versions": versions, "run_seed": args.seed,
        "observables": "first k X Pauli components, then k Z Pauli components; same physical shot",
        "count_rule": "Fjoint=FX+FZ-Fboth+decodefail; invalid recoveries have no assigned sector",
    }
    experiment_id = sha256_text(canonical_json(identity))
    point_dir = args.output / "points" / f"{code_name}_p{p:.8g}_R{rounds}_{experiment_id[:12]}"
    point_dir.mkdir(parents=True, exist_ok=True)
    (point_dir / "circuit.stim").write_text(circuit_text, encoding="utf-8")
    (point_dir / "model.dem").write_text(dem_text, encoding="utf-8")
    stop_settings = {"relative_error": args.relative_error, "min_shots": args.min_shots,
                     "max_shots": args.max_shots, "batch_size_limit": args.batch_size}
    atomic_json(point_dir / "config.json", {"identity": identity, "metadata": metadata,
                                            "stop_settings": stop_settings})
    result_path = point_dir / "result.json"
    if result_path.exists():
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if result["experiment_id"] != experiment_id or result["identity"] != json.loads(canonical_json(identity)):
            raise RuntimeError(f"Incompatible experiment identity in {result_path}")
    else:
        result = {
            "schema_version": SCHEMA_VERSION, "experiment_id": experiment_id,
            "code_name": code_name, "p": p, "rounds": rounds, "k": k,
            "metadata": metadata, "identity": identity, "counts": empty_counts(k),
            "status": "prepared", "created_utc": utc_now(), "elapsed_seconds": 0.0,
            "worker_seconds": 0.0, "ledger": {"next_batch_id": 0, "pending": {}, "completed_batches": 0},
            "decoder_exception_examples": [],
        }
    bounds_file = BASE_DIR / "distance_bounds.json"
    if bounds_file.exists():
        bounds = json.loads(bounds_file.read_text(encoding="utf-8"))
        candidate = bounds.get("codes", bounds).get(code_name)
        if candidate:
            schedule_text = canonical_json(metadata.get("schedule"))
            schedule_match = ("schedule" in candidate and
                              canonical_json(candidate["schedule"]) == schedule_text)
            if candidate.get("schedule_sha256"):
                schedule_match = candidate["schedule_sha256"] == sha256_text(schedule_text)
            mismatches = []
            if candidate.get("rounds") != rounds:
                mismatches.append("rounds")
            candidate_noise = candidate.get("noise") or {}
            candidate_idle = candidate.get("idle_scale", candidate_noise.get("idle_scale")
                                           if isinstance(candidate_noise, dict) else None)
            if candidate_idle != args.idle_scale:
                mismatches.append("idle_scale")
            if not schedule_match:
                mismatches.append("schedule")
            if candidate.get("boundary") and candidate["boundary"] != metadata.get("boundary"):
                mismatches.append("boundary")
            result["distance_bounds"] = (candidate if not mismatches else {
                "status": "not_applicable", "mismatches": mismatches,
                "lower_bound": None, "upper_bound": None, "exact": False,
            })
    result["stop_settings"] = stop_settings
    result["statistics"] = statistics(result["counts"])
    result["updated_utc"] = utc_now()
    atomic_json(result_path, result)
    return point_dir, result, circuit_text, dem_text


def run_point(point_dir, result, circuit_text, dem_text, args):
    result_path = point_dir / "result.json"
    original_elapsed = result.get("elapsed_seconds", 0.0)
    started = time.monotonic()

    def save():
        result["updated_utc"] = utc_now()
        result["elapsed_seconds"] = original_elapsed + time.monotonic() - started
        result["statistics"] = statistics(result["counts"])
        atomic_json(result_path, result)

    reason = stopping_reason(result["counts"], args)
    if reason:
        result["status"] = reason
        save()
        print(f'Skip {result["code_name"]} p={result["p"]:g}: {reason}, N={result["counts"]["N"]:,}', flush=True)
        return False

    context = mp.get_context("spawn")
    print(f'Preparing {result["code_name"]} p={result["p"]:g}: full correlated DEM, fixed beam=2.', flush=True)
    with ProcessPoolExecutor(max_workers=1, mp_context=context) as probe:
        future = probe.submit(_preflight_worker, circuit_text, dem_text, defaults.DECODER_SETTINGS, result["k"])
        while not future.done():
            wait([future], timeout=30)
            if not future.done():
                print(f'[{utc_now()}] Decoder memory/compatibility preflight is still running.', flush=True)
        preflight = future.result()
    workers, memory = effective_workers(args.workers, preflight, args.no_memory_cap)
    result["runtime"] = {"preflight": preflight, **memory}
    print(f'Workers: {workers}/{args.workers}; one CPU thread per worker; '
          f'estimated memory/worker {memory["estimated_bytes_per_worker"]/1024**2:.0f} MiB; resume enabled.', flush=True)
    result["status"] = "running"
    save()
    pool = ProcessPoolExecutor(max_workers=workers, mp_context=context,
                               initializer=_worker_initialize,
                               initargs=(circuit_text, dem_text, defaults.DECODER_SETTINGS, result["k"]))
    futures = {}
    ledger = result["ledger"]
    to_resume = sorted((int(batch_id), int(n)) for batch_id, n in ledger["pending"].items())
    # Pending ids are saved before submission, so a crash can only replay an
    # uncounted batch; already counted batches are never submitted twice.
    interrupted = False
    draining = False
    last_progress = 0.0
    try:
        while True:
            try:
                reason = stopping_reason(result["counts"], args)
                draining = interrupted or reason is not None
                if draining:
                    for future, batch_id in list(futures.items()):
                        if future.cancel():
                            del futures[future]
                    # Keep cancelled/unsubmitted batches pending for a stricter
                    # future run; they have not contributed to any count.
                else:
                    while len(futures) < workers:
                        outstanding_shots = sum(int(n) for n in ledger["pending"].values())
                        remaining = args.max_shots - result["counts"]["N"] - outstanding_shots
                        if to_resume:
                            batch_id, batch_n = to_resume.pop(0)
                        elif remaining > 0:
                            target_remaining = (target_shots_for_next_wave(result["counts"], args)
                                                - result["counts"]["N"] - outstanding_shots)
                            if target_remaining <= 0:
                                break
                            batch_id = ledger["next_batch_id"]
                            ledger["next_batch_id"] += 1
                            active_slots = max(1, workers - len(futures))
                            adaptive_size = max(1, math.ceil(target_remaining / active_slots))
                            batch_n = min(args.batch_size, adaptive_size, remaining)
                            ledger["pending"][str(batch_id)] = batch_n
                            save()
                        else:
                            break
                        seed = point_seed(args.seed, result["experiment_id"], batch_id)
                        futures[pool.submit(_worker_batch, batch_id, batch_n, seed)] = batch_id
                if not futures:
                    break
                finished, _ = wait(futures, timeout=2.0, return_when=FIRST_COMPLETED)
                for future in finished:
                    batch_id = futures.pop(future)
                    batch = future.result()
                    if batch["batch_id"] != batch_id:
                        raise AssertionError("Wrong returned batch id")
                    expected_n = ledger["pending"].pop(str(batch_id))
                    if batch["counts"]["N"] != expected_n:
                        raise AssertionError("Wrong returned batch shot count")
                    add_counts(result["counts"], batch["counts"])
                    ledger["completed_batches"] += 1
                    result["worker_seconds"] += batch["seconds"]
                    for example in batch["exception_examples"]:
                        if example not in result["decoder_exception_examples"] and len(result["decoder_exception_examples"]) < 10:
                            result["decoder_exception_examples"].append(example)
                    save()
                if time.monotonic() - last_progress >= 30:
                    _progress(result, workers, len(futures), draining)
                    last_progress = time.monotonic()
            except KeyboardInterrupt:
                if interrupted:
                    for process in list(getattr(pool, "_processes", {}).values()):
                        process.terminate()
                    result["status"] = "interrupted"
                    save()
                    raise
                interrupted = True
                result["status"] = "draining_after_interrupt"
                save()
                print("Ctrl+C: stop submitting; finish/save running batches. Press Ctrl+C again to stop workers immediately.", flush=True)
        result["status"] = "interrupted" if interrupted else (stopping_reason(result["counts"], args) or "paused")
        save()
        _progress(result, workers, 0, draining)
    except BaseException as exc:
        if result["status"] != "interrupted":
            result["status"] = "error"
            result["last_error"] = type(exc).__name__ + ": " + str(exc)
        save()
        for future in futures:
            future.cancel()
        raise
    finally:
        pool.shutdown(wait=not interrupted, cancel_futures=True)
    return interrupted


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=defaults.WORKERS)
    parser.add_argument("--codes", default=",".join(defaults.CODE_NAMES), help="Comma-separated w7,bb72,surface7")
    parser.add_argument("--p", default=",".join(str(p) for p in defaults.P_VALUES), help="Comma-separated physical error rates")
    parser.add_argument("--rounds", type=int, default=defaults.ROUNDS, help="Override rounds for all codes; default uses each static distance d")
    parser.add_argument("--idle-scale", type=float, default=defaults.IDLE_SCALE)
    parser.add_argument("--relative-error", type=float, default=defaults.RELATIVE_ERROR)
    parser.add_argument("--min-shots", type=int, default=defaults.MIN_SHOTS)
    parser.add_argument("--max-shots", type=int, default=defaults.MAX_SHOTS)
    parser.add_argument("--batch-size", type=int, default=defaults.BATCH_SIZE)
    parser.add_argument("--seed", type=int, default=defaults.RUN_SEED)
    parser.add_argument("--output", type=Path, default=Path(defaults.OUTPUT))
    parser.add_argument("--prepare-only", action="store_true", help="Build circuits/DEMs and identities without Monte Carlo sampling")
    parser.add_argument("--smoke", action="store_true", help="Run all selected codes at p=.003, 16 shots each, at most 2 workers")
    parser.add_argument("--no-memory-cap", action="store_true", help="Disable automatic memory-based worker cap; CPU affinity cap remains")
    args = parser.parse_args()
    args.codes = list(dict.fromkeys(x.strip() for x in args.codes.split(",") if x.strip()))
    if not args.codes or any(code not in defaults.CODE_NAMES for code in args.codes):
        parser.error("--codes must contain w7, bb72, and/or surface7")
    try:
        args.p_values = sorted({float(x) for x in args.p.split(",")}, reverse=True)
    except ValueError:
        parser.error("--p must be comma-separated numbers")
    if any(not 0 < p < 0.5 or not math.isfinite(p) for p in args.p_values):
        parser.error("Each p must satisfy 0 < p < 0.5")
    if min(args.workers, args.min_shots, args.max_shots, args.batch_size) < 1:
        parser.error("workers, rounds, min/max shots and batch size must be positive")
    if args.rounds is not None and args.rounds < 1:
        parser.error("rounds must be positive")
    if not 0 < args.relative_error < 1 or args.idle_scale < 0 or args.seed < 0:
        parser.error("Require 0<relative-error<1, idle-scale>=0 and seed>=0")
    if args.smoke:
        args.p_values = [0.003]
        args.workers = min(args.workers, 2)
        args.min_shots = args.max_shots = 16
        args.batch_size = 4
        if args.output == Path(defaults.OUTPUT):
            args.output = Path("smoke_results")
    args.output = args.output.expanduser().resolve()
    return args


def main():
    args = parse_args()
    versions = {name: version(name) for name in ("stim", "tesseract-decoder", "numpy")}
    print("Joint X/Z same-shot simulation; full correlated DEM; fixed Tesseract beam=2.", flush=True)
    print("Versions: " + canonical_json(versions), flush=True)
    print(f"R={args.rounds if args.rounds is not None else defaults.ROUNDS_BY_CODE}; idle_scale={args.idle_scale:g}; p high -> low; codes={args.codes}", flush=True)
    print(f"Stop: at least {args.min_shots} shots and relative 1sigma <= {100*args.relative_error:g}%; "
          f"maximum {args.max_shots:,} shots; no wall-time cap.", flush=True)
    args.output.mkdir(parents=True, exist_ok=True)
    # Advisory lock is held by this live file object for the whole main call.
    # Plotting is read-only and does not acquire it.
    lock_handle = (args.output / "RUNNING.lock").open("a+", encoding="utf-8")
    if sys.platform != "win32":
        import fcntl
        try:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f"Another simulation is writing {args.output}; stop it or choose another --output")
        lock_handle.seek(0)
        lock_handle.truncate()
        lock_handle.write(f"pid={os.getpid()} started={utc_now()}\n")
        lock_handle.flush()
    index = {"schema_version": SCHEMA_VERSION, "points": [], "updated_utc": utc_now()}
    # Physical error rate is the outer loop, so all codes receive their larger
    # p points before the campaign moves to rarer events.
    for p in args.p_values:
        for code_name in args.codes:
            point_dir, result, circuit_text, dem_text = prepare_point(code_name, p, args, versions)
            index["points"].append(str((point_dir / "result.json").relative_to(args.output)))
            index["updated_utc"] = utc_now()
            atomic_json(args.output / "results_index.json", index)
            if args.prepare_only:
                print(f"Prepared {code_name} p={p:g}: {point_dir}", flush=True)
            elif run_point(point_dir, result, circuit_text, dem_text, args):
                return 130
    print(f"Finished. Results: {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    mp.freeze_support()
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted; completed batches are saved and will resume.", file=sys.stderr)
        raise SystemExit(130)
