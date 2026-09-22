#!/usr/bin/env python3
"""Design, distance-screen and simulate alternative W7 extraction schedules."""
from __future__ import annotations

import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[_name] = "1"

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import signal
import sys


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    temp.replace(path)


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("all", "generate", "screen", "simulate", "plot"), default="all")
    parser.add_argument("--output", type=Path, default=Path("scan_results"))
    parser.add_argument("--workers", type=int, default=380, help="Maximum simulation workers; affinity and RAM caps apply")
    parser.add_argument("--distance-workers", type=int, default=32, help="Memory-intensive independent distance audits")
    parser.add_argument("--candidates", type=int, default=64, help="Number of unique generated candidate error models, including baseline")
    parser.add_argument("--max-circuits", type=int, default=12, help="Maximum qualified circuits to simulate, including baseline")
    parser.add_argument("--seed", type=int, default=20260921)
    parser.add_argument("--sample-seed", type=int, default=20260921, help="Independent Monte Carlo seed; change it in a new output for confirmation runs")
    parser.add_argument("--p", default="0.003,0.0025", help="Only 0.003 and 0.0025, or a subset, are supported")
    parser.add_argument("--distance-seconds", type=float, default=30, help="Heuristic upper-bound search time per candidate")
    parser.add_argument("--proof-seconds", type=float, default=300, help="Total lower-bound proof budget per candidate")
    parser.add_argument("--retry-inconclusive", action="store_true")
    parser.add_argument("--allow-uncertified", action="store_true", help="Exploratory only: permit UB>=6 without proving LB>=6; use a separate output directory")
    parser.add_argument("--relative-error", type=float, default=0.20)
    parser.add_argument("--min-shots", type=int, default=1000)
    parser.add_argument("--max-shots", type=int, default=100_000_000)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--only-ids", default="", help="Comma-separated candidate IDs; baseline is always included")
    args = parser.parse_args()
    try:
        args.p_values = sorted(set(float(value) for value in args.p.split(",")), reverse=True)
    except ValueError:
        parser.error("--p must be 0.003,0.0025, or a subset")
    if not args.p_values or not set(args.p_values) <= {0.003, 0.0025}:
        parser.error("This matched-circuit campaign supports only p=0.003 and p=0.0025")
    if min(args.workers, args.distance_workers, args.candidates, args.max_circuits, args.min_shots, args.max_shots, args.batch_size) < 1:
        parser.error("Worker, candidate and shot counts must be positive")
    if args.min_shots > args.max_shots or not 0 < args.relative_error < 1:
        parser.error("Require min-shots <= max-shots and 0 < relative-error < 1")
    if args.distance_seconds < 0 or args.proof_seconds <= 0:
        parser.error("Distance search budget must be nonnegative and proof budget positive")
    if args.seed < 0 or args.sample_seed < 0:
        parser.error("Seeds must be nonnegative")
    args.output = args.output.expanduser().resolve()
    return args


def ensure_manifest(args):
    path = args.output / "candidates.json"
    request = {"seed": args.seed, "sample_seed": args.sample_seed, "allow_uncertified": args.allow_uncertified}
    settings_path = args.output / "campaign.json"
    if settings_path.exists():
        saved = json.loads(settings_path.read_text())
        if saved["fixed_request"] != request:
            raise ValueError("Seed or certification mode differs from this output directory; use a separate --output.")
    else:
        atomic_json(settings_path, {"fixed_request": request, "code": "w7", "rounds": 7,
                                    "idle_scale": 0, "beam": 2, "p_supported": [0.003, 0.0025]})
    if path.exists():
        manifest = json.loads(path.read_text())
        if args.stage not in ("all", "generate") or len(manifest["candidates"]) >= args.candidates:
            return manifest
    elif args.stage not in ("all", "generate"):
        raise ValueError("No candidates.json found. Run --stage generate or --stage all first.")
    else:
        manifest = None
    from candidates import generate_candidates
    print(f"Generating {args.candidates} distinct W7 schedules; checking extraction order and error-model uniqueness.", flush=True)
    generated = generate_candidates(count=args.candidates, seed=args.seed)
    if manifest is not None:
        mapping = {item["id"]: item["schedule"] for item in generated["candidates"]}
        if any(mapping.get(item["id"]) != item["schedule"] for item in manifest["candidates"]):
            raise ValueError("Candidate generation changed existing schedules; use a new output directory.")
    atomic_json(path, generated)
    print(f"Saved {len(generated['candidates'])} candidates to {path}", flush=True)
    return generated


def choose_candidates(records, args):
    qualified = {
        record["id"]: record for record in records
        if record["status"] == "qualified"
        and record.get("audit", {}).get("upper_bound") is not None
        and record["audit"]["upper_bound"] >= 6
        and (args.allow_uncertified or record["audit"].get("lower_bound", 1) >= 6)
    }
    baseline = qualified.get("w7_baseline")
    if baseline is None:
        raise ValueError("The baseline has not passed the distance gate. Increase --proof-seconds and rerun --stage screen --retry-inconclusive before simulation.")
    path = args.output / "selected_schedules.json"
    if args.only_ids:
        ids = list(dict.fromkeys(["w7_baseline"] + [item.strip() for item in args.only_ids.split(",") if item.strip()]))
        unavailable = [item for item in ids if item not in qualified]
        if unavailable:
            raise ValueError(f"Requested candidate IDs are not qualified: {unavailable}")
        selected = [qualified[item] for item in ids]
    else:
        selected = [baseline]
        if path.exists():
            saved_ids = json.loads(path.read_text()).get("candidate_ids", [])
            selected += [qualified[item] for item in saved_ids if item in qualified and item != "w7_baseline"]
        selected = selected[:args.max_circuits]
        used = {item["id"] for item in selected}
        groups = defaultdict(list)
        for item in records:
            if item["id"] in qualified and item["id"] not in used:
                groups[item.get("family", "other")].append(item)
        # Selection is fixed before Monte Carlo and does not use observed LER.
        while len(selected) < args.max_circuits and any(groups.values()):
            for family in sorted(groups):
                if groups[family] and len(selected) < args.max_circuits:
                    selected.append(groups[family].pop(0))
        if len(selected) < args.max_circuits:
            print(f"Qualified circuits available: {len(selected)}/{args.max_circuits}. All available will be simulated.", flush=True)
    atomic_json(path, {"candidate_ids": [item["id"] for item in selected],
                       "selection_policy": "baseline first; retain previous choices; family-balanced manifest order; no Monte Carlo ranking",
                       "max_circuits": args.max_circuits,
                       "selected_utc": datetime.now(timezone.utc).isoformat(),
                       "allow_uncertified": args.allow_uncertified})
    if len(selected) == 1:
        raise ValueError("Only the baseline qualified; no circuit comparison is possible yet. Increase --candidates or --proof-seconds, then rerun.")
    return selected


def run(args):
    if args.stage == "plot":
        from plot_scan import plot_scan
        plot_scan(args.output)
        return
    import fcntl
    args.output.mkdir(parents=True, exist_ok=True)
    lock = (args.output / "RUNNING.lock").open("a+", encoding="utf-8")
    try:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError(f"Another scan is writing {args.output}; wait for it to exit or use a different output.")
    with lock:
        lock.seek(0)
        lock.truncate()
        lock.write(f"pid={os.getpid()}\n")
        lock.flush()
        manifest = ensure_manifest(args)
        if args.stage == "generate":
            return
        from screening import screen_candidates
        # Cached complete audits are reused. An explicit simulation-only run
        # still checks their identity and distance gate before using them.
        records = screen_candidates(
            manifest=manifest, output=args.output,
            workers=min(args.distance_workers, args.workers), min_distance=6,
            heuristic_seconds=args.distance_seconds, proof_seconds=args.proof_seconds,
            allow_uncertified=args.allow_uncertified,
            retry_inconclusive=args.retry_inconclusive,
        )
        if args.stage == "screen":
            return
        selected = choose_candidates(records, args)
        print(f"Selected {len(selected)} circuits. Same W7 code, R=7, idle=0, joint X/Z, fixed beam=2.", flush=True)
        print(f"p order: {args.p_values}; full-wave stopping; relative 1sigma target={args.relative_error:.0%}.", flush=True)
        from simulation import simulate_point
        candidate_map = {item["id"]: item for item in manifest["candidates"]}
        for p in args.p_values:
            for item in selected:
                result = simulate_point(
                    candidate=candidate_map[item["id"]], audit=item["audit"], p=p,
                    output=args.output, workers=args.workers,
                    min_shots=args.min_shots, max_shots=args.max_shots,
                    relative_error=args.relative_error, batch_size=args.batch_size,
                    run_seed=args.sample_seed,
                )
                if result.get("status") == "interrupted":
                    print("Interrupted. Saved batches will resume with the same command.", flush=True)
                    return
        from plot_scan import plot_scan
        plot_scan(args.output)
        print(f"Finished. Reports and plots: {args.output}", flush=True)


def main():
    signal.signal(signal.SIGINT, signal.default_int_handler)
    args = parse_args()
    try:
        run(args)
    except KeyboardInterrupt:
        print("Interrupted; completed distance audits and simulation batches are saved.", file=sys.stderr, flush=True)
        return 130
    except (ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
