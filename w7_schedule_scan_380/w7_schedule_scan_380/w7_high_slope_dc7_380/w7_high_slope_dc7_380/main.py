#!/usr/bin/env python3
"""Select historical high-slope W7 schedules, decide d_circ=7, then resample."""
from __future__ import annotations

import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
import csv
import json
import math
from pathlib import Path
import sys
from contextlib import contextmanager

import settings as defaults

HERE = Path(__file__).resolve().parent


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


@contextmanager
def writer_lock(output):
    # flock releases itself after termination; a leftover filename is harmless.
    import fcntl
    output.mkdir(parents=True, exist_ok=True)
    with (output / "main.lock").open("a+") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"Another main.py is writing {output}; use a separate --output.") from exc
        try:
            yield
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("all", "select", "distance", "simulate", "plot"), default="all")
    parser.add_argument("--output", type=Path, default=defaults.OUTPUT_DIRECTORY)
    parser.add_argument("--log", type=Path, default=HERE / "data" / "old_scan.log")
    parser.add_argument("--manifest", type=Path, default=HERE / "data" / "measured_candidates.json")
    parser.add_argument("--top", type=int, default=defaults.TOP_CANDIDATES)
    parser.add_argument("--workers", type=int, default=defaults.SIMULATION_WORKERS)
    parser.add_argument("--distance-workers", type=int, default=defaults.DISTANCE_WORKERS)
    parser.add_argument("--screening-workers", type=int, default=defaults.SCREENING_WORKERS)
    parser.add_argument("--proof-seconds", type=float, default=defaults.SCREENING_PROOF_SECONDS)
    parser.add_argument("--heuristic-seconds", type=float, default=defaults.HEURISTIC_SECONDS)
    parser.add_argument("--sat-seconds", type=float, default=defaults.SAT_SECONDS_PER_TASK)
    parser.add_argument("--solver", choices=("minisat22", "glucose3", "glucose4"), default=defaults.SAT_SOLVER)
    parser.add_argument("--shards", type=int, default=defaults.DISTANCE_SHARDS)
    parser.add_argument("--retry-unknown", action="store_true", default=defaults.RETRY_UNKNOWN)
    parser.add_argument("--p", type=float, nargs="+", default=defaults.P_VALUES)
    parser.add_argument("--relative-error", type=float, default=defaults.RELATIVE_1SIGMA)
    parser.add_argument("--min-shots", type=int, default=defaults.MIN_SHOTS)
    parser.add_argument("--max-shots", type=int, default=defaults.MAX_SHOTS)
    parser.add_argument("--batch-size", type=int, default=defaults.BATCH_SIZE)
    parser.add_argument("--sample-seed", type=int, default=defaults.SAMPLE_SEED)
    args = parser.parse_args()
    if min(args.top, args.workers, args.distance_workers, args.screening_workers,
           args.shards, args.min_shots, args.max_shots, args.batch_size) < 1:
        parser.error("Counts and worker limits must be positive.")
    if args.min_shots > args.max_shots:
        parser.error("min-shots must not exceed max-shots.")
    if not math.isfinite(args.relative_error) or not 0 < args.relative_error < 1:
        parser.error("relative-error must be in (0,1).")
    if any(not math.isfinite(p) or not 0 < p < .5 for p in args.p):
        parser.error("p must be finite and in (0,0.5).")
    if any(not math.isfinite(s) or s < 0 for s in (args.proof_seconds, args.heuristic_seconds, args.sat_seconds)):
        parser.error("Time budgets must be finite and nonnegative; 0 disables a search.")
    if args.sample_seed < 0:
        parser.error("sample-seed must be nonnegative.")
    args.p = sorted(set(args.p), reverse=True)
    args.output = args.output.expanduser().resolve()
    return args


def select(args):
    from select_candidates import build_selection
    source_text = args.log.read_text(encoding="utf-8")
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    selected, rows = build_selection(source_text, manifest, top=args.top)
    path = args.output / "selected_candidates.json"
    if path.exists():
        old = json.loads(path.read_text(encoding="utf-8"))
        old_schedules = {c["id"]: c["schedule"] for c in old["candidates"]}
        for candidate in selected["candidates"]:
            if candidate["id"] in old_schedules and old_schedules[candidate["id"]] != candidate["schedule"]:
                raise ValueError("An existing candidate changed schedule. Use a new output directory.")
    atomic_json(path, selected)
    atomic_json(args.output / "historical_ranking.json", rows)
    if rows:
        with (args.output / "historical_ranking.csv").open("w", newline="", encoding="utf-8") as stream:
            fields = list(dict.fromkeys(key for row in rows for key in row))
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    print("Selected: " + ", ".join(c["id"] for c in selected["candidates"]), flush=True)
    print("Historical two-point slopes select candidates only; new samples remain independent.", flush=True)
    return selected["candidates"]


def run(args):
    with writer_lock(args.output):
        if args.stage == "plot":
            from plot_results import write_outputs
            saved = args.output / "distance_records.json"
            records = json.loads(saved.read_text(encoding="utf-8")) if saved.exists() else []
            write_outputs(args.output, distance_records={r["id"]: r for r in records})
            return 0
        candidates = select(args)
        if args.stage == "select":
            return 0
        records_path = args.output / "distance_records.json"
        records = []
        if args.stage in ("all", "distance"):
            from exact_distance import run_distance_campaign
            records = run_distance_campaign(
                candidates, args.output,
                workers=args.distance_workers, seconds_per_task=args.sat_seconds,
                shards=args.shards, screening_workers=args.screening_workers,
                proof_seconds=args.proof_seconds, heuristic_seconds=args.heuristic_seconds,
                retry_unknown=args.retry_unknown,
                solver_name=args.solver,
            )
            atomic_json(records_path, records)
            for record in records:
                audit = record.get("audit", record)
                print(f"Distance {record['id']}: [{audit.get('lower_bound')}, "
                      f"{audit.get('upper_bound')}], exact={audit.get('exact', False)}", flush=True)
            if args.stage == "distance":
                return 0
        elif records_path.exists():
            records = json.loads(records_path.read_text(encoding="utf-8"))
        print(f"Confirmation: R=7, idle=0, joint X/Z, beam=2; p={args.p}; "
              f"relative 1sigma <= {args.relative_error:.0%}.", flush=True)
        from run_simulation import run_selected
        result = run_selected(
            candidates, args.output, workers=args.workers, p_values=args.p,
            relative_error=args.relative_error, min_shots=args.min_shots,
            max_shots=args.max_shots, batch_size=args.batch_size,
            seed=args.sample_seed, distance_records={r["id"]: r for r in records},
        )
        print(f"Results: {args.output}", flush=True)
        return 130 if result.get("interrupted") else 0


def main():
    args = parse_args()
    try:
        return run(args)
    except KeyboardInterrupt:
        print("Interrupted. Completed proofs and committed simulation batches are saved.", flush=True)
        return 130
    except (ValueError, RuntimeError, FileNotFoundError, ImportError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
