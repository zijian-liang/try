"""Independent confirmation of selected W7 circuits with Tesseract beam=2.

Defaults keep the original p=.3%, .25%, R=7 and zero idle noise. The new RNG
seed separates confirmation samples from the historical candidate selection.
Stop precision is evaluated only after complete precommitted shot waves.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import json
import math
import os
from pathlib import Path
from importlib.metadata import version

from simulation import simulate_point, _atomic_json
from plot_results import write_outputs, normalize_distance_records, validate_distance_record


@contextmanager
def _writer_lock(output):
    """POSIX lock is released automatically if the coordinator exits/crashes."""
    import fcntl
    path = Path(output) / ".simulation.lock"
    with path.open("a+", encoding="utf-8") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"Another simulation writer uses {output}") from exc
        stream.seek(0)
        stream.truncate()
        stream.write(str(os.getpid()) + "\n")
        stream.flush()
        try:
            yield
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def run_selected(candidates, output, workers=380, p_values=(.003, .0025),
                 relative_error=.20, max_shots=100_000_000, min_shots=1000,
                 batch_size=64, seed=20260923, distance_records=None):
    """Return summary dict and save point checkpoints, raw counts, plots/slopes.

    distance_records maps candidate ID to its latest audit or {'audit': audit}.
    Distance uncertainty never blocks simulation. Use a separate directory for
    a different candidate set or RNG seed; stop precision can change on resume.
    """
    candidates = list(candidates)
    if not candidates:
        raise ValueError("No candidates selected")
    if min(workers, min_shots, max_shots, batch_size) < 1 or min_shots > max_shots:
        raise ValueError("Positive workers/shots/batch size and min_shots <= max_shots are required")
    if not math.isfinite(relative_error) or not 0 < relative_error < 1 or seed < 0:
        raise ValueError("Require finite 0<relative_error<1 and nonnegative seed")
    identifiers = [item["id"] for item in candidates]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("Duplicate candidate IDs")
    p_values = sorted(set(float(value) for value in p_values), reverse=True)
    if not p_values or not all(math.isfinite(value) and 0 < value < .5 for value in p_values):
        raise ValueError("Require at least one finite physical error probability with 0<p<.5")
    output = Path(output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    if distance_records is None and (output / "distance_records.json").exists():
        distance_records = json.loads((output / "distance_records.json").read_text(encoding="utf-8"))
    distance_records = normalize_distance_records(distance_records)
    for candidate in candidates:
        if candidate["id"] in distance_records:
            validate_distance_record(candidate, distance_records[candidate["id"]])
    campaign_identity = {
        "purpose": "independent_confirmation_after_historical_slope_selection",
        "run_seed": seed, "rounds": 7, "idle_scale": 0.0,
        "candidates": [{"id": item["id"], "schedule": item["schedule"]}
                       for item in sorted(candidates, key=lambda value: value["id"])],
    }
    with _writer_lock(output):
        config_path = output / "simulation_config.json"
        if config_path.exists():
            prior = json.loads(config_path.read_text(encoding="utf-8"))
            if prior["identity"] != campaign_identity:
                raise RuntimeError("Candidate set/seed changed; choose a new output folder")
        for path in (output / "points").glob("*/result.json"):
            prior = json.loads(path.read_text(encoding="utf-8"))
            if prior["candidate_id"] not in identifiers or prior["identity"]["run_seed"] != seed:
                raise RuntimeError(f"Unrelated/selection data found in confirmation output: {path}")
            prior_versions = prior["identity"].get("versions", {})
            if any(version(package) != old_version for package, old_version in prior_versions.items()):
                raise RuntimeError("Simulation dependency versions changed; choose a new output folder")
        _atomic_json(config_path, {
            "identity": campaign_identity, "p_values": p_values,
            "workers": workers, "relative_error": relative_error,
            "min_shots": min_shots, "max_shots": max_shots, "batch_size": batch_size,
        })
        interrupted = False
        try:
            for physical_p in p_values:
                for candidate in candidates:
                    record = distance_records.get(candidate["id"], {})
                    audit = record.get("audit", record)
                    result = simulate_point(candidate, audit, physical_p, output,
                                            workers=workers, min_shots=min_shots,
                                            max_shots=max_shots, relative_error=relative_error,
                                            batch_size=batch_size, run_seed=seed)
                    write_outputs(output, distance_records)
                    if result["status"] == "interrupted":
                        interrupted = True
                        break
                if interrupted:
                    break
        finally:
            summary = write_outputs(output, distance_records)
        summary["interrupted"] = interrupted
        return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path,
                        help="Manifest list or object with a candidates/selected list")
    parser.add_argument("--candidate", action="append", default=[], help="Optional candidate ID filter; repeatable")
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--workers", type=int, default=380)
    parser.add_argument("--p", type=float, nargs="+", default=[.003, .0025])
    parser.add_argument("--relative-error", type=float, default=.20)
    parser.add_argument("--min-shots", type=int, default=1000)
    parser.add_argument("--max-shots", type=int, default=100_000_000)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20260923)
    parser.add_argument("--distance-records", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.candidates.read_text(encoding="utf-8"))
    candidates = manifest if isinstance(manifest, list) else manifest.get("candidates", manifest.get("selected", []))
    if args.candidate:
        candidates = [item for item in candidates if item["id"] in args.candidate]
        missing = set(args.candidate) - {item["id"] for item in candidates}
        if missing:
            raise ValueError(f"Unknown candidates: {sorted(missing)}")
    records = json.loads(args.distance_records.read_text()) if args.distance_records else None
    summary = run_selected(candidates, args.output, args.workers, args.p,
                           args.relative_error, args.max_shots, args.min_shots,
                           args.batch_size, args.seed, records)
    return 130 if summary.get("interrupted") else 0


if __name__ == "__main__":
    raise SystemExit(main())
