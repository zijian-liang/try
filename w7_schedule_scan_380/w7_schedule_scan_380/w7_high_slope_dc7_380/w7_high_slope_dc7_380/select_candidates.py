"""Rank old completed W7 points and select unresolved high-slope schedules.

The two-point exponent is descriptive, not a circuit-distance certificate.
Only ``point_done`` records enter the ranking; provisional progress is ignored.
All computations here use Python's standard library.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re


DATA_DIR = Path(__file__).resolve().parent / "data"
POINT_RE = re.compile(
    r"\]\s+(?P<id>\S+)\s+p=(?P<p>[0-9.eE+-]+)\s+"
    r"N=(?P<N>[\d,]+)\s+Fjoint=(?P<Fjoint>[\d,]+)\s+"
    r"FX=(?P<FX>[\d,]+)\s+FZ=(?P<FZ>[\d,]+)\s+"
    r"Fboth=(?P<Fboth>[\d,]+)\s+decodefail=(?P<decodefail>[\d,]+)"
)
BOUND_RE = re.compile(
    r"^(?P<id>\S+):\s+(?:qualified|rejected|inconclusive)\s+"
    r"\[(?P<lower>\d+),\s*(?P<upper>\d+)\]"
)


def schedule_hash(schedule):
    return hashlib.sha256(
        json.dumps(schedule, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def parse_log(log_text):
    """Return latest completed cumulative points plus historical bound labels.

    Repeated identical lines are harmless. A later completed record can extend
    an earlier one, but regressing counts are rejected rather than combined.
    Historical distance labels remain unverified assertions until rerun by the
    new distance audit; ``certified=True`` in a log alone is not a proof file.
    """
    points, bounds = {}, {}
    for number, line in enumerate(log_text.splitlines(), 1):
        b = BOUND_RE.search(line)
        if b:
            bounds[b["id"]] = [int(b["lower"]), int(b["upper"])]
        if " point_done " not in line:
            continue
        m = POINT_RE.search(line)
        if not m:
            raise ValueError(f"Unrecognized point_done record on line {number}")
        point = {field: int(m[field].replace(",", ""))
                 for field in ("N", "Fjoint", "FX", "FZ", "Fboth", "decodefail")}
        point.update(p=float(m["p"]), source_line=number)
        if not (point["N"] > 0 and 0 <= point["Fjoint"] <= point["N"]):
            raise ValueError(f"Invalid shot/failure counts on line {number}")
        if point["Fjoint"] != point["FX"] + point["FZ"] - point["Fboth"] + point["decodefail"]:
            raise ValueError(f"Joint/X/Z failure counts disagree on line {number}")
        if not (0 <= point["decodefail"] <= point["Fjoint"]):
            raise ValueError(f"Invalid decoder failure count on line {number}")
        if not (0 <= point["Fboth"] <= min(point["FX"], point["FZ"])):
            raise ValueError(f"Invalid overlap counts on line {number}")
        key = (m["id"], point["p"])
        prior = points.get(key)
        if prior is not None:
            fields = ("N", "Fjoint", "FX", "FZ", "Fboth", "decodefail")
            if any(point[field] < prior[field] for field in fields):
                raise ValueError(f"Completed cumulative counts regress for {key}")
            if point["N"] == prior["N"] and any(
                    point[field] != prior[field] for field in fields):
                raise ValueError(f"Conflicting completed counts for {key}")
        points[key] = point
    return points, bounds


def equivalent_rate(point, k=12, rounds=7):
    """Equivalent per-logical-per-round rate and nominal delta-method SE.

    This normalization is NOT a directly measured per-qubit marginal.
    The SE treats the final counts as binomial and does not correct for
    adaptive stopping or selection among many candidate circuits.
    """
    n, failures = point["N"], point["Fjoint"]
    probability = failures / n
    if not 0 < probability < 1:
        raise ValueError("A finite two-point slope requires 0 < Fjoint < N")
    exponent = 1.0 / (k * rounds)
    log_survival = math.log1p(-probability)
    rate = -math.expm1(exponent * log_survival)
    derivative = exponent * math.exp((exponent - 1.0) * log_survival)
    standard_error = derivative * math.sqrt(probability * (1.0 - probability) / n)
    return rate, standard_error


def validate_manifest(manifest):
    """Verify schedule hashes and candidate IDs before using any schedule."""
    by_id = {}
    for candidate in manifest["candidates"]:
        identity = candidate["id"]
        if identity in by_id:
            raise ValueError(f"Duplicate manifest ID: {identity}")
        fingerprint = schedule_hash(candidate["schedule"])
        if fingerprint != candidate["schedule_sha256"]:
            raise ValueError(f"Schedule hash mismatch for {identity}")
        if identity != "w7_baseline" and not identity.endswith("_" + fingerprint[:12]):
            raise ValueError(f"Candidate ID does not match schedule hash: {identity}")
        schedule = candidate["schedule"]
        xs, zs = schedule["x"], schedule["z"]
        if len(xs) != len(zs):
            raise ValueError(f"Unequal X/Z schedule depths for {identity}")
        if any(x is None and z is None for x, z in zip(xs, zs)):
            raise ValueError(f"Empty CNOT layer for {identity}")
        if any(sorted(x for x in schedule[sector] if x is not None) != list(range(7))
               for sector in ("x", "z")):
            raise ValueError(f"Wrong direction coverage for {identity}")
        if any(x is not None and z is not None and (x < 3) == (z < 4)
               for x, z in zip(xs, zs)):
            raise ValueError(f"Colliding CNOT layer for {identity}")
        tx, tz = [xs.index(i) for i in range(7)], [zs.index(i) for i in range(7)]
        if not all((tx[a] < tz[b]) == (tx[b + 3] < tz[a + 4])
                   for a in range(3) for b in range(4)):
            raise ValueError(f"Invalid extraction parity for {identity}")
        by_id[identity] = candidate
    return by_id


def build_selection(log_text, manifest, top=3, p_low=0.0025, p_high=0.003,
                    k=12, rounds=7, include_baseline=True):
    """Return ``(selected_manifest, ranking_rows)`` from completed old points.

    Select the largest ``top`` slopes among prior [6,7] candidates and add
    baseline as a control. Bounds are re-audited by the caller, never inferred
    from the slope. The baseline comes first in the selected manifest.
    """
    if not isinstance(top, int) or isinstance(top, bool) or top < 1:
        raise ValueError("top must be a positive integer")
    if not 0 < p_low < p_high < 1 or k < 1 or rounds < 1:
        raise ValueError("Require 0 < p_low < p_high < 1 and positive k, rounds")
    candidates = validate_manifest(manifest)
    points, bounds = parse_log(log_text)
    names = sorted({identity for identity, p in points if p in (p_low, p_high)})
    rows = []
    for identity in names:
        if (identity, p_low) not in points or (identity, p_high) not in points:
            continue
        if identity not in candidates:
            raise ValueError(f"Measured schedule missing from manifest: {identity}")
        low, high = points[identity, p_low], points[identity, p_high]
        if low["decodefail"] or high["decodefail"]:
            raise ValueError(f"Decode failures in old points for {identity}; inspect explicitly")
        q_low, se_low = equivalent_rate(low, k, rounds)
        q_high, se_high = equivalent_rate(high, k, rounds)
        divisor = math.log(p_high / p_low)
        alpha = math.log(q_high / q_low) / divisor
        alpha_se = math.hypot(se_low / q_low, se_high / q_high) / divisor
        historical_bounds = bounds.get(identity, [None, None])
        row = {
            "id": identity, "alpha": alpha, "alpha_nominal_1sigma": alpha_se,
            "prior_distance_lower": historical_bounds[0],
            "prior_distance_upper": historical_bounds[1],
            "is_baseline": identity == "w7_baseline",
            "eligible_unresolved_6_to_7": historical_bounds == [6, 7],
            "p_low": p_low, "p_high": p_high,
            "k": k, "rounds": rounds,
            "q_low": q_low, "q_high": q_high,
            "q_low_nominal_1sigma": se_low, "q_high_nominal_1sigma": se_high,
        }
        for suffix, point in (("low", low), ("high", high)):
            for field in ("N", "Fjoint", "FX", "FZ", "Fboth", "decodefail", "source_line"):
                row[f"{field}_{suffix}"] = point[field]
        rows.append(row)
    rows.sort(key=lambda row: (-row["alpha"], row["id"]))
    eligible = [row for row in rows if row["eligible_unresolved_6_to_7"] and not row["is_baseline"]]
    if len(eligible) < top:
        raise ValueError(f"Requested top={top}, but only {len(eligible)} measured [6,7] candidates")
    selected_ids = [row["id"] for row in eligible[:top]]
    if include_baseline:
        if "w7_baseline" not in {row["id"] for row in rows}:
            raise ValueError("Baseline control has no complete pair of points")
        selected_ids.insert(0, "w7_baseline")
    by_row = {}
    for rank, row in enumerate(rows, 1):
        row.update(rank=rank, selected=row["id"] in selected_ids)
        by_row[row["id"]] = row
    selected = []
    for identity in selected_ids:
        candidate = dict(candidates[identity])
        row = by_row[identity]
        candidate.update(
            historical_distance_bounds=bounds.get(identity),
            historical_distance_bounds_source="Old log only; fresh audit required",
            historical_slope={"alpha": row["alpha"],
                              "nominal_1sigma": row["alpha_nominal_1sigma"],
                              "p_low": p_low, "p_high": p_high,
                              "rank_among_measured": row["rank"]},
            historical_points=[points[identity, p_high], points[identity, p_low]],
            selection_reason="baseline control" if identity == "w7_baseline"
                             else "descending two-point slope among prior [6,7] candidates",
        )
        selected.append(candidate)
    result = {
        "schema_version": 1,
        "selection": {
            "top": top, "baseline_included": include_baseline,
            "p_low": p_low, "p_high": p_high, "k": k, "rounds": rounds,
            "normalization": "q=1-(1-Fjoint/N)^(1/(k*rounds))",
            "records": "point_done only; latest completed cumulative record per point",
            "uncertainty": "nominal binomial delta-method standard error; no stopping/selection adjustment",
            "warning": "Two-point slope is not a distance certificate. Selection on noisy old data requires independent validation.",
            "log_sha256": hashlib.sha256(log_text.encode()).hexdigest(),
            "source_manifest_sha256": hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest(),
        },
        "candidates": selected,
        "all_measured_ranking": rows,
    }
    return result, rows


def write_selection(output_dir, selected_manifest, ranking_rows):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "candidates.json").write_text(
        json.dumps(selected_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with (output_dir / "ranking.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ranking_rows[0]))
        writer.writeheader()
        writer.writerows(ranking_rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, default=DATA_DIR / "old_scan.log")
    parser.add_argument("--manifest", type=Path, default=DATA_DIR / "measured_candidates.json")
    parser.add_argument("--top", type=int, default=3)
    parser.add_argument("--p-low", type=float, default=0.0025)
    parser.add_argument("--p-high", type=float, default=0.003)
    parser.add_argument("--k", type=int, default=12)
    parser.add_argument("--rounds", type=int, default=7)
    parser.add_argument("--no-baseline", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=DATA_DIR)
    args = parser.parse_args()
    selection, rows = build_selection(
        args.log.read_text(encoding="utf-8"),
        json.loads(args.manifest.read_text(encoding="utf-8")),
        top=args.top, p_low=args.p_low, p_high=args.p_high, k=args.k,
        rounds=args.rounds, include_baseline=not args.no_baseline)
    write_selection(args.output_dir, selection, rows)
    for row in rows:
        marker = "*" if row["selected"] else " "
        print(f"{marker} {row['rank']:2d} {row['id']:30s} "
              f"alpha={row['alpha']:.4f} +/- {row['alpha_nominal_1sigma']:.4f} "
              f"prior d=[{row['prior_distance_lower']},{row['prior_distance_upper']}]")
    print(f"Saved candidates.json and ranking.csv in {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
