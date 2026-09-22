"""Plot independent confirmation data, never the historical selection counts.

Rates are an equivalent normalization of the joint block probability, not
measured individual-logical-qubit marginals. Intervals are nominal Wilson
z=1 intervals transformed monotonically; sequential stopping is not an
anytime-valid coverage guarantee.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def equivalent_rate(probability, k=12, rounds=7):
    if probability == 1:
        return 1.0
    return -math.expm1(math.log1p(-probability) / (k * rounds))


def wilson_interval(n, failures):
    if not n:
        return None, None
    p = failures / n
    denominator = 1 + 1 / n
    center = (p + .5 / n) / denominator
    radius = math.sqrt(p * (1 - p) / n + .25 / n**2) / denominator
    return max(0., center - radius), min(1., center + radius)


def _distance_bounds(record):
    if not record:
        return None, None
    record = record.get("audit", record)
    return record.get("lower_bound"), record.get("upper_bound")


def normalize_distance_records(records):
    if records is None:
        return {}
    if isinstance(records, list):
        return {record.get("id", record.get("candidate_id")): record for record in records}
    if isinstance(records, dict):
        for name in ("records", "candidates", "results"):
            if isinstance(records.get(name), list):
                return normalize_distance_records(records[name])
        return records
    raise ValueError("Distance records must be a mapping or list of candidate records")


def validate_distance_record(candidate, record):
    """Bind a distance label to the current complete reference circuit."""
    if not record:
        return {}
    from circuits import build_circuit
    audit = record.get("audit", record)
    claimed = audit.get("circuit_sha256") or record.get("identity", {}).get("circuit_sha256")
    if not claimed:
        raise ValueError(f"Distance record for {candidate['id']} lacks a circuit hash; rerun the distance stage")
    circuit, _ = build_circuit("w7", .003, rounds=7, idle_scale=0, schedule=candidate["schedule"])
    observed = hashlib.sha256(str(circuit).encode("utf-8")).hexdigest()
    if claimed != observed:
        raise ValueError(f"Distance record circuit hash mismatch for {candidate['id']}; rerun the distance stage")
    return record


def distance_label(lower, upper):
    if lower is not None and lower == upper:
        return rf"$d_{{\rm circ}}={lower}$"
    if lower is not None and upper is not None:
        return rf"${lower}\leq d_{{\rm circ}}\leq {upper}$"
    if lower is not None:
        return rf"$d_{{\rm circ}}\geq {lower}$"
    if upper is not None:
        return rf"$d_{{\rm circ}}\leq {upper}$"
    return "distance pending"


def summarize_point(result, distance_record=None):
    counts = result["counts"]
    n, failures = int(counts["N"]), int(counts["Fjoint"])
    k, rounds = int(result["k"]), int(result["rounds"])
    p = failures / n if n else None
    lo, hi = wilson_interval(n, failures)
    q = equivalent_rate(p, k, rounds) if p is not None else None
    qlo = equivalent_rate(lo, k, rounds) if lo is not None else None
    qhi = equivalent_rate(hi, k, rounds) if hi is not None else None
    q_sigma = None
    if n and 0 < failures < n:
        sigma = math.sqrt(p * (1 - p) / n)
        q_sigma = (1 - p)**(1 / (k * rounds) - 1) * sigma / (k * rounds)
    ledger = result.get("ledger", {})
    wave_complete = not ledger.get("pending") and ledger.get("active_wave") is None
    completed = bool(result.get("point_done")) and wave_complete
    lower, upper = _distance_bounds(distance_record or result.get("audit"))
    return {
        "candidate_id": result["candidate_id"], "physical_p": result["p"],
        "N": n, "Fjoint": failures, "FX": counts["FX"], "FZ": counts["FZ"],
        "Fboth": counts["Fboth"], "decodefail": counts["decodefail"],
        "k": k, "rounds": rounds, "Pjoint": p,
        "equivalent_per_logical_per_round": q,
        "q_wilson_z1_lower": qlo, "q_wilson_z1_upper": qhi,
        "q_delta_method_sigma": q_sigma,
        "relative_joint_sigma": math.sqrt((1 - p) / failures) if failures else None,
        "status": result["status"], "completed": completed,
        "target_precision_reached": result["status"] == "precision_reached",
        "distance_lower_bound": lower, "distance_upper_bound": upper,
        "run_seed": result.get("identity", {}).get("run_seed"),
        "experiment_id": result["experiment_id"],
        "source_result": result.get("result_path"),
    }


def two_point_slopes(rows):
    """Adjacent completed positive-probability pairs; no fit to partial points."""
    groups = {}
    for row in rows:
        groups.setdefault(row["candidate_id"], []).append(row)
    out = []
    for candidate, points in groups.items():
        points.sort(key=lambda item: item["physical_p"])
        for low, high in zip(points, points[1:]):
            entry = {
                "candidate_id": candidate, "p_low": low["physical_p"],
                "p_high": high["physical_p"], "alpha": None,
                "nominal_delta_method_sigma": None,
                "both_targets_reached": low["target_precision_reached"] and high["target_precision_reached"],
                "distance_lower_bound": low["distance_lower_bound"],
                "distance_upper_bound": low["distance_upper_bound"],
                "status": "unavailable_unfinished_point",
            }
            if not low["completed"] or not high["completed"]:
                out.append(entry)
                continue
            if not all(0 < row["Fjoint"] < row["N"] for row in (low, high)):
                entry["status"] = "unavailable_boundary_probability"
                out.append(entry)
                continue
            denominator = math.log(high["physical_p"] / low["physical_p"])
            qlow = low["equivalent_per_logical_per_round"]
            qhigh = high["equivalent_per_logical_per_round"]
            entry.update(
                alpha=math.log(qhigh / qlow) / denominator,
                nominal_delta_method_sigma=math.hypot(
                    low["q_delta_method_sigma"] / qlow,
                    high["q_delta_method_sigma"] / qhigh) / abs(denominator),
                status="estimated_from_independent_confirmation_data",
            )
            out.append(entry)
    return out


def _write_csv(path, rows):
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(output, distance_records=None, make_plot=True):
    output = Path(output)
    if distance_records is None and (output / "distance_records.json").exists():
        distance_records = json.loads((output / "distance_records.json").read_text(encoding="utf-8"))
    distance_records = normalize_distance_records(distance_records)
    rows = []
    for path in sorted((output / "points").glob("*/result.json")):
        result = json.loads(path.read_text(encoding="utf-8"))
        record = distance_records.get(result["candidate_id"], result.get("audit", {}))
        if record:
            validate_distance_record(result["candidate"], record)
            # A matching current reference proof must also describe the old
            # simulation circuit, including any local circuit-builder edits.
            from circuits import build_circuit
            actual, _ = build_circuit("w7", result["p"], rounds=result["rounds"],
                                      idle_scale=result["identity"]["idle_scale"],
                                      schedule=result["candidate"]["schedule"])
            actual_hash = hashlib.sha256(str(actual).encode("utf-8")).hexdigest()
            if actual_hash != result["identity"]["circuit_sha256"]:
                raise ValueError(f"Historical simulation circuit differs from the current builder: {path}")
        rows.append(summarize_point(result, record))
    seen = set()
    for row in rows:
        key = row["candidate_id"], row["physical_p"]
        if key in seen:
            raise ValueError(f"Multiple independent experiments at {key}; use separate output folders")
        seen.add(key)
    rows.sort(key=lambda item: (item["candidate_id"], item["physical_p"]))
    slopes = two_point_slopes(rows)
    summary = {
        "schema_version": 1, "normalization": "q=1-(1-Fjoint/N)^(1/(k*R)); equivalent rate, not marginal",
        "data_policy": "Independent confirmation only; historical selection counts are never pooled.",
        "intervals": "Nominal Wilson z=1, monotonically transformed to q; not anytime-valid intervals.",
        "slope_uncertainty": "Delta method, independent point counts; nominal only under sequential stopping.",
        "rows": rows, "two_point_slopes": slopes,
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "simulation_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    _write_csv(output / "simulation_summary.csv", rows)
    _write_csv(output / "two_point_slopes.csv", slopes)
    if make_plot:
        plot_summary(summary, output)
    return summary


def plot_summary(summary, output):
    rows = [row for row in summary["rows"] if row["completed"] and row["N"] > 0]
    if not rows:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    groups = {}
    for row in rows:
        groups.setdefault(row["candidate_id"], []).append(row)
    figure, ax = plt.subplots(figsize=(8.2, 7.3))
    colors = plt.get_cmap("tab10")
    for index, (candidate, points) in enumerate(groups.items()):
        points.sort(key=lambda row: row["physical_p"])
        color = "#a84d49" if candidate == "w7_baseline" else colors(index % 10)
        first = points[0]
        label = candidate + "  " + distance_label(first["distance_lower_bound"], first["distance_upper_bound"])
        positives = [row for row in points if row["Fjoint"] > 0]
        if positives:
            ax.plot([row["physical_p"] for row in positives],
                    [row["equivalent_per_logical_per_round"] for row in positives],
                    color=color, linewidth=1.3, alpha=.8)
        labelled = False
        for row in points:
            q, lo, hi = (row["equivalent_per_logical_per_round"], row["q_wilson_z1_lower"], row["q_wilson_z1_upper"])
            if row["Fjoint"] == 0:
                ax.errorbar(row["physical_p"], hi, yerr=.28 * hi, uplims=True,
                            color=color, fmt="v", capsize=3, mfc="white",
                            label=label if not labelled else None)
            else:
                ax.errorbar(row["physical_p"], q, yerr=[[q - lo], [hi - q]],
                            color=color, fmt="o", markersize=5.5, capsize=3,
                            mfc=color if row["target_precision_reached"] else "white",
                            label=label if not labelled else None)
            labelled = True
    ax.set_xscale("log")
    ax.set_yscale("log")
    unique_p = sorted({row["physical_p"] for row in rows})
    ax.set_xticks(unique_p)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, position: f"{100 * value:g}%"))
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda value, position: ""))
    ax.set_xlabel("Physical error probability, p")
    ax.set_ylabel("Equivalent logical failure rate per logical qubit per noisy round")
    ax.set_title("Weight-7: independent confirmation", pad=12)
    ax.grid(True, which="major", alpha=.22)
    ax.legend(loc="best", fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    caption = (
        "All circuits: [[56,12,7]], R=7 noisy rounds plus two ideal closing rounds; idle noise=0.\n"
        "Same-shot joint X/Z, full correlated DEM, Tesseract beam=2. q=1-(1-Fjoint/N)^(1/(12 R)).\n"
        "q is an equivalent normalization of the block rate, not a measured per-qubit marginal.\n"
        "Bars: nominal Wilson 1-sigma intervals. Open symbols: precision target not reached.\n"
        "Downward arrows: zero-failure upper limits. Only finished points are shown.\n"
        "Historical selection data are excluded. Sequential stopping does not ensure nominal coverage."
    )
    figure.text(.12, .035, caption, va="bottom", ha="left", fontsize=8.1, linespacing=1.35)
    figure.subplots_adjust(left=.13, right=.97, top=.91, bottom=.26)
    output = Path(output)
    figure.savefig(output / "independent_confirmation.png", dpi=220)
    figure.savefig(output / "independent_confirmation.pdf")
    plt.close(figure)
    (output / "figure_caption.txt").write_text(caption + "\n", encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results"),
                        help="Raw results directory, or figure destination when --summary is used")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--distance-records", type=Path,
                        help="Distance records to validate against raw point circuits")
    source.add_argument("--summary", type=Path,
                        help="Render an existing simulation_summary.json using only matplotlib; "
                             "stored statistics and distance labels are not revalidated")
    args = parser.parse_args(argv)
    if args.summary is not None:
        summary = json.loads(args.summary.read_text(encoding="utf-8"))
        args.output.mkdir(parents=True, exist_ok=True)
        plot_summary(summary, args.output)
        print(f"Rendered {args.output} from {args.summary}; "
              "stored statistics and distance labels were not revalidated")
        return
    records = json.loads(args.distance_records.read_text(encoding="utf-8")) if args.distance_records else None
    summary = write_outputs(args.output, records)
    print(f"Updated {args.output}: {len(summary['rows'])} points, {len(summary['two_point_slopes'])} adjacent slope pairs")


if __name__ == "__main__":
    main()
