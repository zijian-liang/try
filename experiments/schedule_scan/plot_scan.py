#!/usr/bin/env python3
"""Read saved scan snapshots and plot candidate schedules without changing them."""
from __future__ import annotations

import argparse
import csv
import json
import math
import textwrap
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


FINAL_STATUSES = {"precision_reached", "max_shots"}
DEFAULT_P_VALUES = (0.003, 0.0025)
RAW_COUNTS = ("N", "Fjoint", "FX", "FZ", "Fboth", "decodefail", "lowconfidence",
              "invalid_syndrome", "decoder_exceptions")


def wilson(failures: int, shots: int) -> tuple[float, float]:
    """Nominal 68.27% Wilson score interval, z=1."""
    q = failures / shots
    denominator = 1 + 1 / shots
    centre = (q + 0.5 / shots) / denominator
    half = math.sqrt(q * (1 - q) / shots + 0.25 / shots**2) / denominator
    return max(0.0, centre - half), min(1.0, centre + half)


def per_round(q: float, rounds: int) -> float:
    if q >= 1:
        return 1.0
    return -math.expm1(math.log1p(-q) / rounds)


def point_statistics(record: dict) -> dict:
    n, f, rounds = int(record["counts"]["N"]), int(record["counts"]["Fjoint"]), int(record["rounds"])
    if n == 0:
        return {key: None for key in ("P_fail_shot", "p_round", "p_round_low_1sigma", "p_round_high_1sigma", "relative_sigma")}
    q = f / n
    lo, hi = wilson(f, n)
    return {"P_fail_shot": q, "p_round": per_round(q, rounds),
            "p_round_low_1sigma": per_round(lo, rounds), "p_round_high_1sigma": per_round(hi, rounds),
            "relative_sigma": math.sqrt(q * (1 - q) / n) / q if q > 0 else None}


def precision_reached(record: dict) -> bool:
    if record.get("status") == "precision_reached":
        return True
    stop = record.get("stop_settings", {})
    relative = point_statistics(record)["relative_sigma"]
    target = stop.get("relative_error")
    minimum = stop.get("min_shots", 0)
    return target is not None and relative is not None and int(record["counts"]["N"]) >= int(minimum) and relative <= float(target)


def _audit(record: dict) -> dict:
    value = record.get("_effective_audit") or record.get("audit") or record.get("candidate", {}).get("audit") or {}
    if not isinstance(value, dict):
        return {}
    if "codes" in value:
        value = value["codes"].get("w7", {})
    return value if isinstance(value, dict) else {}


def selected_candidates(results: Path) -> set[str] | None:
    path = results / "selected_schedules.json"
    if not path.exists():
        return None
    selection = json.loads(path.read_text(encoding="utf-8"))
    return {str(value) for value in selection["candidate_ids"]}


def update_distance_labels(records: list[dict], results: Path) -> None:
    """Use a stronger, applicable latest report; retain the point's original audit."""
    for record in records:
        path = results / "distance" / (record["candidate_id"] + ".json")
        if not path.exists():
            continue
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            warnings.warn(f"Skipped unreadable distance snapshot {path}: {exc}")
            continue
        newer = report.get("audit", report)
        if not isinstance(newer, dict):
            continue
        if "codes" in newer:
            newer = newer["codes"].get("w7", {})
        if not isinstance(newer, dict):
            continue
        expected_schedule = record.get("candidate", {}).get("schedule", record.get("metadata", {}).get("schedule"))
        idle = newer.get("idle_scale", newer.get("noise", {}).get("idle_scale"))
        if newer.get("rounds") != record["rounds"] or idle != 0 or newer.get("schedule") != expected_schedule:
            warnings.warn(f"Distance snapshot not applicable to {record['candidate_id']}; using saved point audit")
            continue
        old = _audit(record)
        old_lo, old_hi = old.get("lower_bound"), old.get("upper_bound")
        new_lo, new_hi = newer.get("lower_bound"), newer.get("upper_bound")
        lower_values = [x for x in (old_lo, new_lo) if x is not None]
        upper_values = [x for x in (old_hi, new_hi) if x is not None]
        if lower_values and upper_values and max(lower_values) > min(upper_values):
            raise ValueError(f"Conflicting distance certificates for {record['candidate_id']}")
        stronger_lo = new_lo is not None and (old_lo is None or new_lo > old_lo)
        stronger_hi = new_hi is not None and (old_hi is None or new_hi < old_hi)
        if stronger_lo or stronger_hi:
            effective = dict(newer)
            effective["lower_bound"] = max(lower_values) if lower_values else None
            effective["upper_bound"] = min(upper_values) if upper_values else None
            effective["exact"] = bool(lower_values and upper_values and max(lower_values) == min(upper_values))
            record["_effective_audit"] = effective
            record["_distance_update_source"] = str(path)


def _is_baseline(record: dict) -> bool:
    return bool(record.get("candidate", {}).get("is_baseline", record.get("is_baseline", False)))


def _candidate_id(record: dict) -> str:
    return str(record.get("candidate_id") or record.get("candidate", {}).get("id") or "")


def read_records(results: Path) -> list[dict]:
    paths = [results] if results.is_file() else sorted((results / "points").glob("*/result.json"))
    records, seen, common = [], {}, None
    for path in paths:
        try:
            r = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            warnings.warn(f"Skipped unreadable snapshot {path}: {exc}")
            continue
        candidate_id = _candidate_id(r)
        if not candidate_id or "counts" not in r:
            continue
        n = int(r["counts"].get("N", 0))
        counts = {key: int(r["counts"].get(key, 0)) for key in RAW_COUNTS}
        if any(not 0 <= counts[key] <= n for key in RAW_COUNTS):
            raise ValueError(f"Invalid count in {path}")
        if counts["Fjoint"] != counts["FX"] + counts["FZ"] - counts["Fboth"] + counts["decodefail"]:
            raise ValueError(f"Same-shot count identity failed in {path}")
        if counts["Fboth"] > min(counts["FX"], counts["FZ"]):
            raise ValueError(f"Invalid Fboth in {path}")
        r["rounds"] = int(r.get("rounds", r.get("metadata", {}).get("rounds", 7)))
        r["k"] = int(r.get("k", r.get("metadata", {}).get("k", 12)))
        if r["rounds"] != 7 or r["k"] != 12:
            raise ValueError(f"This scan comparison requires R=7 and k=12: {path}")
        identity = r.get("identity", {})
        protocol = {key: identity[key] for key in ("protocol", "idle_scale", "decoder", "decoder_settings", "versions", "count_rule") if key in identity}
        if protocol:
            signature = json.dumps(protocol, sort_keys=True)
            if common is not None and signature != common:
                raise ValueError("Mixed noise or decoder protocols; choose one scan results directory")
            common = signature
        key = (candidate_id, float(r["p"]))
        if key in seen:
            raise ValueError(f"Duplicate candidate/p result: {seen[key]} and {path}; choose one scan campaign")
        seen[key] = path
        r["candidate_id"], r["_path"] = candidate_id, str(path)
        r["_statistics"] = point_statistics(r)
        records.append(r)
    return records


def _p_key(value: float) -> float:
    return round(float(value), 12)


def _distance_label(record: dict) -> str:
    audit = _audit(record)
    lo, hi = audit.get("lower_bound"), audit.get("upper_bound")
    if audit.get("exact") and lo is not None and lo == hi:
        return rf"$d_{{\rm circ}}={int(lo)}$"
    if lo is not None and hi is not None:
        return rf"${int(lo)}\leq d_{{\rm circ}}\leq {int(hi)}$"
    if lo is not None:
        return rf"$d_{{\rm circ}}\geq {int(lo)}$"
    if hi is not None:
        return rf"$d_{{\rm circ}}\leq {int(hi)}$"
    return r"$d_{\rm circ}$ unverified"


def _short_id(value: str) -> str:
    return value if len(value) <= 24 else value[:12] + "..." + value[-8:]


def _included(record: dict, include_incomplete: bool) -> bool:
    return int(record["counts"]["N"]) > 0 and (record.get("status") in FINAL_STATUSES or include_incomplete)


def _baseline_records(records: list[dict]) -> dict:
    baseline_ids = {_candidate_id(r) for r in records if _is_baseline(r)}
    if len(baseline_ids) > 1:
        raise ValueError("Multiple baseline candidate IDs; select one scan campaign")
    return {_p_key(r["p"]): r for r in records if _is_baseline(r)}


def write_tables(records: list[dict], out: Path, include_incomplete: bool, p_values: tuple[float, ...], selection: set[str] | None) -> None:
    plotted_p = {_p_key(p) for p in p_values}
    baseline = _baseline_records([r for r in records if selection is None or r["candidate_id"] in selection])
    rows = []
    for r in records:
        stats, candidate, audit = r["_statistics"], r.get("candidate", {}), _audit(r)
        b = baseline.get(_p_key(r["p"]))
        is_selected = selection is None or r["candidate_id"] in selection
        eligible = is_selected and _included(r, include_incomplete) and _p_key(r["p"]) in plotted_p
        b_eligible = b is not None and _included(b, include_incomplete)
        b_rate = b["_statistics"]["p_round"] if b_eligible else None
        rate = stats["p_round"]
        ratio = rate / b_rate if b_rate is not None and b_rate > 0 and rate is not None else None
        row = {"candidate_id": r["candidate_id"], "family": candidate.get("family", ""),
               "is_baseline": _is_baseline(r), "p": float(r["p"]), "rounds": r["rounds"], "k": r["k"],
               "status": r.get("status", ""), "included_in_plot": eligible,
               "selected_for_campaign": is_selected,
               "precision_reached": precision_reached(r), "zero_failures": int(r["counts"]["Fjoint"]) == 0,
               "distance_lower": audit.get("lower_bound"), "distance_upper": audit.get("upper_bound"),
               "distance_exact": audit.get("exact", False), **{key: r["counts"].get(key, 0) for key in RAW_COUNTS},
               "saved_point_audit": json.dumps(r.get("audit", {}), sort_keys=True),
               "distance_update_source": r.get("_distance_update_source", ""),
               **stats, "ratio_to_baseline": ratio,
               "baseline_N": b["counts"]["N"] if b else None,
               "baseline_Fjoint": b["counts"]["Fjoint"] if b else None,
               "baseline_p_round": b_rate, "baseline_status": b.get("status", "") if b else "missing",
               "logical_X": json.dumps(r["counts"].get("logical_X", [])),
               "logical_Z": json.dumps(r["counts"].get("logical_Z", [])),
               "schedule": json.dumps(candidate.get("schedule", r.get("metadata", {}).get("schedule", {})), sort_keys=True),
               "experiment_id": r.get("experiment_id", ""), "source": r["_path"]}
        rows.append(row)
    if not rows:
        return
    fields = list(rows[0])
    rows.sort(key=lambda row: (-row["p"], not row["is_baseline"], row["candidate_id"]))
    with (out / "comparison.csv").open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    rankings = []
    for p in p_values:
        subset = [row for row in rows if _p_key(row["p"]) == _p_key(p) and row["included_in_plot"]]
        subset.sort(key=lambda row: (row["zero_failures"], row["p_round"], row["candidate_id"]))
        rank = 0
        for row in subset:
            if not row["zero_failures"]:
                rank += 1
            note = "point-estimate ordering only; independently validate before choosing"
            if row["zero_failures"]:
                note = "zero-failure upper limit; no point-estimate rank"
            elif row["status"] not in FINAL_STATUSES:
                note = "provisional ongoing point; independently validate"
            elif not row["precision_reached"]:
                note = "shot cap reached without target precision; descriptive rank only"
            rankings.append({"rank_point_estimate": "" if row["zero_failures"] else rank,
                             "ranking_note": note, **row})
    with (out / "ranking.csv").open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["rank_point_estimate", "ranking_note", *fields])
        writer.writeheader()
        writer.writerows(rankings)


def plot_scan(results: str | Path, out: str | Path | None = None, *,
              include_incomplete: bool = False, p_values: tuple[float, ...] = DEFAULT_P_VALUES) -> dict:
    """Save comparison.pdf/png and raw/descriptive CSV tables from atomic JSON snapshots."""
    results = Path(results)
    out = Path(out) if out is not None else (results.parent if results.is_file() else results) / "plots"
    records = read_records(results)
    selection = selected_candidates(results) if results.is_dir() else None
    update_distance_labels(records, results if results.is_dir() else results.parent.parent.parent)
    out.mkdir(parents=True, exist_ok=True)
    write_tables(records, out, include_incomplete, p_values, selection)
    selected = [r for r in records if (selection is None or r["candidate_id"] in selection) and _included(r, include_incomplete) and _p_key(r["p"]) in {_p_key(p) for p in p_values}]
    if not selected:
        return {"plotted": False, "read_points": len(records), "message": "No completed nonempty points; simulation snapshots unchanged", "out": str(out)}
    baseline = _baseline_records([r for r in records if selection is None or r["candidate_id"] in selection])
    baseline_lower = next((_audit(r).get("lower_bound") for r in baseline.values() if _audit(r).get("lower_bound") is not None), 6)
    groups = {}
    for r in selected:
        groups.setdefault(r["candidate_id"], []).append(r)
    ids = sorted(groups, key=lambda name: (not any(_is_baseline(r) for r in groups[name]), name))
    labels, representatives = [], {}
    for candidate_id in ids:
        candidate_records = groups[candidate_id]
        schedules = {json.dumps(r.get("candidate", {}).get("schedule", r.get("metadata", {}).get("schedule")), sort_keys=True) for r in candidate_records}
        if len(schedules) != 1:
            raise ValueError(f"Candidate ID {candidate_id} refers to multiple schedules")
        representative = max(candidate_records, key=lambda r: (_audit(r).get("lower_bound") or -1, bool(_audit(r).get("exact"))))
        representatives[candidate_id] = representative
        prefix = "Baseline: " if _is_baseline(representative) else ""
        labels.append(prefix + _short_id(candidate_id) + "   " + _distance_label(representative))
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9.5,
                         "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42})
    height = max(4.4, 2.55 + 0.40 * len(ids))
    fig, axes = plt.subplots(1, len(p_values), figsize=(12.4, height), sharey=True, squeeze=False)
    axes = axes[0]
    lookup = {(r["candidate_id"], _p_key(r["p"])): r for r in selected}
    has_uncertified, has_unmet, has_zero = False, False, False
    for ax, p in zip(axes, p_values):
        for row, candidate_id in enumerate(ids):
            r = lookup.get((candidate_id, _p_key(p)))
            if r is None:
                ax.text(0.025, row, "not completed", transform=ax.get_yaxis_transform(), color="#8B9299", fontsize=8)
                continue
            s = r["_statistics"]
            n, f = int(r["counts"]["N"]), int(r["counts"]["Fjoint"])
            lower = _audit(r).get("lower_bound")
            uncertified = lower is None or (baseline_lower is not None and lower < baseline_lower)
            complete_precision = r.get("status") in FINAL_STATUSES and precision_reached(r)
            has_uncertified |= uncertified
            has_unmet |= not complete_precision
            has_zero |= f == 0
            color = "#B84336" if _is_baseline(r) else ("#B57A23" if uncertified else "#2B6B9C")
            marker = "D" if uncertified else "o"
            face = color if complete_precision else "white"
            if f > 0:
                y, lo, hi = s["p_round"], s["p_round_low_1sigma"], s["p_round_high_1sigma"]
                ax.errorbar(y, row, xerr=[[max(0.0, y-lo)], [max(0.0, hi-y)]], fmt=marker,
                            ms=5.5, color=color, markerfacecolor=face, markeredgewidth=1.1,
                            capsize=2.5, elinewidth=1.1, zorder=3)
            else:
                hi = s["p_round_high_1sigma"]
                ax.plot(hi, row, marker="<", markerfacecolor="white", color=color, markersize=6)
                ax.plot(hi / 1.7, row, marker=".", alpha=0)
                ax.annotate("", xy=(hi / 1.65, row), xytext=(hi, row), arrowprops={"arrowstyle": "->", "color": color, "lw": 1.1})
        b = baseline.get(_p_key(p))
        if b and _included(b, include_incomplete) and b["_statistics"]["p_round"] > 0:
            ax.axvline(b["_statistics"]["p_round"], ls="--", lw=0.9, color="#B84336", alpha=0.65)
        ax.set_xscale("log")
        ax.set_xlabel("Joint block failure rate per noisy round")
        ax.set_title(rf"$p={p:g}$", fontsize=12)
        ax.set_yticks(range(len(ids)), labels)
        ax.set_ylim(len(ids) - 0.45, -0.55)
        ax.grid(axis="x", which="major", color="#DCE1E5", lw=0.7)
        ax.grid(axis="x", which="minor", color="#EDF0F2", lw=0.4)
        ax.grid(axis="y", color="#EDF0F2", lw=0.6)
        ax.tick_params(axis="y", length=0, pad=7)
    axes[-1].tick_params(labelleft=False) if len(axes) > 1 else None
    legend = [Line2D([], [], marker="o", color="#B84336", linestyle="--", label="Baseline"),
              Line2D([], [], marker="o", color="#2B6B9C", linestyle="none", label="Target precision reached")]
    if has_unmet:
        legend.append(Line2D([], [], marker="o", color="#2B6B9C", markerfacecolor="white", linestyle="none", label="Target precision unmet / provisional"))
    if has_uncertified:
        legend.append(Line2D([], [], marker="D", color="#B57A23", linestyle="none", label="Baseline-level distance uncertified"))
    if has_zero:
        legend.append(Line2D([], [], marker="<", color="#607080", markerfacecolor="white", linestyle="none", label="Zero-failure upper limit"))
    fig.suptitle("Weight-7 schedule comparison: same-shot joint X/Z failure", fontsize=13, y=0.985)
    fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, 0.087), ncol=min(3, len(legend)), fontsize=8, frameon=False)
    note = "All candidates: [[56,12,7]], k=12, R=7. Rate = 1-(1-Fjoint/N)^(1/7), per block (not per logical qubit). Nominal 1-sigma Wilson intervals; arrows denote upper limits. Candidate order is fixed; ratios and descriptive ranks are in CSV."
    if has_uncertified:
        note += " Schedules without a baseline-level distance certificate are exploratory."
    note += " Multiple-candidate selection requires independent validation; this scan does not prove a best schedule."
    fig.text(0.5, 0.014, textwrap.fill(note, 178), ha="center", va="bottom", fontsize=7.1)
    fig.tight_layout(rect=(0.02, max(0.16, 0.85 / height), 0.995, 0.985))
    fig.savefig(out / "comparison.pdf", bbox_inches="tight")
    fig.savefig(out / "comparison.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return {"plotted": True, "read_points": len(records), "plotted_points": len(selected), "candidates": len(ids), "out": str(out)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path("scan_results"))
    parser.add_argument("--out", type=Path, default=None, help="Default: RESULTS/plots")
    parser.add_argument("--include-incomplete", action="store_true", help="Include nonfinal snapshots with explicit provisional markers")
    parser.add_argument("--p", default="0.003,0.0025", help="Comma-separated panel probabilities")
    args = parser.parse_args()
    p_values = tuple(float(p) for p in args.p.split(",") if p.strip())
    if len(p_values) != 2 or any(not 0 < p < 1 for p in p_values):
        parser.error("Provide two valid physical probabilities, e.g. --p 0.003,0.0025")
    print(json.dumps(plot_scan(args.results, args.out, include_incomplete=args.include_incomplete, p_values=p_values), indent=2))


if __name__ == "__main__":
    main()
