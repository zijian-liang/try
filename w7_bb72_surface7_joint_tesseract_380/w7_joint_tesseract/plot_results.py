#!/usr/bin/env python3
"""Plot saved aggregates only; never touches a running simulation's state."""
from __future__ import annotations

import argparse
import csv
import json
import math
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ORDER = ("w7", "bb72", "surface7")
COLORS = {"w7": "#B94135", "bb72": "#2871A6", "surface7": "#428663"}
MARKERS = {"w7": "o", "bb72": "s", "surface7": "^"}
CODE_LABELS = {
    "w7": r"Weight-7 $[[56,12,7]]$",
    "bb72": r"BB $[[72,12,6]]$",
    "surface7": r"Surface $[[49,1,7]]$",
}


def wilson(successes: int, trials: int) -> tuple[float, float]:
    """Wilson score interval with z=1 (nominal 68.27% confidence)."""
    if trials <= 0 or not 0 <= successes <= trials:
        raise ValueError("Wilson counts must satisfy 0 <= failures <= shots, shots > 0")
    q = successes / trials
    den = 1 + 1 / trials
    centre = (q + 1 / (2 * trials)) / den
    half = math.sqrt(q * (1 - q) / trials + 1 / (4 * trials**2)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def round_rate(shot_probability: float, rounds: int, patches: int = 1) -> float:
    """IBM normalization, after independent-patch block conversion if requested."""
    if rounds <= 0 or patches <= 0:
        raise ValueError("rounds and patches must be positive")
    if shot_probability >= 1:
        return 1.0
    return -math.expm1((patches / rounds) * math.log1p(-shot_probability))


def load_records(results: Path) -> list[dict]:
    paths = [results] if results.is_file() else sorted(results.rglob("result.json"))
    records = []
    seen = {}
    physical_settings = {}
    schedules = {}
    for path in paths:
        record = json.loads(path.read_text(encoding="utf-8"))
        code = record.get("code_name")
        if code not in ORDER:
            continue
        counts = record["counts"]
        n = int(counts["N"])
        if n == 0:
            continue
        for key in ("Fjoint", "FX", "FZ", "Fboth", "decodefail"):
            value = int(counts.get(key, 0))
            if not 0 <= value <= n:
                raise ValueError(f"Invalid {key} count in {path}")
        fx, fz, both, failed = (int(counts.get(k, 0)) for k in ("FX", "FZ", "Fboth", "decodefail"))
        if both > min(fx, fz) or int(counts["Fjoint"]) != fx + fz - both + failed:
            raise ValueError(f"Same-shot count identity failed in {path}: expected Fjoint=FX+FZ-Fboth+decodefail")
        key = (code, float(record["p"]), int(record["rounds"]))
        if key in seen:
            raise ValueError(f"Multiple experiment records for {key}:\n  {seen[key]}\n  {path}\nSelect one campaign directory; incompatible histories must not be merged.")
        seen[key] = path
        identity = record.get("identity", {})
        settings = {key: identity.get(key) for key in ("protocol", "idle_scale", "decoder", "decoder_settings",
                    "decoder_adapter_sha256", "versions", "observables", "count_rule")}
        settings_key = json.dumps(settings, sort_keys=True)
        if physical_settings and settings_key not in physical_settings:
            raise ValueError(f"Mixed noise/decoder protocols in {results}; select one campaign directory.")
        physical_settings[settings_key] = path
        schedule = record.get("metadata", {}).get("schedule")
        schedule_key = json.dumps(schedule, sort_keys=True)
        if code in schedules and schedules[code] != schedule_key:
            raise ValueError(f"Mixed schedules for {code}; select one campaign directory.")
        schedules[code] = schedule_key
        record["_path"] = str(path)
        records.append(record)
    if not records:
        raise ValueError(f"No nonempty result.json records found in {results}")
    return sorted(records, key=lambda r: (ORDER.index(r["code_name"]), int(r["rounds"]), float(r["p"])))


def applicable_distance(code: str, records: list[dict], bounds: dict) -> dict:
    entry = bounds.get("codes", {}).get(code, {})
    if not entry:
        return {}
    expected_rounds = entry.get("rounds", bounds.get("rounds"))
    expected_idle = entry.get("idle_scale", entry.get("noise", {}).get("idle_scale", bounds.get("idle_scale")))
    expected_schedule = entry.get("schedule")
    if expected_rounds is None or expected_idle is None or expected_schedule is None:
        return {}
    for record in records:
        if (record.get("distance_bounds") or {}).get("status") == "not_applicable":
            return {}
        idle = record.get("identity", {}).get("idle_scale", record.get("metadata", {}).get("idle_scale"))
        schedule = record.get("metadata", {}).get("schedule")
        if int(record["rounds"]) != int(expected_rounds) or idle != expected_idle or schedule != expected_schedule:
            return {}
    return entry


def distance_text(entry: dict) -> str:
    lo, hi = entry.get("lower_bound"), entry.get("upper_bound")
    exact = entry.get("exact", False)
    if exact and hi is not None and lo == hi:
        return rf"$d_{{\rm circ}}={int(hi)}$"
    if lo is not None and hi is not None:
        return rf"${int(lo)}\leq d_{{\rm circ}}\leq {int(hi)}$"
    if hi is not None:
        return rf"$d_{{\rm circ}}\leq {int(hi)}$"
    return r"$d_{\rm circ}$: unverified"


def series_label(code: str, entry: dict, patches: int, rounds: int, multi_rounds: bool) -> str:
    name = CODE_LABELS[code]
    if patches == 12:
        name = r"12 independent surface patches ($k=12$)"
    label = name + ", " + distance_text(entry)
    return label + (rf", $R={rounds}$" if multi_rounds else "")


def transformed_point(record: dict, component: str = "joint", patches: int = 1) -> tuple[float, float, float, int]:
    counts = record["counts"]
    n = int(counts["N"])
    if component == "joint":
        f = int(counts["Fjoint"])
    elif component in ("X", "Z"):
        # Invalid/aborted decoding has no certified logical residual. Count it
        # as operational sector failure, without inventing an X/Z residual.
        f = int(counts["F" + component]) + int(counts.get("decodefail", 0))
    else:
        raise ValueError(component)
    lo, hi = wilson(f, n)
    r = int(record["rounds"])
    return round_rate(f / n, r, patches), round_rate(lo, r, patches), round_rate(hi, r, patches), f


def draw_series(ax, records: list[dict], code: str, bounds: dict, component: str,
                patches: int, errors: bool, multi_rounds: bool) -> None:
    rounds_values = sorted({int(r["rounds"]) for r in records if r["code_name"] == code})
    for rounds in rounds_values:
        subset = [r for r in records if r["code_name"] == code and int(r["rounds"]) == rounds]
        subset.sort(key=lambda r: float(r["p"]))
        pts = [(float(r["p"]), *transformed_point(r, component, patches)) for r in subset]
        if not pts:
            continue
        valid = [(p, y, lo, hi) for p, y, lo, hi, f in pts if f > 0]
        label = series_label(code, applicable_distance(code, subset, bounds), patches, rounds, multi_rounds)
        color, marker = COLORS[code], MARKERS[code]
        if valid:
            x, y, lo, hi = zip(*valid)
            yerr = [[max(0.0, yy - ll) for yy, ll in zip(y, lo)], [max(0.0, hh - yy) for yy, hh in zip(y, hi)]] if errors else None
            ax.errorbar(x, y, yerr=yerr, color=color, marker=marker, label=label,
                        lw=1.6, ms=5.5, capsize=2.7, elinewidth=1, markeredgewidth=0.8)
        else:
            ax.plot([], [], color=color, marker=marker, label=label)
        for p, y, lo, hi, f in pts:
            if f == 0:
                # Zero-failure samples establish an upper limit, not y=0 on log axes.
                ax.plot(p, hi, marker="v", markersize=6, color=color, markerfacecolor="white", linestyle="none")
                ax.plot(p, hi / 1.7, alpha=0, marker=".", linestyle="none")
                ax.annotate("", xy=(p, hi / 1.65), xytext=(p, hi),
                            arrowprops={"arrowstyle": "->", "color": color, "lw": 1.1})


def decorate(ax, title: str) -> None:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Physical error probability $p$")
    ax.set_ylabel(r"Block failure rate per noisy round")
    ax.set_title(title, fontsize=11)
    ax.grid(which="major", color="#D9DDE1", lw=0.65)
    ax.grid(which="minor", color="#ECEFF1", lw=0.35)
    ax.legend(loc="best", fontsize=8, framealpha=0.96)


def save_figure(fig, out: Path, stem: str) -> None:
    fig.savefig(out / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(out / f"{stem}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_table(records: list[dict], out: Path) -> None:
    fields = ["code", "patches", "p", "rounds", "k", "N", "Fjoint", "FX", "FZ", "Fboth", "decodefail",
              "lowconfidence", "P_fail_shot", "p_round", "p_round_low_1sigma", "p_round_high_1sigma",
              "X_operational_p_round", "Z_operational_p_round", "zero_failures", "experiment_id", "source"]
    with (out / "plot_summary.csv").open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()
        for r in records:
            c = r["counts"]
            for patches in ([1, 12] if r["code_name"] == "surface7" else [1]):
                y, lo, hi, f = transformed_point(r, "joint", patches)
                shot_q = -math.expm1(patches * math.log1p(-int(c["Fjoint"]) / int(c["N"]))) if int(c["Fjoint"]) < int(c["N"]) else 1.0
                row = {"code": r["code_name"], "patches": patches, "p": r["p"], "rounds": r["rounds"],
                       "k": int(r.get("k", 1)) * patches, "P_fail_shot": shot_q,
                       "p_round": y, "p_round_low_1sigma": lo, "p_round_high_1sigma": hi,
                       "X_operational_p_round": transformed_point(r, "X", patches)[0],
                       "Z_operational_p_round": transformed_point(r, "Z", patches)[0],
                       "zero_failures": f == 0, "experiment_id": r.get("experiment_id", ""), "source": r["_path"]}
                # Counts always refer to sampled single-code shots. For patches=12,
                # probabilities are analytically derived, never fabricated counts.
                row.update({key: c.get(key, 0) for key in ("N", "Fjoint", "FX", "FZ", "Fboth", "decodefail", "lowconfidence")})
                writer.writerow(row)


def make_plots(records: list[dict], out: Path, bounds: dict) -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42})
    multi_rounds = len({int(r["rounds"]) for r in records}) > 1
    rounds_note = "" if multi_rounds else f"; R={records[0]['rounds']} noisy rounds"
    has_zero = any(int(r["counts"]["Fjoint"]) == 0 for r in records)
    has_decodefail = any(int(r["counts"].get("decodefail", 0)) > 0 for r in records)
    base_note = "1-sigma Wilson intervals, transformed monotonically" + rounds_note + "."
    if has_zero:
        base_note += " Open downward markers: zero-failure upper limits."
    if has_decodefail:
        base_note += " Decoder failures included; see CSV."
    for equal_k in (False, True):
        for errors in (True, False):
            fig, ax = plt.subplots(figsize=(7.8, 5.4))
            for code in ORDER:
                draw_series(ax, records, code, bounds, "joint", 12 if equal_k and code == "surface7" else 1, errors, multi_rounds)
            title = "Joint X/Z block failure: equal k=12" if equal_k else "Joint X/Z block failure: one code block"
            decorate(ax, title)
            note = base_note if errors else "Point estimates" + rounds_note + "." + (" Open downward markers: zero-failure upper limits." if has_zero else "")
            fig.text(0.5, 0.01, textwrap.fill(note, 112), ha="center", va="bottom", fontsize=7.2)
            fig.tight_layout(rect=(0, 0.085, 1, 1))
            stem = ("joint_equal_k12" if equal_k else "joint_one_block") + ("_1sigma" if errors else "_points")
            save_figure(fig, out, stem)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    for ax, component in zip(axes, ("X", "Z")):
        for code in ORDER:
            draw_series(ax, records, code, bounds, component, 1, True, multi_rounds)
        decorate(ax, f"Logical {component} component: one code block" + ("\n(operational upper bound)" if has_decodefail else ""))
        ax.set_ylabel("Sector failure rate per noisy round")
    sector_note = "Same noisy shots; X then Z. " + base_note
    if has_decodefail:
        sector_note += " Sector curves conservatively include every decoder failure."
    fig.text(0.5, 0.01, textwrap.fill(sector_note, 168), ha="center", va="bottom", fontsize=7.2)
    fig.tight_layout(rect=(0, 0.11, 1, 1))
    save_figure(fig, out, "sectors_X_then_Z_1sigma")
    write_table(records, out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path("results"), help="Results root or a single result.json")
    parser.add_argument("--out", type=Path, default=None, help="Output directory; default RESULTS/plots")
    parser.add_argument("--distance-audit", type=Path, default=Path(__file__).with_name("distance_bounds.json"))
    args = parser.parse_args()
    bounds = json.loads(args.distance_audit.read_text(encoding="utf-8")) if args.distance_audit.exists() else {}
    records = load_records(args.results)
    out = args.out or ((args.results.parent if args.results.is_file() else args.results) / "plots")
    out.mkdir(parents=True, exist_ok=True)
    make_plots(records, out, bounds)
    print(f"Read {len(records)} points. Wrote 5 PDF/PNG figures and plot_summary.csv to {out}")


if __name__ == "__main__":
    main()
