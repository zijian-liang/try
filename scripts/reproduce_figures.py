#!/usr/bin/env python3
"""Replot saved observations without launching or modifying a simulation."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS = ("joint_memory", "schedule_scan", "high_slope")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment", choices=("all", *EXPERIMENTS), default="all")
    parser.add_argument(
        "--output", type=Path, default=ROOT / "generated" / "figures",
        help="Destination root; one subdirectory is created per experiment.",
    )
    args = parser.parse_args()
    selected = EXPERIMENTS if args.experiment == "all" else (args.experiment,)
    destination = args.output.resolve()
    # Keep the saved observations immutable even if an output path is mistyped.
    data_root = (ROOT / "data").resolve()
    if destination == data_root or data_root in destination.parents or destination in data_root.parents:
        parser.error("Choose an output directory separate from the saved data tree.")
    for name in selected:
        source = ROOT / "experiments" / name
        data = ROOT / "data" / name
        output = destination / name
        output.mkdir(parents=True, exist_ok=True)
        if name == "high_slope":
            command = [sys.executable, str(source / "plot_results.py"),
                       "--summary", str(data / "simulation_summary.json"),
                       "--output", str(output)]
        else:
            script = "plot_results.py" if name == "joint_memory" else "plot_scan.py"
            command = [sys.executable, str(source / script),
                       "--results", str(data), "--out", str(output)]
        print(f"Reproducing {name} -> {output}", flush=True)
        subprocess.run(command, cwd=source, check=True)
    print(f"Finished. Figures are in {destination}")


if __name__ == "__main__":
    main()
