# Saved observations and plotting inputs

These snapshots accompany the supplied figures. They preserve the original JSON and CSV bytes, including counts, identities, settings, and saved audit records. Absolute paths inside historical records describe the original runs; the plotting commands locate inputs through the current directory layout.

| Dataset | Point snapshots | Contents |
| --- | ---: | --- |
| `joint_memory/` | 27 | Three codes at nine probabilities; point results, result index, and plot summary |
| `schedule_scan/` | 24 | Twelve schedules at two probabilities; selection and campaign manifests, the twelve corresponding distance snapshots, and comparison/ranking tables |
| `high_slope/` | 8 | Four independently evaluated schedules at two probabilities; point results, distance records, selection/configuration, and summary/slopes |

Each point is stored at `points/<original-point-id>/result.json`. Generated `model.dem`, `circuit.stim`, locks, solver tasks, and virtual environments are omitted. These are plotting and inspection snapshots, not full restart checkpoints.

The [manifest](manifest.json) maps each copied data or figure asset to its source path relative to the original `weight_7_BB_code/` directory and records its SHA256. The mapping covers the imported assets only; later generated output has its own destination under `generated/`.

## Reproduction and interpretation

Run `python scripts/reproduce_figures.py` from the repository root after installing `requirements-plot.txt`. Joint-memory and scan plots read point snapshots; high-slope plot-only mode reads the supplied `simulation_summary.json`. That mode reproduces saved labels and intervals and does not recertify circuit distances.

The high-slope experiment's historical candidate-selection inputs are stored separately in `experiments/high_slope/data/`. Do not combine their counts with the independent confirmation data here. Ranking and slopes are descriptive. Adaptive stopping and candidate selection require care when interpreting nominal error bars.
