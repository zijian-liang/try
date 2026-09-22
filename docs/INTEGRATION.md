# Integrating these experiments

## Directory mapping

| Previous archive path | Current path |
| --- | --- |
| `w7_bb72_surface7_joint_tesseract_380/w7_joint_tesseract/` | `experiments/joint_memory/` |
| `w7_schedule_scan_380/w7_schedule_scan_380/` | `experiments/schedule_scan/` |
| `w7_schedule_scan_380/w7_schedule_scan_380/w7_high_slope_dc7_380/w7_high_slope_dc7_380/` | `experiments/high_slope/` |

The three experiment directories are complete local modules, not an installed Python package. Keep each directory intact when moving it. Launch its CLI from that directory, or invoke its script by path in a separate Python process. The programs use local imports such as `circuits`, `decoder`, and `simulation`; importing all three into one interpreter would create module-name collisions.

## What to preserve

- Keep the JSON inputs beside the corresponding code, including seed schedules, baseline distance evidence, and static distance records.
- Keep `experiments/high_slope/data/`: its historical log, measured candidates, and physical witnesses are program inputs. They are separate from the independent confirmation observations under the repository-level `data/high_slope/`.
- Move `data/` and `artifacts/` with the experiments if the larger project should retain these saved observations and figures. Update `scripts/reproduce_figures.py` if their relative layout changes.
- Preserve `licenses/IBM_SOURCE_LICENSE.txt` and the attribution in `PROVENANCE.md`.
- Give each new simulation its own output directory. Saved plotting snapshots omit generated circuits, detector-error models, and task state and are not complete restart checkpoints.

## Deliberate boundaries

Duplicated archive folders and repeated provenance/validation documents have been consolidated. Shared-looking Python modules remain local to each experiment: their byte hashes can be part of experiment identities and distance certificates. Moving or deduplicating those modules into a common package should be a separate, tested change with fresh validation of any affected evidence.

`generated/` is ignored by Git and is the default destination for regenerated figures. The supplied images in `artifacts/` and observations in `data/` are tracked explicitly. No machine-specific execution path is required by the run scripts; historical path strings in saved records are provenance only.

## Checks after moving

```bash
python scripts/reproduce_figures.py
python -m unittest discover -s experiments/high_slope -p 'test_plot_results.py'
```

Use the circuit and distance validation commands in each experiment's README before starting a new simulation. Plot reproduction checks saved data handling, not the physical correctness of a new simulation or the validity of historical distance proofs.
