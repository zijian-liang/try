# W7 high-slope confirmation

Select W7 schedules using completed historical measurements, audit their circuit
distance, and collect an independent confirmation sample. The physical model is
`[[56,12,7]]`, seven noisy rounds and two ideal closing rounds, zero idle noise,
same-shot joint X/Z failures, the full correlated detector-error model (DEM),
and Tesseract with beam 2. Preparation and measurement fail with probability
`p`; each CNOT is followed by two-qubit depolarization with total probability `p`.

## Selection and interpretation

The default selection contains the three largest historical two-point slopes
among previously unresolved `[6,7]` candidates, plus the original W7 baseline.
Only final `point_done` records in `data/old_scan.log` enter the ranking.

| Candidate | Historical slope, nominal 1 sigma | Previously established distance status |
| --- | ---: | --- |
| `parallel9_a2b5cacfbb84` | 6.6687 +/- 1.3864 | `[6,7]`, unresolved |
| `parallel9_4ec1ca467df4` | 6.0208 +/- 1.5513 | `[6,7]`, unresolved |
| `parallel8_4c4fb6c67ebc` | 5.8322 +/- 1.2733 | 6, physical witness and exclusion through weight 5 |
| `w7_baseline` | 5.5930 +/- 1.3025 | 6, baseline control |

The third candidate remains in confirmation as an additional control. Historical
slopes are computed at `p = 0.003` and `0.0025` using

```text
q = 1 - (1 - Fjoint/N) ** (1 / (12 * 7))
alpha = log(q_0.003 / q_0.0025) / log(0.003 / 0.0025)
```

Here `q` is an equivalent normalization of the joint block failure probability;
it is not a measured failure probability for an individual logical qubit.
Selecting large noisy slopes creates selection bias. The listed uncertainties
do not adjust for candidate selection or adaptive stopping, and neither the
ranking nor a steep slope establishes distance 7 or superiority to the baseline.
Confirmation uses a separate seed and never pools historical selection counts.

## Run a new campaign

Use Linux or WSL: the coordinator uses `fcntl` to prevent concurrent writers.
Use the existing Python environment that successfully runs the pinned Tesseract
version, or install the complete requirements in a new environment.

```bash
cd experiments/high_slope
python3 -m pip install -r requirements.txt
python3 -u main.py --output results
```

Defaults are in [settings.py](settings.py): three candidates plus baseline,
`p = [0.003, 0.0025]`, seed `20260923`, a minimum of 1,000 shots, maximum of
100,000,000 shots, and a 20% nominal relative standard-error target. Simulation
and SAT worker ceilings are 380, automatically reduced for CPU and memory
availability; screening has a separate ceiling of 32. Use smaller `--workers`,
`--distance-workers`, and `--screening-workers` limits as needed.

```bash
# Select candidates without distance solving or sampling.
python3 main.py --stage select --output results

# Audit distance; retry unresolved tasks with a larger per-task budget.
python3 -u main.py --stage distance --output results
python3 -u main.py --stage distance --output results --sat-seconds 3600 --retry-unknown

# Collect new independent samples using the current distance records.
python3 -u main.py --stage simulate --output results

# Rebuild summaries and figures from the raw point records.
python3 main.py --stage plot --output results
```

Repeating a command resumes completed work. Use a new output directory when
changing the seed or candidate set, and run only one coordinator per directory.
The first interrupt allows the committed sampling wave to be saved. Every wave
is included in full before testing the stopping rule. Decoder exceptions are
recorded separately and counted as joint failures. A shot limit does not imply
that the precision target has been reached.

## Distance certificates and limits

Historical log labels are not reused as proof. Each audit rebuilds the complete
`R=7`, idle-free circuit and DEM, verifies small-weight lower bounds, and checks
upper-bound witnesses against distinct physical fault locations. The bundled
`data/known_physical_witnesses.json` contains the baseline and third-candidate
witnesses, which are checked again before use.

The exact stage partitions the at-most-six-fault SAT decision problem over both
CSS projections, 12 logical bits per projection, and minimum selected response
indices. The default 16 shards per logical bit give at most 384 tasks per circuit.
These projections relax the complete correlated fault problem; they do not split
the Monte Carlo decoder into independent CSS decoders.

- Excluding all faults of weight at most six in every covering shard, together
  with a verified seven-location witness, establishes circuit distance 7.
- A verified six-location witness and complete exclusion through weight five
  establish circuit distance 6.
- A timeout, error, or incomplete shard remains unresolved. A projected SAT
  witness is an upper bound only after mapping to the full physical fault model.

The default 300 seconds per SAT task does not guarantee a conclusive answer.
Solver UNSAT answers are used as exact decisions; this package does not include
an independent DRAT proof checker. No bundled candidate is claimed to have
certified circuit distance 7. Confirmation can continue while distance remains
unresolved, with the bounds displayed explicitly.

## Figures and compact data

Published figures are in [artifacts/high_slope](../../artifacts/high_slope/), with
their saved summary and raw point records in [data/high_slope](../../data/high_slope/).
To render the saved summary without simulation dependencies or another audit:

```bash
# Run from this experiment directory; this mode requires matplotlib only.
python3 plot_results.py --summary ../../data/high_slope/simulation_summary.json --output ../../artifacts/high_slope
```

This mode preserves the stored statistics and distance labels; it does not
recertify them. `--summary` and `--distance-records` are mutually exclusive.
For raw records, `python3 plot_results.py --output results` keeps circuit-hash
validation and reconstructs the summary before plotting.

The raw-record and full-campaign workflows produce `simulation_summary.json/.csv`,
`two_point_slopes.csv`, and `independent_confirmation.png/.pdf`, plus
`figure_caption.txt`. The `--summary` mode only writes the PNG, PDF, and caption;
it does not rewrite the summary or slope tables. The plot uses transformed nominal Wilson
one-sigma intervals, open markers for unmet precision targets, and upper-limit
arrows for zero-failure points. Only completed points are plotted. These
intervals are not anytime-valid confidence sequences under adaptive stopping.

Historical selection inputs stay in this experiment's `data/` directory. The
launcher notice on the first line of `old_scan.log` has been translated to
English; all parsed counts, bounds, line numbers, rankings, and schedules are
unchanged. This translation changes the file hash. The original log hash remains
in `data/measured_candidates.json` as provenance, not as the translated file's
checksum.

The W7 supports and inherited IBM circuit conventions are documented in the
repository provenance notes. Preserve [IBM's source license](../../licenses/IBM_SOURCE_LICENSE.txt)
when integrating the applicable upstream-derived material into another project.
