# Joint X/Z circuit-level memory experiment

Compare the weight-7 `[[56,12,7]]`, IBM BB `[[72,12,6]]`, and rotated surface
`[[49,1,7]]` codes using full noisy-circuit sampling in Stim and joint detector
error model (DEM) decoding in Tesseract. Each sampled shot contains both X and Z
error components. The primary outcome is any logical error or decoder failure
in that shot.

## Install and run

Run these commands from this directory on Linux:

```bash
bash install.sh
source .venv/bin/activate
python main.py --smoke
```

The smoke run uses 16 shots per code at `p=0.003` and at most two workers. Its
output is a functional check, not a logical error-rate estimate. For a campaign,
choose a worker count appropriate for the host:

```bash
bash run_linux.sh --workers 8
```

Both shell scripts select Python in this order: `QEC_PYTHON`, `PYTHON_BIN`,
`bin/python` in `VENV_DIR` (or the local `.venv`), then `python3`.
`install.sh` creates or reuses `VENV_DIR` (default `.venv`) and installs there.
An explicit interpreter lets the runner use a shared environment when this
experiment is integrated into another project:

```bash
QEC_PYTHON=/path/to/shared-env/bin/python bash run_linux.sh --workers 8
```

Use the same `VENV_DIR` for installation and execution when choosing a custom
environment location. Relative environment paths are resolved from this
experiment directory.

The configured default is 376 workers for the original large server. Each
worker uses one CPU thread; the runner caps concurrency using CPU affinity and
measured decoder memory use. A large Tesseract search queue can make individual
shots slow. The queue limit is not a time limit.

Select codes and error probabilities, or resume with a tighter stopping target:

```bash
python main.py --codes w7 --p 0.006,0.004,0.002 --workers 8
python main.py --workers 8 --relative-error 0.05 --max-shots 300000000
```

Use `--output results_custom` to select a separate campaign directory. Completed
batches are reused only when the saved experiment identity matches. Worker
count, batch size, shot limit, and precision target may change when resuming;
changes to the circuit, noise, rounds, or decoder define a different experiment.
The output-directory lock prevents simultaneous writers. The first Ctrl+C
stops new work and saves running batches; a second stops workers, leaving
unfinished batches to be recalculated on the next run.

## Experiment settings

Defaults are in [campaign_config.py](campaign_config.py); dependency constraints
and the pinned Stim/Tesseract versions are in [requirements.txt](requirements.txt).

| Setting | Default |
| --- | --- |
| Codes | `w7,bb72,surface7` |
| Noisy rounds `R` | Static distance `d`: 7, 6, and 7, respectively |
| Closing boundary | Two additional noiseless syndrome cycles |
| Physical error probabilities | `0.006, 0.005, 0.004, 0.0035, 0.003, 0.0025, 0.002, 0.0015, 0.001` |
| Two-qubit gates | After each ideal CNOT, probability `p` of a uniformly selected nonidentity two-qubit Pauli |
| Preparation / measurement | Opposite-eigenstate preparation / readout flip, each with probability `p` |
| Idle noise | `IDLE_SCALE=0`; otherwise data-qubit depolarization with probability `IDLE_SCALE*p` |
| Decoder | Joint undecomposed X/Z DEM; `det_beam=2`, `beam_climbing=False`, `num_det_orders=20`, `pqlimit=200000` |
| Minimum / maximum shots per point | 1,000 / 100,000,000 |
| Precision target | At least one failure and relative estimated binomial standard error `sqrt(P*(1-P)/N)/P <= 0.10` |
| Randomness | Fixed campaign seed and independent batch seeds recorded for resumption |

X/Z preparation and readout use native `RX/MX` and `R/M` operations, including
for the surface code; no additional noisy Hadamard gates are inserted. Idle
positions follow the IBM data-qubit timing convention. IBM's original noise
model uses idle probability `p`, whereas this campaign defaults to zero idle
noise and uses Tesseract instead of separate X/Z BP+OSD decoding. These results
are therefore not an exact reproduction of the IBM paper's simulation.

BB72 has `ell=m=6`, `A=x^3+y+y^2`, and `B=y^3+x+x^2`, with CNOT orders
`sX=['idle',1,4,3,5,0,2]` and `sZ=[3,5,0,1,2,4,'idle']`. Z ancillas are prepared
at the end of each cycle for the next one. Their preparation before the first
noisy cycle is ideal; preparation errors at the end of the last noisy cycle
remain present in the clean closing cycles. Detectors compare consecutive
syndromes, including the initial and terminal boundaries.

The weight-7 code uses `x^2=y^14=1`, `A=1+y+x*y^10`,
`B=1+y^3+x*y^11+x*y^2`, `HX=[A|B]`, and `HZ=[B^T|A^T]`. Its retained schedule has
eight CNOT layers. Its construction and current schedule are defined in
[circuits.py](circuits.py).

## Logical outcomes and normalization

Ideal encoded Bell references make all logical X/Z Pauli-frame components
simultaneously observable through commuting extended observables. The
references are noiseless simulation bookkeeping and do not count as hardware.
For residual `E=X(e_X)Z(e_Z)`, the logical X component is measured by `L_Z e_X`
and the logical Z component by `L_X e_Z`. Either component may fail, or both may
fail in the same shot.

| Count | Meaning |
| --- | --- |
| `N` | Sampled shots, not sampled logical qubits |
| `FX`, `FZ` | Shots with a nonzero logical X or Z residual after a valid recovery |
| `Fboth` | Shots with both residual components nonzero |
| `decodefail` | Shots without a valid decoder recovery |
| `Fjoint` | `FX + FZ - Fboth + decodefail` |
| `lowconfidence` | Decoder low-confidence flags; these shots remain included |
| `logical_X`, `logical_Z` | Known residual-error counts for individual logical qubits |

The decoder is heuristic and does not guarantee optimal recovery. Invalid
recoveries and decoder exceptions contribute to joint operational failure, but
are not assigned invented X/Z residuals. Sector plots conservatively count
`FX+decodefail` or `FZ+decodefail` and are labeled operational upper bounds when
decoder failures occur.

For `P=Fjoint/N`, the plotted per-noisy-round rate follows the convention
`p_round=1-(1-P)^(1/R)`. This is a normalization of final block failure, not an
exact inference of an independent logical transition process.

| Code | Data qubits | Syndrome ancillas | Hardware total | Logical qubits |
| --- | ---: | ---: | ---: | ---: |
| Weight-7 | 56 | 56 | 112 | 12 |
| BB72 | 72 | 72 | 144 | 12 |
| Surface, one patch | 49 | 48 | 97 | 1 |
| Surface, 12 patches | 588 | 576 | 1,164 | 12 |

The equal-`k=12` comparison analytically converts the measured single-surface-
patch probability to `P_12=1-(1-P_surface)^12`, assuming independent patches.
The resulting per-round rate is `1-(1-P_surface)^(12/R)`. Shot counts remain
those of the single sampled patch; within-patch X/Z correlations are preserved.

Error bars transform the endpoints of Wilson intervals with `z=1` (nominal
68.27% coverage). Zero observed failures produce an upper-limit marker, not a
known zero rate. Points that reach the shot cap before the precision target
remain available with their stopping status. The estimated relative standard
error controls adaptive stopping; the plotted intervals are not
sequentially valid confidence sequences. No low-error power-law fit is inferred
from zero-failure samples.

## Retained plots and regeneration

Saved figures are in [artifacts/joint_memory](../../artifacts/joint_memory/), and
compact aggregate plotting inputs are in [data/joint_memory](../../data/joint_memory/).
From this directory, regenerate the figures without running Monte Carlo:

```bash
python plot_results.py --results ../../data/joint_memory --out ../../artifacts/joint_memory
```

For new campaign output, use `python plot_results.py --results results`. The
plotter reads `result.json` files recursively, or accepts a single result file.
It rejects duplicate points and incompatible noise/decoder settings instead
of silently merging them. Compact committed plot inputs are for analysis and
plot regeneration; resume simulations from a complete local output directory
containing its batch records.

Outputs include PDF/PNG pairs for `joint_one_block_1sigma`,
`joint_one_block_points`, `joint_equal_k12_1sigma`, `joint_equal_k12_points`, and
`sectors_X_then_Z_1sigma`, plus `plot_summary.csv`. For derived 12-patch CSV rows,
only probabilities and intervals are transformed; raw counts still describe
the single patch.

## Distance evidence and validation

The retained [static_distance_results.json](static_distance_results.json) and
[distance_bounds.json](distance_bounds.json) report these certified distances
for the saved schedules, default rounds, zero idle noise, and two clean closing
cycles:

| Code | Static distance | CNOT layers / cycle | Time steps / cycle | CNOTs / cycle | Circuit distance |
| --- | ---: | ---: | ---: | ---: | ---: |
| Weight-7 | 7 | 8 | 10 | 392 | 6 |
| BB72 | 6 | 7 | 8 | 432 | 6 |
| Surface | 7 | 4 | 6 | 168 | 7 |

The weight-7 and BB72 lower bounds use complete low-weight fault exclusion;
explicit six-fault logical witnesses establish the matching upper bounds.
The surface lower bound uses exact graph search on complete CSS projections,
with an explicit seven-fault witness for the upper bound. Projection is a
relaxation used in the distance proof, not a replacement of joint Monte Carlo
sampling by independent X/Z experiments. Neither static distance nor a fitted
Monte Carlo slope proves circuit distance. The weight-7 schedule was selected
from a finite search and is not claimed to be globally optimal.

Distance labels are drawn from `distance_bounds.json` and apply only to matching
rounds, idle noise, and schedules. A fault witness alone proves an upper bound;
an equality requires a matching proved lower bound. The retained BB72 equality
is this repository's audit result, not a claim that IBM's published upper bound
alone established equality.

Optional checks, without production Monte Carlo:

```bash
python validate_circuits.py --output validation_circuits.json
python static_distance.py --output static_distance_recheck.json
python distance_audit.py --codes w7 bb72 surface7 --seconds 0 --projected-five-seconds 120 --output distance_bounds_recheck.json
```

Circuit validation compares Stim with independent classical Pauli propagation.
Its optional `--ibm-source /path/to/BivariateBicycleCodes` checks the original
IBM operation ordering when that source is available separately. Exact distance
audits can require about 1 GB of memory and substantial runtime. Exhausting a
search budget is not a proof; inspect the reported bounds and completion flags.

## Provenance

BB72 construction, scheduling, and closing-boundary conventions follow IBM's
[BivariateBicycleCodes](https://github.com/sbravyi/BivariateBicycleCodes).
The weight-7 polynomial supports were transcribed from the originally supplied
code-construction reference and are recorded explicitly above and in the source.
Sampling uses [Stim](https://github.com/quantumlib/Stim); joint decoding uses
[Tesseract](https://github.com/quantumlib/tesseract-decoder).

Reference: S. Bravyi et al., *High-threshold and low-overhead fault-tolerant quantum
memory*, Nature 627, 778-782 (2024),
[doi:10.1038/s41586-024-07107-7](https://doi.org/10.1038/s41586-024-07107-7).
The retained IBM source license is [IBM_SOURCE_LICENSE.txt](../../licenses/IBM_SOURCE_LICENSE.txt).
