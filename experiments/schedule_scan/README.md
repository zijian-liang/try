# Weight-7 schedule scan

Compare stabilizer-extraction schedules for the same `[[56,12,7]]` code under
fixed noise and decoder settings. The original eight-layer schedule is
`w7_baseline`; alternatives use eight, nine, or fourteen CNOT layers. Each
candidate passes a circuit-distance screen before Monte Carlo simulation.
X and Z failures come from the same physical shot and one joint decoding of
the full correlated detector error model (DEM).

Published figures are in [artifacts/schedule_scan](../../artifacts/schedule_scan/),
with compact plotting inputs in [data/schedule_scan](../../data/schedule_scan/).
These saved snapshots are separate from new runs in `scan_results/`.

## Environment and execution

Use Linux and Python 3.12 with the versions in `requirements.txt`. The scan
uses Unix file locks and process controls. From the repository root:

```bash
cd experiments/schedule_scan
bash install.sh
bash run_linux.sh --workers 8 --distance-workers 2
```

The installer creates a local `.venv`. Set `QEC_PYTHON` (or the legacy
`PYTHON_BIN`) to choose the installer interpreter. To reuse an environment
that already has the dependencies, skip installation and run:

```bash
QEC_PYTHON=/path/to/environment/bin/python bash run_linux.sh --workers 8
```

`run_linux.sh` prefers `QEC_PYTHON`, then `PYTHON_BIN`, then the environment in
`VENV_DIR` or the local `.venv`, and finally `python3`. Without `VENV_DIR`, it
also checks the former `.venv-qec312` directory before falling back to `python3`.
Choose worker counts for the available CPU and RAM. The CLI defaults request
380 simulation workers and 32 distance workers; resource checks may reduce
these counts, and each worker uses one numerical-library thread.

The full workflow generates candidates, screens every candidate, simulates
selected schedules at `p=0.003` then `p=0.0025`, and writes plots. To run stages
separately after installation:

```bash
.venv/bin/python main.py --stage generate --candidates 64
.venv/bin/python main.py --stage screen --workers 8 --distance-workers 2
.venv/bin/python main.py --stage simulate --workers 8 --max-circuits 12
```

For an inconclusive distance screen, increase the proof budget and retry:

```bash
.venv/bin/python main.py --stage screen --workers 8 --distance-workers 2 \
  --proof-seconds 1200 --retry-inconclusive
```

Reusing the same command and output directory resumes completed work. During
simulation, the first SIGINT finishes and saves the committed wave; a second
SIGINT terminates workers, and unfinished batches are replayed on resume.
Run only one writer per output directory. Worker count, batch size, and target
precision can change on resume. Physical settings, decoder identity, seeds,
or certification mode require compatible identities or a new `--output`.

## Experimental settings

| Setting | CLI default or fixed value |
|---|---|
| Code | `[[56,12,7]]`; 56 data and 56 syndrome qubits |
| Candidates | 64: 30 eight-layer (including baseline), 30 nine-layer, 4 fourteen-layer |
| Simulated schedules | Up to 12, including baseline; selected before simulation by schedule family |
| Circuit | 7 noisy rounds plus 2 ideal closing rounds |
| Noisy locations | 2,744 CNOTs, 392 preparations, 392 measurements; no idle noise |
| Noise | Two-qubit depolarizing probability `p` after CNOT; preparation/readout flips with probability `p` |
| Physical probabilities | `0.003,0.0025`; `--p` accepts either value or both |
| Decoder | Tesseract `0.1.1.dev20260910235247`, beam 2, `beam_climbing=False` |
| Distance budgets | 30 seconds heuristic upper-bound search and 300 seconds lower-bound proof per candidate |
| Stopping | At least 1,000 shots, at least one joint failure, relative binomial one-sigma error at most 20%; cap 100,000,000 shots |
| Batch size | 64; stopping is checked only after a complete committed wave |
| Seeds | Candidate and sampling seeds both 20260921; decoder seed 2384753 |

Other fixed decoder options are `num_det_orders=20`, `no_revisit_dets=True`,
`merge_errors=True`, `pqlimit=200000`, and `sparsify_errors=False`.
The CLI's current relative-error default is **0.20**; use
`--relative-error 0.10` or `0.05` for tighter precision. Saved snapshots record
their actual stopping settings, which can differ from current defaults.

## Distance gate

An upper bound of seven alone does not establish distance at least six.
The pipeline keeps `lower_bound`, `upper_bound`, and `exact` separate:

1. Check edge coverage, qubit conflicts, extraction order, and noiseless
   detectors and logical observables.
2. Build the full `R=7`, `idle=0`, `p=0.003` circuit DEM and remove duplicate
   labeled DEMs. This is not a complete test of graph isomorphism.
3. Search for a physical undetectable logical-fault witness. A witness of
   weight at most five rejects the candidate. Every accepted upper-bound
   witness is checked against allowed physical faults with zero detectors
   and nonzero logical response.
4. Prove a lower bound using complete CSS-projected exclusions through
   weight four, then weight five when the required premises hold. A timeout
   remains inconclusive.
5. By default, simulate only candidates with certified `lower_bound >= 6`
   and a physical upper-bound witness. `[6,6]` is exact six; `[6,7]` remains
   an interval. The baseline and at least one alternative must qualify.

Both supported positive noise probabilities have the same fault support,
so one distance screen applies to both. Historical seed orders and witnesses
are search aids; old lower bounds are not silently adopted. The
`--allow-uncertified` option, used in a new output directory, permits
exploratory candidates with upper bound at least six without a matching
lower-bound certificate. Such results cannot establish that all compared
circuits have distance at least six.

The candidate search is finite and translation invariant. It does not prove
a globally optimal schedule or guarantee finding circuit distance seven.

## Results and plotting

Each new run writes the following under `scan_results/`:

- `candidates.json`, `screening_summary.json`, and `selected_schedules.json`:
  generated schedules, screening outcomes, and the fixed comparison set.
- `distance/<ID>.json`: bounds, physical witnesses, proof status, and identity.
- `points/<ID>_p..._<hash>/`: circuit, DEM, configuration, and `result.json`.
- `plots/comparison.png`, `comparison.pdf`, `comparison.csv`, and `ranking.csv`:
  figures, counts, statistics, and descriptive rankings.

To plot a new run or reproduce the included figures without running simulations:

```bash
.venv/bin/python plot_scan.py --results scan_results
.venv/bin/python plot_scan.py --results ../../data/schedule_scan \
  --out ../../artifacts/schedule_scan
```

By default only finalized points appear. A point that hits `max_shots` without
meeting precision is marked accordingly. `--include-incomplete` also displays
clearly marked provisional snapshots. Plotting leaves input records unchanged.

All schedules use `k=12` and `R=7`. The plotted joint block failure rate per
noisy round is `1 - (1 - Fjoint/N)^(1/7)`, without a per-logical-qubit
normalization. Records retain X, Z, simultaneous, low-confidence, and decoder
failure counts, with `Fjoint = FX + FZ - Fboth + decodefail`. Difficult shots
are not discarded; invalid recoveries and decoder exceptions count as joint
failures.

Error bars transform the Wilson interval with `z=1`. They are nominal
one-sigma summaries, with no strict coverage guarantee under adaptive stopping
or multiple comparisons. Rankings and ratios are descriptive. Selecting the
lowest observed schedule introduces selection bias; confirm it against the
baseline with independent sampling in a new output directory, for example:

```bash
.venv/bin/python main.py --workers 8 --distance-workers 2 \
  --output confirm_results --sample-seed 20260922 \
  --only-ids parallel9_a2b5cacfbb84 --min-shots 100000 --max-shots 100000
```

Substitute a candidate ID from the generated manifest, and retain the original
`--seed` and `--candidates` settings. A lower observed failure rate measures
performance of the schedule with this decoder; it does not isolate fault
propagation, fault multiplicity, and decoder behavior as separate mechanisms.

## Validation and provenance

```bash
.venv/bin/python validate_scan.py
.venv/bin/python validate_scan.py --simulation-only \
  --output validation_simulation_resume.json
```

The default circuit validator covers five representative schedules, ideal
behavior, all 15 two-qubit Pauli outcomes at representative CNOT locations,
random multifault histories, and terminal Z-preparation faults. It compares
Stim detector/logical responses with independent classical Pauli propagation.
The simulation check exercises two-worker decoding and an 8-to-16-shot resume;
those diagnostic shots are not performance estimates. Use `--manifest` to
validate every candidate in an existing manifest.

The weight-7 checks use `HX=[A|B]`, `HZ=[B^T|A^T]` on `Z2 x Z14`, with
`A=1+y+x*y^10` and `B=1+y^3+x*y^11+x*y^2`; each check matrix has rank 22.
Ideal reference qubits track commuting encoded Bell observables for joint
X/Z bookkeeping and add no physical hardware or noisy locations.

Circuit boundaries and the joint Pauli-frame convention derive from IBM's
`BivariateBicycleCodes` implementation. The retained IBM source license is
[licenses/IBM_SOURCE_LICENSE.txt](../../licenses/IBM_SOURCE_LICENSE.txt).
This campaign uses zero idle noise and is not an exact replication of the
original paper's Figure 3. Runtime helpers retain `seed_schedules.json`,
`baseline_distance.json`, and `static_distance_results.json` as necessary
inputs; they should move with this experiment when it is integrated elsewhere.
