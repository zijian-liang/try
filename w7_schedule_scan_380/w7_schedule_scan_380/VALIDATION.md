# Validation report

All checks below passed for the supplied W7 schedule-scan package. The distance
section records completed proofs; the remaining sections are functional circuit
and execution checks. Diagnostic shots are not performance estimates.

## Independent circuit checks

The default validator selects actual candidates from the package generator:
the baseline, another eight-layer schedule, one nine-layer schedule, and both
fourteen-layer serial orientations. Every circuit uses seven noisy rounds,
zero idle noise, and two noiseless tail rounds.

| Candidate | CNOT layers | Fault histories | Result |
|---|---:|---:|---|
| `w7_baseline` | 8 | 288 | PASS |
| `parallel9_a2b5cacfbb84` | 9 | 288 | PASS |
| `parallel8_215be26722e3` | 8 | 288 | PASS |
| `serial14_x_then_z_08640a11162c` | 14 | 288 | PASS |
| `serial14_z_then_x_7e76e2e76408` | 14 | 288 | PASS |

For each candidate, validation includes:

- 64 no-fault shots with zero detectors and zero logical labels.
- Independent reconstruction of all 392 CNOT edges per round from the actual
  circuit, checking supports, orientations, and simultaneous-qubit conflicts.
- Exactly 2,744 noisy CNOT locations, 392 preparation faults, and 392 readout
  faults; no idle-noise locations; 504 detectors and 24 simultaneous logical labels.
- All fifteen two-qubit Pauli outcomes at fourteen representative first-round
  CNOT locations, covering all seven directions and both check orientations.
- 50 random histories of 1–12 physical faults at distinct locations.
- All 28 last-noisy-cycle Z-ancilla preparation faults, exercising the terminal boundary.
- Detector and logical responses compared against a separate classical X/Z
  Pauli-frame propagator. Logical labels are also checked directly using the
  logical-operator matrices. Histories exercising both logical components are included.

This totals 1,440 deterministic fault-history checks over the five schedules.
The machine-readable results and circuit hashes are in `validation_scan.json`.

## Actual two-worker decoding and resume

The real Tesseract decoder ran at p=0.003, R=7, idle=0, fixed beam=2:

1. Two workers completed exactly eight shots.
2. Increasing the minimum/maximum shots to sixteen, changing batch size from
   four to eight, and tightening relative-error settings from 0.10 to 0.05
   resumed the same experiment and completed exactly sixteen cumulative shots.
3. Repeating the completed request performed no further sampling or counting.

Experiment identity remained unchanged. Every completed wave was fully counted;
no pending batches or partial active wave remained. The same-shot identity
`Fjoint = FX + FZ - Fboth + decodefail` held. There were no decoder failures in
these sixteen diagnostic shots. This tiny count makes no logical-error-rate claim.
See `validation_simulation_resume.json` for counts, versions, and wave records.

Additional actual signal tests on the simulation module passed. A first SIGINT
after an eight-shot wave was committed drained all eight shots and returned
an interrupted status with no pending batches. A second-SIGINT test stopped
after four of eight shots; resuming with a changed worker count and even a
new shot cap below the pending wave completed the four pending shots first,
leaving eight total shots, no duplication and no pending batches.

## End-to-end command line

The actual main.py path was run with the baseline and a distinct nine-layer
candidate, using their freshly verified distance reports. Both p=0.003 and
p=0.0025 completed with two workers and eight diagnostic shots per point.
The program produced distance summaries, circuit/DEM files, count ledgers,
comparison plots and CSVs. Repeating the command reused both distance reports
and all four completed points without sampling further shots. The report
`validation_cli.json` contains the small-run identities and completion state.
These are functional checks, not new schedule-performance measurements.

## Reproduce

```bash
python validate_scan.py
python validate_scan.py --simulation-only --output validation_simulation_resume.json
```

For an existing candidate manifest, `python validate_scan.py --manifest candidates.json`
validates every listed candidate; this can take longer than the five-case default.

## Actual distance certificates

The first eight generated candidates were screened on their full R=7,
p=0.003, idle=0 circuits. All eight completed exact exclusion of physical
logical fault sets of weights up to five through the CSS projection method.
Every upper-bound witness was checked against supported physical faults.

| Candidate | Certified interval | Exact? |
|---|---|---|
| w7_baseline | [6,6] | Yes |
| parallel8_215be26722e3 | [6,7] | No |
| parallel8_4c4fb6c67ebc | [6,7] | No |
| parallel8_64ed9c3703ea | [6,7] | No |
| parallel9_12e4cb50666f | [6,7] | No |
| parallel9_29757d3627fd | [6,7] | No |
| parallel9_a2b5cacfbb84 | [6,7] | No |
| parallel9_ad6ecb11d763 | [6,7] | No |

The full records, source hashes and physical witnesses are included in
`certificate_examples.json`. These are eight validated candidates, not a claim
that all 64 defaults or every generated candidate pass the distance screen.
Fresh runs perform their own screening and do not silently import these records
as a cache. The baseline completed in 56.2 seconds; alternatives took
76.6–83.1 seconds each on the validation host with two concurrent workers.
These timings are observations, not server runtime guarantees.

A real repeated screening request reused a completed certificate and suppressed
a duplicate full DEM. Interrupting active screening left no living descendant
worker or heuristic-search processes.
