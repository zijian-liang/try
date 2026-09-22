# Local verification (2026-09-21)

The production defaults are unchanged: 380-worker ceiling, R=7, idle=0,
Tesseract beam=2, p=.003/.0025, min1000 shots, relative1sigma<=20%.
Local smoke-run caps were intentionally smaller and are recorded separately.

The first two ranked candidates remain [6,7]. The third-ranked candidate,
parallel8_4c4fb6c67ebc, is now exactly6: SAT supplied a real six-location witness
and fresh complete exclusion of weights<=5 supplied the lower bound. Original
baseline is exactly6. Full records and the new physical witness are included.

A total of290 randomly generated small decision partitions matched brute-force
enumeration; exact6/exact7 synthetic cases and timeout UNKNOWN also passed.
Same-budget distance resume scheduled no new tasks. The bundled actual witness
was rechecked at SAT budget0 and skipped all unnecessary SAT tasks.

Real same-shot Stim/Tesseract tests used384 shots at .003 and256 at .0025,
zero joint failures and zero decoder exceptions. These counts are INSUFFICIENT
for a performance claim, did not meet20% precision, and produce no slope.
Interrupted sampling resumed a committed wave without dropping or duplicating
shots. Another resume performed no sampling. Plots were rendered and checked.

Historical selection samples, local smoke samples, and the user's future
confirmation campaign are separate. Refer to validation/validation_summary.json
and README_zh.md for details. No candidate is claimed exact7 in this package.
