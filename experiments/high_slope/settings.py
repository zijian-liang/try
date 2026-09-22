"""Edit these defaults, then run: python main.py.

This is a NEW independent confirmation campaign. It never writes old results.
"""
from pathlib import Path

TOP_CANDIDATES = 3  # Largest historical two-point slopes; plus baseline.
OUTPUT_DIRECTORY = Path(__file__).resolve().parent / "results"

# Physical model is deliberately fixed: W7 [[56,12,7]], R=7, idle=0,
# simultaneous same-shot X/Z, original boundaries, full DEM, Tesseract beam=2.
P_VALUES = [0.003, 0.0025]  # Descending: 0.30%, 0.25%.
SIMULATION_WORKERS = 380  # Automatically limited by available CPUs and RAM.
RELATIVE_1SIGMA = 0.20
MIN_SHOTS = 1000
MAX_SHOTS = 100_000_000
BATCH_SIZE = 64
SAMPLE_SEED = 20260923  # Independent of the historical selection campaign.

# Lower-bound screening first, then exact weight <=6 SAT decision tasks.
SCREENING_WORKERS = 32
SCREENING_PROOF_SECONDS = 300
HEURISTIC_SECONDS = 10  # Finds upper bounds only.
DISTANCE_WORKERS = 380
DISTANCE_SHARDS = 16  # 2 CSS sectors x 12 logical bits x 16 = 384 tasks/code.
SAT_SECONDS_PER_TASK = 300  # A timeout is UNKNOWN, never a proof of distance 7.
SAT_SOLVER = "minisat22"  # Also supported: "glucose3", "glucose4".
RETRY_UNKNOWN = False  # Raise time budget, or set True, to retry unresolved tasks.
