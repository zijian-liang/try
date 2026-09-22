"""Editable defaults for the joint X/Z, full circuit-level campaign."""

CODE_NAMES = ("w7", "bb72", "surface7")
P_VALUES = (0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.005, 0.006)
ROUNDS = None  # None: use the static code distance for each code.
ROUNDS_BY_CODE = {"w7": 7, "bb72": 6, "surface7": 7}
IDLE_SCALE = 0.0
WORKERS = 376
RELATIVE_ERROR = 0.10
MIN_SHOTS = 1000
MAX_SHOTS = 100_000_000
BATCH_SIZE = 64
RUN_SEED = 20260920
OUTPUT = "results"

# This is a fixed beam=2 experiment. Beam climbing is deliberately disabled.
# The source DEM contains both CSS sectors and is never graph-decomposed.
DECODER_SETTINGS = {
    "det_beam": 2,
    "beam_climbing": False,
    "no_revisit_dets": True,
    "merge_errors": True,
    "pqlimit": 200_000,
    "sparsify_errors": False,
    "num_det_orders": 20,
    "seed": 2384753,
}
