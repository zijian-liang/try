"""Tesseract adapter that checks every recovery and never discards a shot."""
from __future__ import annotations

import os

for _thread_variable in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import numpy as np
import stim
from tesseract_decoder import tesseract


def _mask(indices):
    answer = 0
    for i in indices:
        answer ^= 1 << int(i)
    return answer


def bool_mask(bits):
    return int.from_bytes(np.packbits(bits, bitorder="little").tobytes(), "little")


def empty_counts(k):
    return {
        "N": 0, "Fjoint": 0, "FX": 0, "FZ": 0, "Fboth": 0,
        "decodefail": 0, "lowconfidence": 0,
        "invalid_syndrome": 0, "decoder_exceptions": 0,
        "logical_X": [0] * k, "logical_Z": [0] * k,
    }


class CheckedTesseract:
    """Recover from both detector families, checking syndrome consistency.

    A low-confidence prediction is accepted if it is a valid recovery; its
    logical residual is counted normally. Invalid recoveries and exceptions
    are recorded separately as operational decoder failures. They contribute
    to Fjoint, but are never invented as X, Z, or both logical residuals.
    """

    def __init__(self, dem: stim.DetectorErrorModel, settings: dict):
        self.decoder = tesseract.TesseractConfig(dem=dem, **settings).compile_decoder()
        self.num_detectors = dem.num_detectors
        self.num_observables = dem.num_observables
        errors = self.decoder.errors
        self.det_masks = [_mask(e.symptom.detectors) for e in errors]
        self.obs_masks = [_mask(e.symptom.observables) for e in errors]

    def decode_checked(self, syndrome):
        target = bool_mask(syndrome)
        if target == 0:
            return 0, True, False, None
        try:
            self.decoder.decode_to_errors(syndrome)
            lowconfidence = bool(self.decoder.low_confidence_flag)
            predicted_detectors = 0
            predicted_observables = 0
            for error_id in self.decoder.predicted_errors_buffer:
                predicted_detectors ^= self.det_masks[error_id]
                predicted_observables ^= self.obs_masks[error_id]
            valid = predicted_detectors == target
            return predicted_observables, valid, lowconfidence, None if valid else "invalid_syndrome"
        except (RuntimeError, ValueError, IndexError) as exc:
            # No fallback decoder, no postselection, and no silently successful
            # zero prediction. Record the cause in the point's diagnostics.
            return 0, False, False, type(exc).__name__ + ": " + str(exc)[:400]

    def count_batch(self, detectors, observables, k):
        if observables.shape[1] != 2 * k:
            raise ValueError(f"Expected 2*k={2*k} same-shot observables, got {observables.shape}")
        counts = empty_counts(k)
        counts["N"] = len(detectors)
        exception_examples = []
        sector_mask = (1 << k) - 1
        for syndrome, truth in zip(detectors, observables):
            prediction, valid, lowconfidence, reason = self.decode_checked(syndrome)
            counts["lowconfidence"] += int(lowconfidence)
            if not valid:
                counts["decodefail"] += 1
                counts["Fjoint"] += 1
                if reason == "invalid_syndrome":
                    counts["invalid_syndrome"] += 1
                else:
                    counts["decoder_exceptions"] += 1
                    if reason not in exception_examples and len(exception_examples) < 3:
                        exception_examples.append(reason)
                continue
            residual = bool_mask(truth) ^ prediction
            residual_x = residual & sector_mask
            residual_z = (residual >> k) & sector_mask
            x, z = bool(residual_x), bool(residual_z)
            counts["FX"] += int(x)
            counts["FZ"] += int(z)
            counts["Fboth"] += int(x and z)
            counts["Fjoint"] += int(x or z)
            for j in range(k):
                counts["logical_X"][j] += (residual_x >> j) & 1
                counts["logical_Z"][j] += (residual_z >> j) & 1
        return counts, exception_examples
