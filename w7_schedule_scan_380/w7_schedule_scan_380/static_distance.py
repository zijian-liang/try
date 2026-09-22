"""Exact small static CSS distance verification by a meet-in-the-middle search.

This is a code-distance check, NOT a circuit-distance certificate.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np


def _columns(h, logical):
    h = np.asarray(h, dtype=np.uint8)
    logical = np.asarray(logical, dtype=np.uint8)
    syndromes = [sum(int(h[r, q]) << r for r in range(len(h))) for q in range(h.shape[1])]
    labels = [sum(int(logical[r, q]) << r for r in range(len(logical))) for q in range(h.shape[1])]
    return syndromes, labels


def verify_sector(h, logical):
    """Prove no logical of weight <=6, or return a counterexample; find <=7.

    All subsets of size <=3 are paired at equal syndrome. Any support of
    weight <=6 splits this way. Distinct logical labels distinguish a logical
    operator from a stabilizer. For a weight-7 witness, pair size4 with <=3.
    """
    s_cols, l_cols = _columns(h, logical)
    n = len(s_cols)
    entries = {}
    pairs = {}
    below_six = False
    best = None
    for weight in range(4):
        for subset in itertools.combinations(range(n), weight):
            s = l = mask = 0
            for q in subset:
                s ^= s_cols[q]
                l ^= l_cols[q]
                mask ^= 1 << q
            if s in pairs and pairs[s][0] != l:
                below_six = True
            if weight <= 2 and s not in pairs:
                pairs[s] = (l, mask)
            if s in entries:
                old_l, old_mask = entries[s]
                if old_l != l:
                    candidate = old_mask ^ mask
                    if best is None or candidate.bit_count() < best.bit_count():
                        best = candidate
            else:
                entries[s] = (l, mask)
    if best is not None:
        # The minimum encountered relative to one representative need not be
        # globally minimal. Claim only the counterexample upper bound.
        exact_six = best.bit_count() == 6 and not below_six
        return {"lower_bound": 6 if exact_six else 1, "upper_bound": best.bit_count(), "exact": exact_six,
                "witness": [q for q in range(n) if best >> q & 1],
                "method": "exhaustive exclusion <=5; explicit weight6 logical" if exact_six else "exhaustive <=3 + <=3; counterexample found"}
    for subset in itertools.combinations(range(n), 4):
        s = l = mask = 0
        for q in subset:
            s ^= s_cols[q]
            l ^= l_cols[q]
            mask ^= 1 << q
        if s in entries and entries[s][0] != l:
            best = entries[s][1] ^ mask
            assert best.bit_count() == 7
            return {"lower_bound": 7, "upper_bound": 7, "exact": True,
                    "witness": [q for q in range(n) if best >> q & 1],
                    "method": "exhaustive exclusion <=6; explicit weight7 logical"}
    return {"lower_bound": 8, "upper_bound": None, "exact": False,
            "method": "exhaustive exclusion <=7"}


def main():
    from circuits import get_code
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--codes", nargs="+", default=["w7", "bb72", "surface7"])
    parser.add_argument("--output", default="static_distance_results.json")
    args = parser.parse_args()
    result = {}
    for name in args.codes:
        c = get_code(name)
        def field(key):
            return c[key] if isinstance(c, dict) else getattr(c, key)
        result[name] = {"X": verify_sector(field("hz"), field("lz")),
                        "Z": verify_sector(field("hx"), field("lx"))}
        print(name, json.dumps(result[name]), flush=True)
    Path(args.output).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
