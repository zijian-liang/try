"""Generate reproducible, distinct W7 direction-order candidates.

This is a finite search over translation-invariant, unflagged extraction
circuits.  The CNOT count stays 392 per cycle.  Eight/nine-layer parallel
orders dominate; a few fourteen-layer X-then-Z or Z-then-X orders provide
serial controls.  Historical search results select seed orders only and
are never treated as current circuit-distance certificates.

Each accepted schedule passes edge-coverage/collision/extraction-parity
checks, a noiseless Stim sample, and canonical full-DEM deduplication at
R=7, p=.003, idle_scale=0.  Deduplication removes identical labelled DEMs;
it does not solve graph isomorphism or prove global schedule optimality.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random

from circuits import build_circuit, get_code, schedule_layers


GENERATOR_VERSION = "w7_direction_orders_v1_r7_dem"
DEM_ROUNDS = 7
DEM_PROBABILITY = 0.003


def schedule_hash(schedule):
    return hashlib.sha256(json.dumps(schedule,sort_keys=True,separators=(",",":")).encode()).hexdigest()


def canonical_dem_hash(circuit):
    """Hash the complete undecomposed, labelled detector-error distribution.

    Equal response columns combine using the parity probability formula.
    Rounding at 13 significant digits suppresses floating-point association
    noise; it is not a coarse structural/syndrome-degree approximation.
    Detector coordinates and emission order do not affect this fingerprint.
    """
    dem = circuit.detector_error_model(decompose_errors=False).flattened()
    parity_factors = {}
    for op in dem:
        if op.type != "error":
            continue
        detectors = logicals = 0
        for target in op.targets_copy():
            if target.is_relative_detector_id():
                detectors ^= 1 << target.val
            elif target.is_logical_observable_id():
                logicals ^= 1 << target.val
            elif target.is_separator():
                raise ValueError("Unexpected decomposition separator in full joint DEM")
        key = (detectors,logicals)
        parity_factors[key] = parity_factors.get(key,1.0)*(1-2*op.args_copy()[0])
    rows = [(hex(d),hex(l),format((1-factor)/2,".13g"))
            for (d,l),factor in sorted(parity_factors.items())]
    payload = [dem.num_detectors,dem.num_observables,rows]
    return hashlib.sha256(json.dumps(payload,separators=(",",":")).encode()).hexdigest()


def _fast_legal(schedule):
    """Algebraic W7 check for a translation-invariant direction schedule."""
    xs,zs = schedule["x"],schedule["z"]
    if len(xs) != len(zs) or any(x is None and z is None for x,z in zip(xs,zs)):
        return False
    if sorted(x for x in xs if x is not None) != list(range(7)):
        return False
    if sorted(z for z in zs if z is not None) != list(range(7)):
        return False
    for x,z in zip(xs,zs):
        if x is not None and z is not None:
            # XA and ZB use A data; XB and ZA use B data.
            if (x<3) == (z<4):
                return False
    tx = [xs.index(i) for i in range(7)]
    tz = [zs.index(i) for i in range(7)]
    return all((tx[a]<tz[b]) == (tx[b+3]<tz[a+4])
               for a in range(3) for b in range(4))


def _history_by_depth(seed):
    manifest = json.loads(Path(__file__).with_name("seed_schedules.json").read_text())
    grouped = {8:[],9:[]}
    for item in manifest["seeds"]:
        schedule = item["schedule"]
        depth = len(schedule["x"])
        if depth not in grouped or not _fast_legal(schedule):
            continue
        # Stronger historical rescreening first, then its still-unproven
        # heuristic upper bound.  A seeded hash gives a stable tie ordering.
        random_key = hashlib.sha256(f"{seed}:{schedule_hash(schedule)}".encode()).hexdigest()
        priority = (-int(item.get("historical_stronger_rescreen",False)),
                    -int(item.get("historical_heuristic_upper_bound",0)),random_key)
        grouped[depth].append((priority,schedule))
    return {depth:[schedule for _,schedule in sorted(items,key=lambda item:item[0])]
            for depth,items in grouped.items()}


def _parallel_stream(depth, seed, histories):
    yield from histories
    rng = random.Random(f"{GENERATOR_VERSION}:{seed}:parallel:{depth}")
    unsuccessful = 0
    while True:
        xs = list(range(7))+[None]*(depth-7)
        rng.shuffle(xs)
        a_slots = [t for t,x in enumerate(xs) if x is None or x<3]
        z_a_slots = rng.sample(a_slots,3)
        b_slots = [t for t,x in enumerate(xs) if (x is None or x>=3) and t not in z_a_slots]
        if len(b_slots)<4:
            continue
        z_b_slots = rng.sample(b_slots,4)
        zs = [None]*depth
        for t,direction in zip(z_a_slots,rng.sample(range(4,7),3)):
            zs[t] = direction
        for t,direction in zip(z_b_slots,rng.sample(range(4),4)):
            zs[t] = direction
        schedule = {"x":xs,"z":zs}
        if _fast_legal(schedule):
            unsuccessful = 0
            yield schedule
        else:
            unsuccessful += 1
            if unsuccessful > 2_000_000:
                raise RuntimeError(f"Unable to find another legal {depth}-layer schedule")


def _serial_stream(seed):
    rng = random.Random(f"{GENERATOR_VERSION}:{seed}:serial")
    number = 0
    while True:
        xs = rng.sample(range(7),7)
        zs = rng.sample(range(7),7)
        if number % 2 == 0:
            yield {"x":xs+[None]*7,"z":[None]*7+zs}
        else:
            yield {"x":[None]*7+xs,"z":zs+[None]*7}
        number += 1


def _family(schedule):
    depth = len(schedule["x"])
    if depth != 14:
        return f"parallel{depth}"
    return "serial14_x_then_z" if schedule["x"][0] is not None else "serial14_z_then_x"


def generate_candidates(count=64, seed=20260921):
    """Return a manifest; candidate IDs depend on family+schedule, not count.

    Results for a smaller count are a prefix of a larger request at the same
    seed.  At count=64 there are 30 eight-layer schedules (baseline included),
    30 nine-layer schedules and 4 fourteen-layer serial controls.  Every
    returned candidate has a distinct canonical full DEM at the fixed R7
    reference model.  This fingerprint is not a distance result.
    """
    if not isinstance(count,int) or isinstance(count,bool) or count<1:
        raise ValueError("count must be a positive integer, including baseline")
    if not isinstance(seed,int) or isinstance(seed,bool):
        raise ValueError("seed must be an integer")
    code = get_code("w7")
    histories = _history_by_depth(seed)
    streams = {depth:_parallel_stream(depth,seed,histories[depth]) for depth in (8,9)}
    streams[14] = _serial_stream(seed)
    candidates = []
    seen_schedules = set()
    seen_dems = set()
    rejected = Counter()

    def accept(schedule, baseline=False):
        fingerprint = schedule_hash(schedule)
        if fingerprint in seen_schedules:
            rejected["identical_schedule"] += 1
            return False
        seen_schedules.add(fingerprint)
        if not _fast_legal(schedule):
            raise ValueError("Candidate failed its W7 algebraic legality check")
        schedule_layers(code,schedule)  # independent full-matrix validation
        noisy,metadata = build_circuit("w7",DEM_PROBABILITY,rounds=DEM_ROUNDS,
                                       idle_scale=0.0,schedule=schedule)
        dem_fingerprint = canonical_dem_hash(noisy)
        if dem_fingerprint in seen_dems:
            rejected["identical_labelled_dem"] += 1
            return False
        ideal,_ = build_circuit("w7",0,rounds=DEM_ROUNDS,idle_scale=0.0,schedule=schedule)
        detectors,observables = ideal.compile_detector_sampler(seed=int(fingerprint[:8],16)).sample(
            1,separate_observables=True)
        if detectors.any() or observables.any():
            raise ValueError("Candidate has non-deterministic/noiseless syndrome or logical labels")
        seen_dems.add(dem_fingerprint)
        family = "baseline8" if baseline else _family(schedule)
        candidates.append({
            "id":"w7_baseline" if baseline else f"{family}_{fingerprint[:12]}",
            "code":"w7","family":family,"schedule":schedule,
            "schedule_sha256":fingerprint,"is_baseline":baseline,
            "dem_sha256":dem_fingerprint,"cnot_layers":metadata["cnot_layers_per_cycle"],
            "cnots_per_cycle":metadata["cnots_per_cycle"],
            "current_distance_status":"requires this campaign's distance screening",
        })
        return True

    if not accept(code["schedule"],baseline=True):
        raise AssertionError("Baseline was not accepted")
    parallel_position = 0
    while len(candidates)<count:
        position = len(candidates)
        if position % 15 == 0:
            depth = 14
        else:
            depth = 9 if parallel_position % 2 == 0 else 8
            parallel_position += 1
        attempts = 0
        while not accept(next(streams[depth])):
            attempts += 1
            if attempts > 10000:
                raise RuntimeError(f"Exhausted distinct {depth}-layer candidates; request fewer candidates")
    return {
        "schema_version":1,"generator_version":GENERATOR_VERSION,"seed":seed,
        "candidates":candidates,
        "family_counts":dict(Counter(item["family"] for item in candidates)),
        "semantic_deduplication":{
            "method":"canonical full undecomposed labelled DEM hash; parity-merge duplicate response columns",
            "rounds":DEM_ROUNDS,"p":DEM_PROBABILITY,"idle_scale":0.0,
            "rejected":dict(rejected),"does_not_identify_all_graph_isomorphisms":True,
        },
        "scope":"Finite seeded search of legal translation-invariant 8/9-layer no-flag orders, plus 14-layer serial controls; 392 CNOTs per cycle. Not an exhaustive search or a proof of schedule optimality.",
        "distance_warning":"Historical upper bounds select seed orders only. Every candidate requires current R7 distance screening; no historical upper bound is a lower bound or a current certificate.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count",type=int,default=64)
    parser.add_argument("--seed",type=int,default=20260921)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    payload = json.dumps(generate_candidates(args.count,args.seed),indent=2)+"\n"
    if args.output:
        args.output.write_text(payload)
    else:
        print(payload,end="")


if __name__ == "__main__":
    main()

