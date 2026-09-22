"""Independent validation of generated weight-7 schedule-scan circuits.

Run ``python validate_scan.py`` for the default representative schedules.
Use ``--manifest candidates.json`` to validate a saved scan manifest instead.
Every check uses R=7, idle noise zero by default.  Fault-history validation
compares Stim with a separate classical x/z Pauli-frame implementation and
checks logical labels directly from the binary logical-operator matrices.
No Monte Carlo logical-error-rate claim is made by this validation.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import tempfile
import time
import numpy as np
import stim
from circuits import build_circuit, get_code

def fault_locations(circuit):
    """Enumerate physical channel locations, with mutually exclusive Paulis."""
    result = []
    for index,op in enumerate(circuit):
        args = op.gate_args_copy()
        if not args or not args[0]:
            continue
        targets = [t.value for t in op.targets_copy()]
        if op.name == "DEPOLARIZE2":
            for a,b in zip(targets[::2],targets[1::2]):
                result.append((index, (a,b), tuple(x+y for x in "IXYZ" for y in "IXYZ" if x+y != "II"), False, "CNOT"))
        elif op.name == "DEPOLARIZE1":
            result += [(index,(q,),tuple("XYZ"),False,"idle") for q in targets]
        elif op.name in ("X_ERROR","Y_ERROR","Z_ERROR"):
            result += [(index,(q,),(op.name[0],),False,"preparation") for q in targets]
        elif op.name in ("M","MX"):
            result += [(index,(q,),(("X" if op.name == "M" else "Z"),),True,"readout") for q in targets]
    return result

def inject_history(circuit, locations, history):
    """Replace all stochastic channels by one specified physical fault history.

    X_ERROR(1), rather than an ideal X gate, keeps the intended no-fault
    reference when Stim compiles its detector sampler.
    """
    by_instruction = {}
    for location,pauli in history.items():
        index,qubits,choices,before,kind = locations[location]
        if pauli not in choices:
            raise ValueError("Fault is not allowed at this location")
        by_instruction.setdefault(index,[]).extend(zip(qubits,pauli))
    out = stim.Circuit()
    noise_names = {"DEPOLARIZE1","DEPOLARIZE2","X_ERROR","Y_ERROR","Z_ERROR"}
    for index,op in enumerate(circuit):
        faults = by_instruction.get(index,())
        if op.name in noise_names or op.name in ("M","MX"):
            for q,pauli in faults:
                if pauli != "I":
                    out.append(pauli+"_ERROR",[q],1)
            if op.name in ("M","MX"):
                out.append(op.name,op.targets_copy(),0)
        else:
            out.append(op)
    return out

def classical_frame(circuit, n_data, lx, lz):
    """Propagate x/z bit vectors without any Stim simulation machinery."""
    xs = np.zeros(circuit.num_qubits,dtype=np.uint8)
    zs = np.zeros_like(xs)
    records = []
    detectors = []
    observables = np.zeros(circuit.num_observables,dtype=np.uint8)
    analytic_observables = None
    for op in circuit:
        name = op.name
        ts = op.targets_copy()
        qs = [t.value for t in ts]
        if name in ("H", "S", "S_DAG"):
            if name == "H":
                old = xs[qs].copy()
                xs[qs] = zs[qs]
                zs[qs] = old
            else:
                zs[qs] ^= xs[qs]
        elif name == "CX":
            for a,b in zip(qs[::2],qs[1::2]):
                xs[b] ^= xs[a]
                zs[a] ^= zs[b]
        elif name in ("R","RX"):
            xs[qs] = 0
            zs[qs] = 0
        elif name in ("X_ERROR","Y_ERROR","Z_ERROR"):
            assert op.gate_args_copy() == [1.0]
            if name[0] in ("X","Y"):
                xs[qs] ^= 1
            if name[0] in ("Z","Y"):
                zs[qs] ^= 1
        elif name in ("M","MX"):
            records.extend(map(int,xs[qs] if name == "M" else zs[qs]))
        elif name == "MPP":
            if analytic_observables is None:
                analytic_observables = np.concatenate(((lz@xs[:n_data])&1,(lx@zs[:n_data])&1))
            # Every generated MPP instruction contains one logical/ref product.
            parity = 0
            for target in ts:
                if target.is_combiner:
                    continue
                q = target.value
                if target.is_x_target or target.is_y_target:
                    parity ^= int(zs[q])
                if target.is_z_target or target.is_y_target:
                    parity ^= int(xs[q])
            records.append(parity)
        elif name in ("DETECTOR","OBSERVABLE_INCLUDE"):
            parity = 0
            for target in ts:
                assert target.is_measurement_record_target
                parity ^= records[target.value]
            if name == "DETECTOR":
                detectors.append(parity)
            else:
                observables[int(op.gate_args_copy()[0])] ^= parity
        elif name in ("TICK","QUBIT_COORDS","SHIFT_COORDS"):
            pass
        else:
            raise ValueError(f"Unsupported operation in independent frame check: {name}")
    if not np.array_equal(observables,analytic_observables):
        raise AssertionError("Logical/ref measurements disagree with direct residual Pauli logical labels")
    return np.array(detectors,dtype=np.uint8),observables


def check_first_cycle_edges(circuit, code):
    """Reconstruct check supports from actual first-cycle CNOT instructions."""
    n, nx = code["n"], len(code["hx"])
    xanc = set(range(n, n+nx))
    zanc = set(range(n+nx, n+nx+len(code["hz"])))
    started = False
    edges = Counter()
    cnot_batches = 0
    first_readout_instruction = None
    for index, op in enumerate(circuit):
        targets = [t.value for t in op.targets_copy()]
        if not started:
            if op.name == "RX" and set(targets) == xanc:
                started = True
            continue
        if op.name in ("M", "MX"):
            first_readout_instruction = index
            break
        if op.name != "CX":
            continue
        cnot_batches += 1
        if len(targets) != len(set(targets)):
            raise AssertionError("One CNOT layer reuses a physical qubit")
        for control, target in zip(targets[::2], targets[1::2]):
            if control in xanc and 0 <= target < n:
                edges[("X", control-n, target)] += 1
            elif target in zanc and 0 <= control < n:
                edges[("Z", target-n-nx, control)] += 1
            else:
                raise AssertionError("Invalid ancilla/data CNOT orientation")
    expected = Counter()
    for sector, matrix in (("X", code["hx"]), ("Z", code["hz"])):
        for check, datum in zip(*np.nonzero(matrix)):
            expected[(sector, int(check), int(datum))] += 1
    if edges != expected:
        raise AssertionError("Actual CNOT edges disagree with the check matrices")
    if sum(edges.values()) != 392 or first_readout_instruction is None:
        raise AssertionError("Unexpected number of weight-7 check edges")
    return {"edges_per_cycle": sum(edges.values()),
            "nonempty_cnot_layers": cnot_batches,
            "first_readout_instruction": first_readout_instruction}


def validate_schedule_case(label, schedule, *, random_histories=50, seed=730219,
                           rounds=7):
    started = time.monotonic()
    code = get_code("w7")
    ideal, ideal_meta = build_circuit("w7", 0, rounds=rounds,
                                      idle_scale=0, schedule=schedule)
    ds, obs = ideal.compile_detector_sampler(seed=seed).sample(
        64, separate_observables=True)
    if ds.any() or obs.any():
        raise AssertionError(f"{label}: nonzero no-fault detector or observable")
    edge_check = check_first_cycle_edges(ideal, code)
    circuit, meta = build_circuit("w7", 0.003, rounds=rounds,
                                  idle_scale=0, schedule=schedule)
    circuit.detector_error_model(decompose_errors=False)
    noisy_edge_check = check_first_cycle_edges(circuit, code)
    locations = fault_locations(circuit)
    counts = Counter(item[4] for item in locations)
    expected_counts = {"CNOT": rounds*392, "preparation": rounds*56,
                       "readout": rounds*56}
    if dict(counts) != expected_counts:
        raise AssertionError(f"{label}: noise-location mismatch {dict(counts)}")
    if circuit.num_detectors != 56*(rounds+2) or circuit.num_observables != 24:
        raise AssertionError("Incorrect detector or simultaneous logical-label count")
    if meta["noiseless_tail_rounds"] != 2 or meta["idle_scale"] != 0:
        raise AssertionError("Incorrect campaign boundary or idle model")

    rng = np.random.default_rng(seed)
    histories = []
    for _ in range(random_histories):
        count = int(rng.integers(1, 13))
        chosen = rng.choice(len(locations), size=count, replace=False)
        histories.append({int(i): str(rng.choice(locations[int(i)][2]))
                          for i in chosen})

    # Each direction/layer and both CNOT orientations are exercised.  At a
    # selected physical location the fifteen outcomes are mutually exclusive.
    representatives = {}
    for i, (instruction, qubits, choices, before, kind) in enumerate(locations):
        if kind != "CNOT" or instruction >= noisy_edge_check["first_readout_instruction"]:
            continue
        sector = "X" if qubits[0] >= code["n"] else "Z"
        representatives.setdefault((instruction, sector), i)
    if len(representatives) != 14:
        raise AssertionError("Expected seven CNOT directions in each CSS sector")
    for location in representatives.values():
        if len(locations[location][2]) != 15:
            raise AssertionError("A CNOT channel lost a Pauli outcome")
        histories.extend({location: pauli} for pauli in locations[location][2])

    z_ancillas = set(meta["physical_qubits"]["Z_ancillas"])
    z_preps = [i for i, item in enumerate(locations)
               if item[4] == "preparation" and item[1][0] in z_ancillas]
    last_prep = max(locations[i][0] for i in z_preps)
    last_z_preps = [i for i in z_preps if locations[i][0] == last_prep]
    histories.extend({i: "X"} for i in last_z_preps)
    checked_faults = both_components = 0
    for index, history in enumerate(histories):
        injected = inject_history(circuit, locations, history)
        expected_d, expected_o = classical_frame(
            injected, code["n"], code["lx"], code["lz"])
        actual_d, actual_o = injected.compile_detector_sampler(
            seed=seed+index+1).sample(1, separate_observables=True)
        if not np.array_equal(expected_d, actual_d[0]):
            raise AssertionError(f"{label}: detector mismatch, history {index}: {history}")
        if not np.array_equal(expected_o, actual_o[0]):
            raise AssertionError(f"{label}: logical-label mismatch, history {index}: {history}")
        if expected_o[:code["k"]].any() and expected_o[code["k"]:].any():
            both_components += 1
        checked_faults += len(history)
    if not both_components:
        raise AssertionError(f"{label}: no test exercised both logical X/Z components")
    return {
        "passed": True, "label": label, "schedule": meta["schedule"],
        "rounds": rounds, "idle_scale": 0, "noiseless_tail_rounds": 2,
        "no_fault_samples": 64, "random_histories": random_histories,
        "cnot_locations_with_all_15_outcomes": len(representatives),
        "single_cnot_fault_histories": 15*len(representatives),
        "last_noisy_z_preparation_histories": len(last_z_preps),
        "total_histories": len(histories),
        "physical_faults_in_histories": checked_faults,
        "histories_with_both_logical_components": both_components,
        "physical_noise_location_counts": dict(counts),
        "edges_per_cycle": edge_check["edges_per_cycle"],
        "cnot_layers_per_cycle": meta["cnot_layers_per_cycle"],
        "num_detectors": circuit.num_detectors, "num_observables": circuit.num_observables,
        "ideal_circuit_sha256": hashlib.sha256(str(ideal).encode()).hexdigest(),
        "noisy_circuit_sha256": hashlib.sha256(str(circuit).encode()).hexdigest(),
        "seconds": round(time.monotonic()-started, 3),
        "method": "independent classical x/z Pauli propagation versus deterministic Stim fault injection; direct logical-matrix parity cross-check",
    }


def representative_cases():
    """Select real generated candidates, including both serial orientations."""
    from candidates import generate_candidates
    manifest = generate_candidates(count=31, seed=20260921)
    selected = {}
    for item in manifest["candidates"]:
        selected.setdefault(item["family"], item)
    expected = {"baseline8", "parallel8", "parallel9",
                "serial14_x_then_z", "serial14_z_then_x"}
    if set(selected) != expected:
        raise AssertionError("Generated validation set missed a schedule family")
    return [(item["id"], item["schedule"]) for item in selected.values()]


def validate_simulation_resume():
    """Real two-worker Tesseract smoke test, including tighter resume limits."""
    from simulation import simulate_point
    candidate = {"id": "validation_baseline", "family": "baseline8",
                 "is_baseline": True, "code": "w7",
                 "schedule": get_code("w7")["schedule"]}
    audit = json.loads(Path(__file__).with_name("baseline_distance.json").read_text())["codes"]["w7"]
    with tempfile.TemporaryDirectory(prefix="w7_scan_validation_") as directory:
        first = simulate_point(candidate, audit, p=0.003, output=Path(directory),
                               workers=2, min_shots=8, max_shots=8,
                               relative_error=0.1, batch_size=4, run_seed=912860)
        # JSON round trip ensures the following comparisons cannot accidentally
        # observe an in-place mutation of the first run's result dictionary.
        first = json.loads(json.dumps(first))
        second = simulate_point(candidate, audit, p=0.003, output=Path(directory),
                                workers=2, min_shots=16, max_shots=16,
                                relative_error=0.05, batch_size=8, run_seed=912860)
        second = json.loads(json.dumps(second))
        third = simulate_point(candidate, audit, p=0.003, output=Path(directory),
                               workers=2, min_shots=16, max_shots=16,
                               relative_error=0.05, batch_size=8, run_seed=912860)
        if first["counts"]["N"] != 8 or second["counts"]["N"] != 16:
            raise AssertionError("Resume did not increase N from exactly 8 to exactly 16")
        if first["identity"] != second["identity"] or first["experiment_id"] != second["experiment_id"]:
            raise AssertionError("Changing runtime stopping/batch settings changed experiment identity")
        if second["counts"] != third["counts"] or second["ledger"] != third["ledger"]:
            raise AssertionError("A completed point was sampled or counted again on resume")
        for result in (first, second, third):
            counts = result["counts"]
            if counts["Fjoint"] != counts["FX"]+counts["FZ"]-counts["Fboth"]+counts["decodefail"]:
                raise AssertionError("Same-shot X/Z union count mismatch")
            if result["ledger"]["pending"] or result["ledger"]["active_wave"] is not None:
                raise AssertionError("A completed run retained a pending partial wave")
            if not result["point_done"]:
                raise AssertionError("Exact shot-cap completion was not recorded")
            if result["runtime"]["workers"] != 2:
                raise AssertionError("This smoke test requires two actual workers")
            if sum(wave["planned_shots"] for wave in result["wave_history"]) != counts["N"]:
                raise AssertionError("Committed-wave shot totals disagree with accumulated shots")
        for key, value in first["counts"].items():
            newer = second["counts"][key]
            if isinstance(value, list):
                if any(b < a for a, b in zip(value, newer)):
                    raise AssertionError("Resume lost per-observable counts")
            elif newer < value:
                raise AssertionError("Resume lost accumulated counts")
        return {
            "passed": True, "actual_workers": 2, "p": 0.003,
            "rounds": 7, "idle_scale": 0, "decoder": "tesseract",
            "decoder_settings": second["identity"]["decoder_settings"],
            "experiment_id": second["experiment_id"],
            "counts_at_8": first["counts"], "counts_at_16": second["counts"],
            "batch_size_change": [4, 8], "relative_error_change": [0.1, 0.05],
            "identity_unchanged": True, "completed_resume_idempotent": True,
            "complete_wave_accounting": True, "same_shot_union_count_consistent": True,
            "wave_history_at_16": second["wave_history"],
            "versions": second["identity"]["versions"],
            "note": "Functional smoke test only: sixteen shots do not estimate a useful logical error rate.",
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-histories", type=int, default=50)
    parser.add_argument("--seed", type=int, default=730219)
    parser.add_argument("--rounds", type=int, default=7)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--simulation-smoke", action="store_true",
                        help="Also run real two-worker Tesseract 8-to-16-shot resume checks")
    parser.add_argument("--simulation-only", action="store_true",
                        help="Run only the real simulation smoke test; choose a separate --output")
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("validation_scan.json"))
    args = parser.parse_args()
    if args.random_histories < 1:
        parser.error("--random-histories must be positive")
    cases = [] if args.simulation_only or args.manifest else representative_cases()
    if args.manifest:
        content = json.loads(args.manifest.read_text())
        entries = content.get("candidates", content) if isinstance(content, dict) else content
        if isinstance(entries, dict):
            entries = list(entries.values())
        cases = [(str(item.get("candidate_id", item.get("id", i))), item["schedule"])
                 for i, item in enumerate(entries)]
    report = {"schema_version": 1, "passed": False,
              "scope": "circuit correctness only; no circuit-distance or Monte Carlo precision claim",
              "stim_version": stim.__version__, "numpy_version": np.__version__,
              "circuits": {}}
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    for i, (label, schedule) in enumerate(cases):
        result = validate_schedule_case(label, schedule,
                                        random_histories=args.random_histories,
                                        seed=args.seed+10000*i, rounds=args.rounds)
        report["circuits"][label] = result
        args.output.write_text(json.dumps(report, indent=2)+"\n")
        print(f"{label}: PASS; {result['total_histories']} histories, "
              f"{result['cnot_layers_per_cycle']} CNOT layers", flush=True)
    if args.simulation_smoke or args.simulation_only:
        report["simulation_resume_smoke"] = validate_simulation_resume()
    report["passed"] = True
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(f"Saved {args.output.resolve()}", flush=True)


if __name__ == "__main__":
    main()
