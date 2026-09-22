"""Independent Pauli-frame checks for the joint-memory Stim circuits.

Run ``python validate_circuits.py``.  An optional ``--ibm-source`` directory
containing the author's decoder_setup.py also checks the BB72 operation cycle
against the original source.  Normal validation has no source-ZIP dependency.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path
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


def validate_one(name, random_histories, seed):
    started = time.time()
    code = get_code(name)
    ideal,meta = build_circuit(name,0,rounds=code["static_distance"])
    ds,obs = ideal.compile_detector_sampler(seed=seed).sample(64,separate_observables=True)
    if ds.any() or obs.any():
        raise AssertionError(f"{name}: nonzero no-fault detector or observable")
    # Include data idles in this validation even though campaign default is 0.
    circuit,meta = build_circuit(name,0.003,rounds=code["static_distance"],idle_scale=1)
    circuit.detector_error_model(decompose_errors=False)
    locations = fault_locations(circuit)
    rng = np.random.default_rng(seed)
    histories = []
    for _ in range(random_histories):
        count = int(rng.integers(1,13))
        chosen = rng.choice(len(locations),size=count,replace=False)
        histories.append({int(i):str(rng.choice(locations[int(i)][2])) for i in chosen})
    # Exercise all 15 exclusive two-qubit outcomes at one CNOT location.
    cnot_location = next(i for i,item in enumerate(locations) if item[4] == "CNOT")
    histories.extend({cnot_location:pauli} for pauli in locations[cnot_location][2])
    # In particular, check every last-noisy-cycle Z-preparation fault: it is
    # measured in the first clean cycle and must not disappear at the boundary.
    z_ancillas = set(meta["physical_qubits"]["Z_ancillas"])
    z_preps = [i for i,item in enumerate(locations) if item[4] == "preparation" and item[1][0] in z_ancillas]
    final_prep_instruction = max(locations[i][0] for i in z_preps)
    final_z_preps = [i for i in z_preps if locations[i][0] == final_prep_instruction]
    histories.extend({i:"X"} for i in final_z_preps)
    both_sector_histories = 0
    checked_faults = 0
    for number,history in enumerate(histories):
        injected = inject_history(circuit,locations,history)
        expected_d,expected_o = classical_frame(injected,code["n"],code["lx"],code["lz"])
        actual_d,actual_o = injected.compile_detector_sampler(seed=seed+number+1).sample(1,separate_observables=True)
        if not np.array_equal(expected_d,actual_d[0]) or not np.array_equal(expected_o,actual_o[0]):
            raise AssertionError(f"{name}: Stim/frame mismatch for history {number}: {history}")
        if expected_o[:code["k"]].any() and expected_o[code["k"]:].any():
            both_sector_histories += 1
        checked_faults += len(history)
    if not both_sector_histories:
        raise AssertionError(f"{name}: test failed to exercise simultaneous logical X/Z labels")
    return {
        "passed":True,"random_histories":random_histories,"total_histories":len(histories),
        "physical_faults_in_histories":checked_faults,"all_cnot_pauli_outcomes":15,
        "last_noisy_z_preparation_histories":len(final_z_preps),
        "histories_with_both_logical_components":both_sector_histories,
        "no_fault_samples":64,"rounds":code["static_distance"],"idle_scale_tested":1,
        "schedule":meta["schedule"],"cnot_layers_per_cycle":meta["cnot_layers_per_cycle"],
        "validation_circuit_sha256":hashlib.sha256(str(circuit).encode()).hexdigest(),
        "seconds":round(time.time()-started,3),
        "method":"independent classical x/z Pauli propagation vs deterministic Stim fault injection; direct logical matrix cross-check",
    }


def verify_original_ibm(source):
    """Reconstruct the original IBM cycle, without importing its decoder stack."""
    path = Path(source)
    if path.is_dir():
        path = path/"decoder_setup.py"
    source_text = path.read_text()
    namespace = {"np":np,"ell":6,"m":6,"n":72,"n2":36,
                 "a1":3,"a2":1,"a3":2,"b1":3,"b2":1,"b3":2}
    for node in ast.parse(source_text).body:
        if isinstance(node,ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0],ast.Name):
            if node.targets[0].id in ("sX","sZ"):
                namespace[node.targets[0].id] = ast.literal_eval(node.value)
    first = source_text.index("# cyclic shift matrices")
    end = source_text.index("# number of logical qubits",first)
    exec(compile(source_text[first:end],str(path),"exec"),namespace)
    first = source_text.index("# Give a name to each qubit")
    end = source_text.index("# full syndrome measurement circuit",first)
    exec(compile(source_text[first:end],str(path),"exec"),namespace)
    code = get_code("bb72")
    assert np.array_equal(code["hx"],namespace["hx"])
    assert np.array_equal(code["hz"],namespace["hz"])

    def physical(node):
        kind,i = node
        return {"data_left":i,"data_right":36+i,"Xcheck":72+i,"Zcheck":108+i}[kind]

    original = [(operation[0],)+tuple(physical(q) for q in operation[1:]) for operation in namespace["cycle"]]
    # The original emits fixed-size batches for each of its eight time slots.
    lengths = [36+36+36]+[72]*5+[36+36+36]+[72+36+36]
    expected_slots = []
    offset = 0
    for length in lengths:
        expected_slots.append(Counter(original[offset:offset+length]))
        offset += length
    assert offset == len(original)
    circuit,meta = build_circuit("bb72",0.003,rounds=1,idle_scale=1)
    actual_slots = []
    current = []
    started = False
    for op in circuit:
        if not started:
            if op.name != "RX":
                continue
            started = True
        qs = [t.value for t in op.targets_copy()]
        if op.name == "TICK":
            actual_slots.append(Counter(current))
            current = []
            if len(actual_slots) == 8:
                break
        elif op.name == "CX":
            current += [("CNOT",a,b) for a,b in zip(qs[::2],qs[1::2])]
        elif op.name in ("R","RX","M","MX"):
            label = {"R":"PrepZ","RX":"PrepX","M":"MeasZ","MX":"MeasX"}[op.name]
            current += [(label,q) for q in qs]
        elif op.name == "DEPOLARIZE1":
            current += [("IDLE",q) for q in qs]
    if actual_slots != expected_slots:
        raise AssertionError("BB72 operation/time-slot structure differs from original IBM source")
    return {"passed":True,"code":"bb72","source_file":"decoder_setup.py",
            "matched_time_slots":8,"matched_operations":len(original),
            "cnots":sum(op[0]=="CNOT" for op in original),
            "data_idle_locations":sum(op[0]=="IDLE" for op in original),
            "method":"execute original shift-matrix and cycle-construction sections at ell=m=6; compare every operation in each time slot, including data idles"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-histories",type=int,default=128)
    parser.add_argument("--seed",type=int,default=730219)
    parser.add_argument("--ibm-source",type=Path)
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name("validation_circuits.json"))
    args = parser.parse_args()
    result = {"schema_version":1,"circuits":{}}
    for i,name in enumerate(("w7","bb72","surface7")):
        result["circuits"][name] = validate_one(name,args.random_histories,args.seed+i*10000)
        print(name,json.dumps(result["circuits"][name]),flush=True)
    if args.ibm_source:
        result["original_ibm_cycle"] = verify_original_ibm(args.ibm_source)
        print("original_ibm_cycle",json.dumps(result["original_ibm_cycle"]),flush=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(f"Saved {args.output}",flush=True)


if __name__ == "__main__":
    main()
