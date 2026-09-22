"""CSS memory circuits with simultaneous logical X/Z Pauli-frame labels.

The k noiseless reference qubits are simulation bookkeeping, not hardware.
Encoded Bell pairs make both logical/ref XX and ZZ observables deterministic
and commuting.  Sampling the circuit (not its DEM) preserves the exclusive
15-outcome depolarizing channel at each two-qubit gate.
"""
from __future__ import annotations

import copy
from functools import lru_cache
from typing import Any

import numpy as np
import stim


def gf2_rref(a):
    a = np.array(a, dtype=np.uint8, copy=True) & 1
    pivots = []
    row = 0
    for col in range(a.shape[1]):
        matches = np.flatnonzero(a[row:, col])
        if not len(matches):
            continue
        other = row + int(matches[0])
        a[[row, other]] = a[[other, row]]
        for r in np.flatnonzero(a[:, col]):
            if r != row:
                a[r] ^= a[row]
        pivots.append(col)
        row += 1
        if row == len(a):
            break
    return a, pivots


def gf2_rank(a):
    return len(gf2_rref(a)[1])


def gf2_nullspace(a):
    a, pivots = gf2_rref(a)
    free = [j for j in range(a.shape[1]) if j not in pivots]
    result = np.zeros((len(free), a.shape[1]), dtype=np.uint8)
    for i, j in enumerate(free):
        result[i, j] = 1
        for r, col in enumerate(pivots):
            result[i, col] = a[r, j]
    return result


def _inverse(a):
    n = len(a)
    reduced, pivots = gf2_rref(np.hstack((a, np.eye(n, dtype=np.uint8))))
    if pivots[:n] != list(range(n)):
        raise ValueError("Singular logical pairing matrix")
    return reduced[:, n:]


def _quotient_basis(kernel, stabilizers):
    reduced, pivots = gf2_rref(stabilizers)
    basis = reduced[:len(pivots)].copy()
    selected = []
    rank = len(pivots)
    for v in kernel:
        extended = np.vstack((basis, v))
        new_rank = gf2_rank(extended)
        if new_rank > rank:
            selected.append(v.copy())
            basis = extended
            rank = new_rank
    return np.array(selected, dtype=np.uint8).reshape((-1, kernel.shape[1]))


def css_logicals(hx, hz):
    if np.any((hx @ hz.T) & 1):
        raise ValueError("Noncommuting CSS checks")
    lx = _quotient_basis(gf2_nullspace(hz), hx)
    raw_lz = _quotient_basis(gf2_nullspace(hx), hz)
    pairing = (lx @ raw_lz.T) & 1
    lz = (_inverse(pairing).T @ raw_lz) & 1
    if not np.array_equal((lx @ lz.T) & 1, np.eye(len(lx), dtype=np.uint8)):
        raise AssertionError("Logical operators are not canonically paired")
    return lx, lz


def _translation_matrix(shape, offset):
    l, m = shape
    out = np.zeros((l*m, l*m), dtype=np.uint8)
    for i in range(l):
        for j in range(m):
            out[i*m+j, ((i+offset[0]) % l)*m + (j+offset[1]) % m] = 1
    return out


def _bicycle(name, shape, a_offsets, b_offsets, schedule, distance):
    a_parts = [_translation_matrix(shape, a) for a in a_offsets]
    b_parts = [_translation_matrix(shape, b) for b in b_offsets]
    a = np.bitwise_xor.reduce(a_parts)
    b = np.bitwise_xor.reduce(b_parts)
    hx = np.hstack((a, b))
    hz = np.hstack((b.T, a.T))
    m = len(a)
    x_neighbors = np.array([
        [int(np.flatnonzero(p[row])[0]) for p in a_parts]
        + [m + int(np.flatnonzero(p[row])[0]) for p in b_parts]
        for row in range(m)], dtype=int)
    z_neighbors = np.array([
        [int(np.flatnonzero(p[:, row])[0]) for p in b_parts]
        + [m + int(np.flatnonzero(p[:, row])[0]) for p in a_parts]
        for row in range(m)], dtype=int)
    return dict(name=name, n=2*m, hx=hx, hz=hz,
                x_neighbors=x_neighbors, z_neighbors=z_neighbors,
                schedule=schedule, static_distance=distance,
                group_shape=list(shape), a_offsets=[list(v) for v in a_offsets],
                b_offsets=[list(v) for v in b_offsets])


def _surface7():
    # Reuse Stim's standard hook-safe rotated-code CNOT ordering.  Native
    # RX/MX replaces its ideal H-R/H-M basis changes, matching IBM's primitives.
    template = stim.Circuit.generated("surface_code:rotated_memory_x", distance=7, rounds=1)
    instructions = list(template.flattened())
    data = [t.value for op in instructions if op.name == "RX" for t in op.targets_copy()]
    ancillas = [t.value for op in instructions if op.name == "R" for t in op.targets_copy()]
    x_anc = next([t.value for t in op.targets_copy()] for op in instructions if op.name == "H")
    z_anc = [q for q in ancillas if q not in x_anc]
    data_map = {q:i for i, q in enumerate(data)}
    x_map = {q:i for i, q in enumerate(x_anc)}
    z_map = {q:i for i, q in enumerate(z_anc)}
    hx = np.zeros((len(x_anc), len(data)), dtype=np.uint8)
    hz = np.zeros((len(z_anc), len(data)), dtype=np.uint8)
    layers = []
    for op in instructions:
        if op.name != "CX":
            continue
        ts = [t.value for t in op.targets_copy()]
        layer = []
        for a,b in zip(ts[::2], ts[1::2]):
            if a in x_map:
                check, datum = x_map[a], data_map[b]
                hx[check, datum] = 1
                layer.append(["X", check, datum])
            elif b in z_map:
                check, datum = z_map[b], data_map[a]
                hz[check, datum] = 1
                layer.append(["Z", check, datum])
            else:
                raise AssertionError("Unexpected surface-code CNOT")
        layers.append(layer)
    return dict(name="surface7", n=49, hx=hx, hz=hz,
                schedule={"layers":layers}, static_distance=7,
                construction="Stim rotated distance-7 CNOT ordering; native RX/MX")


@lru_cache(maxsize=None)
def _get_code_cached(name):
    if name == "w7":
        code = _bicycle(name, (2,14), [(0,0),(0,1),(1,10)],
                         [(0,0),(0,3),(1,11),(1,2)],
                         # Selected eight-layer candidate from schedule search.
                         # The original seven-layer candidate had exact circuit
                         # distance 5; its audit is archived in distance7layer.json.
                         {"x":[0,None,3,6,5,4,1,2], "z":[6,5,0,1,3,2,None,4]}, 7)
    elif name == "bb72":
        code = _bicycle(name, (6,6), [(3,0),(0,1),(0,2)],
                         [(0,3),(1,0),(2,0)],
                         {"x":[None,1,4,3,5,0,2], "z":[3,5,0,1,2,4,None]}, 6)
    elif name == "surface7":
        code = _surface7()
    else:
        raise ValueError(f"Unknown code {name!r}; use w7, bb72, surface7")
    code["lx"], code["lz"] = css_logicals(code["hx"], code["hz"])
    code["k"] = len(code["lx"])
    code["num_x_checks"] = len(code["hx"])
    code["num_z_checks"] = len(code["hz"])
    code["rank_x"] = gf2_rank(code["hx"])
    code["rank_z"] = gf2_rank(code["hz"])
    code["hardware_qubits"] = code["n"] + len(code["hx"]) + len(code["hz"])
    expected = {"w7":(56,12,112), "bb72":(72,12,144), "surface7":(49,1,97)}[name]
    assert (code["n"], code["k"], code["hardware_qubits"]) == expected
    return code


def get_code(name: str) -> dict[str, Any]:
    """Return an independent copy of matrices, logicals and schedule metadata."""
    return copy.deepcopy(_get_code_cached(name))


def schedule_layers(code, schedule=None):
    """Return layers of [sector, check_index, data_index] operations."""
    schedule = code["schedule"] if schedule is None else schedule
    if "layers" in schedule:
        layers = copy.deepcopy(schedule["layers"])
    else:
        xs, zs = schedule["x"], schedule["z"]
        if len(xs) != len(zs):
            raise ValueError("X and Z schedules must have equal layer counts")
        layers = []
        for xd, zd in zip(xs, zs):
            layer = []
            if xd is not None and xd != "idle":
                layer += [["X", i, int(q)] for i,q in enumerate(code["x_neighbors"][:, int(xd)])]
            if zd is not None and zd != "idle":
                layer += [["Z", i, int(q)] for i,q in enumerate(code["z_neighbors"][:, int(zd)])]
            layers.append(layer)
    validate_schedule(code, layers)
    return layers


def validate_schedule(code, layers):
    """Check edge coverage, hardware conflicts and CSS extraction parity."""
    n = code["n"]
    nx, nz = len(code["hx"]), len(code["hz"])
    tx = np.full((nx,n), -1, dtype=int)
    tz = np.full((nz,n), -1, dtype=int)
    for t,layer in enumerate(layers):
        used = set()
        for sector,check,datum in layer:
            if sector not in ("X", "Z"):
                raise ValueError("Schedule sector must be X or Z")
            anc = n+check if sector == "X" else n+nx+check
            if datum in used or anc in used:
                raise ValueError(f"Qubit collision at CNOT layer {t}")
            used.update((datum, anc))
            arr = tx if sector == "X" else tz
            if arr[check,datum] != -1:
                raise ValueError("Repeated check/data edge")
            arr[check,datum] = t
    if not np.array_equal(tx >= 0, code["hx"].astype(bool)):
        raise ValueError("X schedule does not cover the X check matrix exactly")
    if not np.array_equal(tz >= 0, code["hz"].astype(bool)):
        raise ValueError("Z schedule does not cover the Z check matrix exactly")
    for i in range(nx):
        for j in range(nz):
            overlap = (tx[i] >= 0) & (tz[j] >= 0)
            if int(np.sum(tx[i,overlap] < tz[j,overlap])) % 2:
                raise ValueError(f"Noncommuting extraction order: X check {i}, Z check {j}")


def _pauli_product(row, basis, length, ref=None):
    p = stim.PauliString(length)
    for q in np.flatnonzero(row):
        p[int(q)] = basis
    if ref is not None:
        p[ref] = basis
    return p


@lru_cache(maxsize=None)
def _bell_initialization(name):
    code = _get_code_cached(name)
    n,k = code["n"],code["k"]
    gens = []
    for matrix,basis in ((code["hx"],"X"),(code["hz"],"Z")):
        independent,pivots = gf2_rref(matrix)
        gens += [_pauli_product(row,basis,n+k) for row in independent[:len(pivots)]]
    for i in range(k):
        gens.append(_pauli_product(code["lx"][i],"X",n+k,n+i))
        gens.append(_pauli_product(code["lz"][i],"Z",n+k,n+i))
    prep = stim.Tableau.from_stabilizers(gens).to_circuit()
    physical = stim.Circuit()
    for op in prep:
        targets = [t.value if t.value < n else code["hardware_qubits"]+t.value-n
                   for t in op.targets_copy()]
        physical.append(op.name,targets,op.gate_args_copy())
    return physical


def build_circuit(code_name: str, p: float, rounds: int | None = None,
                  idle_scale: float = 0.0, schedule=None):
    """Return (circuit, JSON-serializable metadata).

    Noise: CNOT DEPOLARIZE2(p); preparation/readout flips p; data-only
    idle DEPOLARIZE1(idle_scale*p), following IBM's primitive operations.
    There are exactly ``rounds`` noisy cycles followed by two clean cycles.
    Z ancillas start perfectly initialized.  Each cycle resets Z ancillas
    at its end, so the last noisy preparation is seen by the first clean
    cycle, matching the supplied IBM Pauli-frame experiment boundary.
    """
    if rounds is None:
        rounds = _get_code_cached(code_name)["static_distance"]
    if rounds < 1 or not isinstance(rounds,int):
        raise ValueError("rounds must be a positive integer")
    if not 0 <= p <= 1 or not 0 <= idle_scale*p <= 1:
        raise ValueError("Invalid noise probability")
    code = _get_code_cached(code_name)
    chosen_schedule = copy.deepcopy(code["schedule"] if schedule is None else schedule)
    layers = schedule_layers(code,chosen_schedule)
    n,k,nx,nz = code["n"],code["k"],len(code["hx"]),len(code["hz"])
    xanc = list(range(n,n+nx))
    zanc = list(range(n+nx,n+nx+nz))
    c = _bell_initialization(code_name).copy()
    c.append("R",zanc)  # IBM starts with no initial Z-preparation fault.
    c.append("TICK")
    last_record = {}
    detector_groups = {"X":[], "Z":[]}
    measurement_count = 0
    hardware = code["hardware_qubits"]
    ibm_timing = (code_name == "bb72" and chosen_schedule == code["schedule"])

    def idle(active_data, probability):
        if probability and idle_scale:
            targets = [q for q in range(n) if q not in active_data]
            if targets:
                c.append("DEPOLARIZE1",targets,probability*idle_scale)

    def cx_layer(layer, probability):
        pairs = []
        active_data = set()
        for sector,check,datum in layer:
            pairs.extend((n+check,datum) if sector == "X" else (datum,n+nx+check))
            active_data.add(datum)
        if pairs:
            c.append("CX",pairs)
            if probability:
                c.append("DEPOLARIZE2",pairs,probability)
        idle(active_data,probability)

    def measure(sector, probability, round_index):
        nonlocal measurement_count
        targets = xanc if sector == "X" else zanc
        c.append("MX" if sector == "X" else "M",targets,probability)
        old_count = measurement_count
        measurement_count += len(targets)
        for i,_ in enumerate(targets):
            current = old_count+i
            key = (sector,i)
            terms = [stim.target_rec(current-measurement_count)]
            if key in last_record:
                terms.append(stim.target_rec(last_record[key]-measurement_count))
            detector_groups[sector].append(c.num_detectors)
            c.append("DETECTOR",terms,[i,0 if sector == "X" else 1,round_index])
            last_record[key] = current

    for r in range(rounds+2):
        probability = p if r < rounds else 0.0
        c.append("RX",xanc)
        if probability:
            c.append("Z_ERROR",xanc,probability)
        if ibm_timing:
            for t,layer in enumerate(layers):
                if t == 6:
                    measure("Z",probability,r)
                cx_layer(layer,probability)
                c.append("TICK")
            measure("X",probability,r)
            idle(set(),probability)
            c.append("R",zanc)
            if probability:
                c.append("X_ERROR",zanc,probability)
            c.append("TICK")
        else:
            idle(set(),probability)
            c.append("TICK")
            for layer in layers:
                cx_layer(layer,probability)
                c.append("TICK")
            measure("Z",probability,r)
            measure("X",probability,r)
            idle(set(),probability)
            c.append("TICK")
            # This reset occupies the same preparation slot as the next
            # cycle's RX.  No artificial extra idle time is introduced.
            c.append("R",zanc)
            if probability:
                c.append("X_ERROR",zanc,probability)

    # First k bits diagnose a residual logical X; next k diagnose logical Z.
    for sector,matrix,basis,start in (("X",code["lz"],"Z",0),
                                     ("Z",code["lx"],"X",k)):
        for i,row in enumerate(matrix):
            qubits = list(map(int,np.flatnonzero(row)))+[hardware+i]
            targets = []
            factory = stim.target_z if basis == "Z" else stim.target_x
            for q in qubits:
                if targets:
                    targets.append(stim.target_combiner())
                targets.append(factory(q))
            c.append("MPP",targets)
            c.append("OBSERVABLE_INCLUDE",[stim.target_rec(-1)],start+i)

    metadata = {
        "code":code_name,"n":n,"k":k,"static_distance":code["static_distance"],
        "rank_x":code["rank_x"],"rank_z":code["rank_z"],
        "num_x_checks":nx,"num_z_checks":nz,"hardware_qubits":hardware,
        "bookkeeping_reference_qubits":k,"stim_qubits":c.num_qubits,
        "rounds":rounds,"noisy_rounds":rounds,"noiseless_tail_rounds":2,
        "p":p,"idle_scale":idle_scale,"idle_noise_targets":"data only",
        "cnot_layers_per_cycle":len(layers),
        "time_slots_per_cycle":8 if ibm_timing else len(layers)+2,
        "cnots_per_cycle":sum(map(len,layers)),"schedule":chosen_schedule,
        "num_detectors":c.num_detectors,"num_observables":c.num_observables,
        "observable_order":[f"logical_X_{i}" for i in range(k)]+[f"logical_Z_{i}" for i in range(k)],
        "observable_groups":{"X":list(range(k)),"Z":list(range(k,2*k))},
        "detector_groups":detector_groups,
        "group_conventions":"observable X/Z means residual Pauli component; detector X/Z means measured check basis (opposite sensitivity)",
        "physical_qubits":{"data":list(range(n)),"X_ancillas":xanc,"Z_ancillas":zanc,"ideal_references":list(range(hardware,hardware+k))},
        "experiment":"joint Pauli-frame memory; ideal encoded Bell references",
        "noise_model":"IBM depolarizing primitives; preparation and measurement p; CNOT p; idle_scale*p on idle data",
        "boundary":"initial encoded Bell state and initial Z prep perfect; Z prep at previous cycle end; two clean cycles; ideal logical/ref MPP",
        "distance_status":"not certified by circuit construction; see distance audit output",
    }
    return c, metadata
