from dataclasses import dataclass
import stim
from full_clifford_sim._layer_translate import to_z_basis_interaction_circuit, to_optimized_circuit
from full_clifford_sim._noise import GidneyNoiseModel


def set_unique(llst):
    flat = []
    for sublist in llst:
        if sublist:
            for item in sublist:
                if item is not None:
                    flat.append(item)
    return (set(flat))  # unique elements


def insert_circuit_errs(orig_circ : stim.Circuit, 
                        p: float,
                        valid: bool=1,
                        convert_nac: bool = False,
                        cirq = False,
                        ) -> stim.Circuit:


    if not valid:
        return orig_circ  

    active_qubits = []
    for lin in orig_circ:
        active_qubits.append([i.qubit_value for i in lin.targets_copy()])

    if cirq:
        model = GidneyNoiseModel.cirq_uniform_depolarizing(p)
        mc = model.noisy_circuit(orig_circ, system_qubits=set_unique(active_qubits))
        mc.append("TICK")
        return mc

    if convert_nac:
        orig_circ = to_z_basis_interaction_circuit(orig_circ)
        model = GidneyNoiseModel.neutralatom(p)
        mc = model.noisy_circuit(orig_circ, system_qubits=set_unique(active_qubits))
        return mc
    
    else:
        model = GidneyNoiseModel.uniform_depolarizing(p)
        mc = model.noisy_circuit(orig_circ, system_qubits=set_unique(active_qubits))
        mc.append("TICK")
        return mc
