# == Some notes ==
# Want to reproduce the results from https://arxiv.org/pdf/2411.05628 on using QC for (2+1)-QED

# External modules
from modules import *

# My modules
from classes import *
from circuit_helpers import *
from plot_helpers import *
from QED_hamiltonian import *

parameters = {
    'L_x': 3,
    'L_y': 2,
    'gauge_truncation': 1,
    'n_fermion_layers' : 1,
    'shots': 1e6,
    'dynamical_links': DYNAMICAL_LINKS, # hardset atm - will generate some good candidates for different lattices
    'm': 1.0,
    'g': 1.0,
    'a': 1.0,
    'max_iters': 1000,
}

print("RUNNING")


def run_and_print(g):
    parameters['g'] = g

    circuit, observables, thetas, total_thetas, n_qubits = initiate_circuit_observables(parameters)
    thetas_values = [np.random.uniform(0,1)]*total_thetas

    hamiltonian = generate_qed_hamiltonian(parameters)

    def thetas_only_wrapper(thetas_values):
        cost = qed_vqe(thetas_values, thetas, hamiltonian, circuit, observables, parameters['shots'])
        print(cost)
        return cost

    mini = scipy.optimize.minimize(thetas_only_wrapper, thetas_values, method = "COBYLA")
    print(mini)

    def get_state_counts(thetas_values, thetas, circuit, observables, n_qubits, shots):
        param_dict = dict(zip(thetas, thetas_values))
        circuit_values = circuit.assign_parameters(param_dict)

        return observables.full_z(circuit_values, n_qubits, shots)

    counts = get_state_counts(mini.x, thetas, circuit, observables, n_qubits, 512)
    print("Counts with shots = 512: ")
    print(counts)
    param_dict = dict(zip(thetas, mini.x))
    circuit_values = circuit.assign_parameters(param_dict)
    energy = observables.energy(circuit_values, hamiltonian, parameters['shots'])
    print(f"Energy at g = {g}: {energy}")
    return energy

g_values = np.linspace(0.001,3,10)
for g in g_values:
    print('-'*10)
    _ = run_and_print(g)


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Plan for what I need to calculate:
# - Groundstate histogram at different g's so that they match the KJ paper
# - Static potential at different values for g
# -- Here, we keep the `distance' constant, while increasing g seperates the sites (by inversally changing a)`
# - Flux tubes
#
# - Then want to introduce BEF

# Also need to introduce the fermion-gauge entanglement
# May be good to swap optimizers, just using COBYLA out of habit