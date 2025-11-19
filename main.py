# == Some notes ==
# Want to reproduce the results from https://arxiv.org/pdf/2411.05628 on using QC for (2+1)-QED

# External modules
from modules import *

# My modules
from classes import *
from circuit_helpers import *
from QED_hamiltonian import *

parameters = {
    'L_x': 3,
    'L_y': 2,
    'gauge_truncation': 1,
    'n_fermion_layers' : 3,
    'shots': 10000,
    'dynamical_links': DYNAMICAL_LINKS, # hardset atm - will generate some good candidates for different lattices
    'm': 1.0,
    'g': 1.0,
    'a': 1.0,
    'max_iters': 1000,
}

print("RUNNING")

# Now need to build a nice VQE
circuit, observables, thetas, total_thetas, n_qubits = initiate_circuit_observables(parameters)
thetas_values = [np.random.uniform(0,1)]*total_thetas

hamiltonian = generate_qed_hamiltonian(parameters)
hamiltonian.latex_plot(save=True)





'''
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

print(get_state_counts(mini.x, thetas, circuit, observables, n_qubits, parameters['shots']))

    
'''