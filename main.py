# == Some notes ==
# Want to reproduce the results from https://arxiv.org/pdf/2411.05628 on using QC for (2+1)-QED

# External modules
from modules import *

# My modules
from classes import *
from circuit_helpers import *
from QED_hamiltonian import *

DYNAMICAL_LINKS = [((0,0),1)]
# self.n_dynamical_links = self.n_links - (self.n_fermion_qubits - 1)
# 2x2 : any one link
# 3x2 : any two links
# 3x3 : 

parameters = {
    'L_x': 2,
    'L_y': 2,
    'gauge_truncation': 1,
    'n_fermion_layers' : 1,
    'shots': 10000,
    'dynamical_links': DYNAMICAL_LINKS, # hardset atm - will generate some good candidates for different lattices
    'm': 1.0,
    'g': 1.0,
    'a': 1.0
}

print("RUNNING")

# Now need to build a nice VQE
circuit, observables, thetas, total_thetas = initiate_circuit_observables(parameters['L_x'],parameters['L_y'],parameters['n_fermion_layers'],parameters['gauge_truncation'])
thetas_values = [1.0]*total_thetas

print(circuit.draw())

hamiltonian = generate_qed_hamiltonian(parameters)

def thetas_only_wrapper(thetas_values):
    return qed_vqe(thetas_values, thetas, hamiltonian, circuit, observables, parameters['shots'])

print(thetas_only_wrapper(thetas_values))

    
