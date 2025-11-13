# == Some notes ==
# Want to reproduce the results from https://arxiv.org/pdf/2411.05628 on using QC for (2+1)-QED

# Other modules
from modules import *

# My modules
from classes import *
# from QED_hamiltonian import *
from gauge_helpers import *
from circuit_helpers import *

def mass_term_n(hamiltonian, lattice, n, m): 
    gauge_string = 'I'*lattice.n_gauge_qubits
    fermion_before = 'I'*n
    fermion_after = 'I'*(lattice.n_fermion_qubits - n - 1)

    coordinates = lattice.get_coordinates(n)
    term = gauge_string + fermion_before + 'Z' + fermion_after
    coeff = m*(-1)**(coordinates[0] + coordinates[1]) / 2

    hamiltonian.add_term(term, coeff)
    hamiltonian.add_term('I'*lattice.n_qubits, coeff)

    return hamiltonian

def electric_field_term_n_direction(hamiltonian, lattice, n, direction, g):
    index = lattice.labels[n]
    link_index = lattice.link_indexing[(index, direction)]
    coeff = g**2 / 2

    fermion_string = 'I' * lattice.n_fermion_qubits

    gauge_before = 'I' * (link_index * lattice.qubits_per_gauge)
    gauge_after = 'I' * (lattice.n_gauge_qubits - (link_index + 1) * lattice.qubits_per_gauge)

    Es = {
        2: {
            "IZ" : -0.5,
            "ZI" : -0.5
        },
        3: {
            "IIZ" : -1.5,
            "IZI" : -1.0,
            "ZII" : -0.5
        }
    }

    for key in list(Es[lattice.qubits_per_gauge]):
        term = gauge_before + key + gauge_after + fermion_string
        value = Es[lattice.qubits_per_gauge][key]
        hamiltonian.add_term(term, value*coeff)

    h_sq = copy.copy(hamiltonian)
    hamiltonian.hamiltinian = h_sq.multiply_hamiltonians(hamiltonian)
    return hamiltonian

def electric_field_term_n(hamiltonian, lattice, n, g):
    for direction in lattice.directions[n]:
        E_temp = Hamiltonian(lattice.n_qubits)
        E_temp = electric_field_term_n_direction(E_temp, lattice, n, direction, g)
        hamiltonian.add_hamiltonians(E_temp)
    return hamiltonian

def U_term_n(hamiltonian, lattice, n, direction, g, a): 
    index = lattice.labels[n]
    link_index = lattice.link_indexing[(index, direction)]
    coeff = -1/(2*(a**2)*(g**2))

    fermion_string = 'I' * lattice.n_fermion_qubits

    gauge_before = 'I' * (link_index * lattice.qubits_per_gauge)
    gauge_after = 'I' * (lattice.n_gauge_qubits - (link_index + 1) * lattice.qubits_per_gauge)

    Us = {
            2: {
        "XI": 1,
        "XX": 1,
        "YI": -1j,
        "YX": -1j
        },        
     
            3:  {
        "IIX": 0.5,
        "IXX": 0.25,
        "XXX": 0.25,
        "IYX": 0.25j, 
        "XYX": 0.25j,
        "YII": 0.5j,
        "YXI": 0.25j,
        "YXX": 0.25j,
        "YYI": 0.25,
        "YYX": -0.25
        }
        }
    

    for key in list(Us[lattice.qubits_per_gauge]):
        term = gauge_before + key + gauge_after + fermion_string
        value = Us[lattice.qubits_per_gauge][key]
        hamiltonian.add_term(term, value*coeff*(1/(2**(n-1))))
    
    return hamiltonian

def magnetic_term_n(hamiltonian, lattice, n, g, a):
    index = lattice.labels[n]
    ns_directions = [
        (lattice.reverse_labels[(index[0],index[1])], 1),
        (lattice.reverse_labels[(index[0]+1,index[1])], 2),
        (lattice.reverse_labels[(index[0],index[1]+1)], 1),
        (lattice.reverse_labels[(index[0],index[1])], 2)
    ]

    Us = [Hamiltonian(lattice.n_qubits) for _ in range(4)]

    for i in range(4): Us[i] = U_term_n(Us[i], lattice, ns_directions[i][0], ns_directions[i][1], g, a)
    Us[0].hamiltonian = Us[0].to_conjugate()
    Us[1].hamiltonian = Us[1].to_conjugate()

    P_n = Hamiltonian(lattice.n_qubits)
    P_n.add_term('I'*lattice.n_qubits, 1.0)

    for i in range(4): P_n.multiply_hamiltonians[Us[i]]

    P_n.to_conjugate()

    P_n_dagger = Hamiltonian(lattice.n_qubits)
    P_n_dagger.hamiltonian = P_n.conjugate

    hamiltonian.add_hamiltonians(P_n)
    hamiltonian.add_hamiltonians(P_n_dagger)

    return hamiltonian




def kinetic_term_n(hamiltonian, lattice, n, a, index):
    return hamiltonian
 

L_x = 2
L_y = 2
gauge_truncation = 1
m = 1.0
g = 1.0
a = 1.0

# Gauge qubits are ordered according to the lists in directions, i.e, site order, then 1 and/or 2 (dependant on site)

lattice = Lattice(L_x,L_y,gauge_truncation)
hamiltonian = Hamiltonian(lattice.n_qubits)

# Mass term - WORKS
for n in range(lattice.n_fermion_qubits):
    hamiltonian = mass_term_n(hamiltonian, lattice, n, m)

# Electric field term - WORKS
for n in range(lattice.n_links):
    hamiltonian = electric_field_term_n(hamiltonian, lattice, n, g)

# Magnetic field term - WORKS (only have space for n = 0 with L_x = L_y = 2)
for n in range(lattice.possible_plaquettes):
    magnetic_term_n(hamiltonian, lattice, n, g, a)
    
hamiltonian.cleanup() # Ends with loads of 0 + 0j otherwise

# Kinetic term
# This will just take a lot of multiplications
# Will also need to make a raising and lowering operator function maybe


    
    
