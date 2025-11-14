from classes import *

def mass_term_n(hamiltonian, lattice, n, m): 
    gauge_string = 'I'*lattice.n_dynamical_gauge_qubits
    fermion_before = 'I'*n
    fermion_after = 'I'*(lattice.n_fermion_qubits - n - 1)

    coordinates = lattice.get_coordinates(n)
    term = gauge_string + fermion_before + 'Z' + fermion_after
    coeff = m*(-1)**(coordinates[0] + coordinates[1]) / 2

    hamiltonian.add_term(term, coeff)
    hamiltonian.add_term('I'*lattice.n_qubits, coeff)

    return hamiltonian

def electric_field_term_n_direction(hamiltonian, lattice, n, direction, g, dynamical_links):
    index = lattice.labels[n]
    link_index = lattice.link_indexing[(index, direction)]
    coeff = g**2 / 2
    if (index, direction) in dynamical_links:
        fermion_string = 'I' * lattice.n_fermion_qubits

        gauge_before = 'I' * ((link_index-1) * lattice.qubits_per_gauge)
        gauge_after = 'I' * (lattice.n_dynamical_gauge_qubits - (link_index + 1) * lattice.qubits_per_gauge)

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
        hamiltonian.multiply_hamiltonians(hamiltonian)
    else:
        hamiltonian.add_term('I'*lattice.n_qubits, 0)

    return hamiltonian

def electric_field_term_n(hamiltonian, lattice, n, g, dynamical_links):
    for direction in lattice.directions[n]:
        E_temp = Hamiltonian(lattice.n_qubits)
        E_temp = electric_field_term_n_direction(E_temp, lattice, n, direction, g, dynamical_links)
        hamiltonian.add_hamiltonians(E_temp)
    return hamiltonian

def U_term_n(hamiltonian, lattice, n, direction, dynamical_links): 
    index = lattice.labels[n]
    if (index,direction) in dynamical_links:
        link_index = lattice.link_indexing[(index, direction)]
        
        fermion_string = 'I' * lattice.n_fermion_qubits

        gauge_before = 'I' * (link_index * lattice.qubits_per_gauge)
        gauge_after = 'I' * (lattice.n_dynamical_gauge_qubits - (link_index + 1) * lattice.qubits_per_gauge)

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
            hamiltonian.add_term(term, value*(1/(2**(n-1))))  
        
    else:
        hamiltonian.add_term('I'*lattice.n_qubits, 1.0)
    return hamiltonian

def magnetic_term_n(hamiltonian, lattice, n, g, a, dynamical_links):
    index = lattice.labels[n]
    ns_directions = [
        (lattice.reverse_labels[(index[0],index[1])], 1),
        (lattice.reverse_labels[(index[0]+1,index[1])], 2),
        (lattice.reverse_labels[(index[0],index[1]+1)], 1),
        (lattice.reverse_labels[(index[0],index[1])], 2)
    ]

    Us = [Hamiltonian(lattice.n_qubits) for _ in range(4)]
    coeff = -1/(2*(a**2)*(g**2))
    for i in range(4): Us[i] = U_term_n(Us[i], lattice, ns_directions[i][0], ns_directions[i][1], dynamical_links)
    Us[0].hamiltonian = Us[0].to_conjugate()
    Us[1].hamiltonian = Us[1].to_conjugate()

    P_n = Hamiltonian(lattice.n_qubits)
    P_n.add_term('I'*lattice.n_qubits, 1.0)

    for i in range(4): P_n.multiply_hamiltonians(Us[i])

    P_n.to_conjugate()

    P_n_dagger = Hamiltonian(lattice.n_qubits)
    P_n_dagger.hamiltonian = P_n.conjugate

    hamiltonian.add_hamiltonians(P_n)
    hamiltonian.add_hamiltonians(P_n_dagger)

    hamiltonian.hamiltonian = multiply_hamiltonian_by_constant(hamiltonian.hamiltonian, coeff)

    return hamiltonian

def creation_operator_n_dict(hamiltonian, lattice, n):
    gauge_string = 'I'*lattice.n_dynamical_gauge_qubits
    fermions_before = 'Z'*(n)
    fermions_after = 'I'*(lattice.n_fermion_qubits - (n+1))
    hamiltonian.add_term(gauge_string + fermions_before + 'X' + fermions_after, (1j)**(n-1))
    hamiltonian.add_term(gauge_string + fermions_before + 'Y' + fermions_after, (-1j)*(1j)**(n-1))
    return hamiltonian

def annihilation_operator_n(hamiltonian, lattice, n):
    gauge_string = 'I'*lattice.n_dynamical_gauge_qubits
    fermions_before = 'Z'*(n)
    fermions_after = 'I'*(lattice.n_fermion_qubits - (n+1))
    hamiltonian.add_term(gauge_string + fermions_before + 'X' + fermions_after, (-1j)**(n-1))
    hamiltonian.add_term(gauge_string + fermions_before + 'Y' + fermions_after, (1j)*(1j)**(n-1))
    return hamiltonian

def kinetic_subterm_n(lattice, n, direction, dynamical_links):
    indices = lattice.labels[n]

    if direction == 1:
        indices_mu = (indices[0]+1, indices[1])
    elif direction == 2:
        indices_mu = (indices[0], indices[1]+1)

    n_mu = lattice.reverse_labels[indices_mu]

    hamiltonian = Hamiltonian(lattice.n_qubits)
    hamiltonian.add_term('I'*lattice.n_qubits, 1.0)

    creation_hamiltonian = Hamiltonian(lattice.n_qubits)
    U_hamiltonian = Hamiltonian(lattice.n_qubits)
    annihilation_hamiltonian = Hamiltonian(lattice.n_qubits)

    creation_hamiltonian = creation_operator_n_dict(creation_hamiltonian, lattice, n)
    U_hamiltonian = U_term_n(U_hamiltonian, lattice, n, direction, dynamical_links)
    annihilation_hamiltonian = annihilation_operator_n(annihilation_hamiltonian, lattice, n_mu)

    hamiltonian.multiply_hamiltonians(creation_hamiltonian)
    hamiltonian.multiply_hamiltonians(U_hamiltonian)
    hamiltonian.multiply_hamiltonians(annihilation_hamiltonian)


    return hamiltonian

def kinetic_term_n(hamiltonian, lattice, n, a, dynamical_links):
    # This works out from both directions for a given n
    temp_hamiltonian = Hamiltonian(lattice.n_qubits)
    temp_hamiltonian_dagger = Hamiltonian(lattice.n_qubits)
    indices = lattice.labels[n]
    coeffs = {
        1: 1j/(2*a),
        2: -1/(2*a)*((-1)**(indices[0] + indices[1]))
    }

    for direction in lattice.directions[n]:
        temp_hamiltonian = kinetic_subterm_n(lattice, n, direction, dynamical_links)
        temp_hamiltonian.to_conjugate()
        temp_hamiltonian_dagger.hamiltonian = temp_hamiltonian.conjugate

        if direction == 1: temp_hamiltonian_dagger.hamiltonian = multiply_hamiltonian_by_constant(temp_hamiltonian_dagger.hamiltonian, -1)

        temp_hamiltonian.add_hamiltonians(temp_hamiltonian_dagger)
        temp_hamiltonian.hamiltonian = multiply_hamiltonian_by_constant(temp_hamiltonian.hamiltonian, coeffs[direction])

        hamiltonian.add_hamiltonians(temp_hamiltonian)
    return hamiltonian
 
def generate_qed_hamiltonian(parameters):
    lattice = Lattice(parameters['L_x'],parameters['L_y'],parameters['gauge_truncation'])
    hamiltonian = Hamiltonian(lattice.n_qubits)
    
    # Mass term
    for n in range(lattice.n_fermion_qubits):
        hamiltonian = mass_term_n(hamiltonian, lattice, n, parameters['m'])

    # Electric field term
    for n in range(lattice.n_links):
        hamiltonian = electric_field_term_n(hamiltonian, lattice, n, parameters['g'], parameters['dynamical_links'])

    # Magnetic field term
    for n in lattice.plaquettes:
        hamiltonian = magnetic_term_n(hamiltonian, lattice, n, parameters['g'], parameters['a'], parameters['dynamical_links'])
    
    # Kinetic term
    for n in range(lattice.n_fermion_qubits):
        hamiltonian = kinetic_term_n(hamiltonian, lattice, n, parameters['a'], parameters['dynamical_links'])
        
    hamiltonian.cleanup() # Ends with loads of 0 + 0j otherwise

    return hamiltonian

