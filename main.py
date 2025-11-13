# == Some notes ==
# Want to reproduce the results from https://arxiv.org/pdf/2411.05628 on using QC for (2+1)-QED

# External modules
from modules import *

# My modules
from classes import *
from circuit_helpers import *
from QED_hamiltonian import *

parameters = {
    'L_x': 2,
    'L_y': 2,
    'gauge_truncation': 1,
    'm': 1.0,
    'g': 1.0,
    'a': 1.0
}


# Believe I have the right Hamiltonian now
# Need to now check their simulation parameters and set up my own



hamiltonian = generate_qed_hamiltonian(parameters)
hamiltonian.latex_print()

    
    
