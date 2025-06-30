from asro_calculator import asro_conditional_probabilities, asro_warren_cowley
from ase.io import read
import numpy as np

#================#
# Segregated fcc #
#================#

elements1 = ['Fe', 'Ni']
distances1 = [3.57/np.sqrt(2), 3.57]

atoms1 = read('structures/fcc_seg_feni.xyz')

P_pq_1 = asro_conditional_probabilities(atoms1, elements1, distances1)
wc_pq_1 = asro_warren_cowley(atoms1, elements1, distances1)

print('\nFor the test segregated FeNi supercell, the conditional pair probabilities on first and second coordination shells are:\n')
print(P_pq_1, '\n')
print('And the Warren-Cowley parameters are:\n')
print(wc_pq_1, '\n')

#================#
# Segregated bcc #
#================#

elements2 = ['Cu', 'Zn']
distances2 = [3.0*np.sqrt(3)*0.5, 3.0]

atoms2 = read('structures/bcc_seg_cuzn.xyz')

P_pq_2 = asro_conditional_probabilities(atoms2, elements2, distances2)
wc_pq_2 = asro_warren_cowley(atoms2, elements2, distances2)

print('\nFor the test segregated CuZn supercell, the conditional pair probabilities on first and second coordination shells are:\n')
print(P_pq_2, '\n')
print('And the Warren-Cowley parameters are:\n')
print(wc_pq_2, '\n')

