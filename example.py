"""
An example of how to call the routines applied to a CrFeCoNi special quasirandom structure.

C. D. Woodgate, University of Bristol, 2025
"""

from asro_calculator import asro_conditional_probabilities, asro_warren_cowley
from ase.io import read
import numpy as np

# Elements is the list of elements in the cell
elements1 = ['Cr', 'Fe', 'Co', 'Ni']

# Distances is the set of lattice distances
distances1 = [3.57/np.sqrt(2), 3.57]

# Read the ASE Atoms object
atoms1 = read('structures/a1_crfeconi_supercell.xyz')

# Calculate the pair probabilities. 'tol' controls the tolerance on distance evaluation
P_pq_1 = asro_conditional_probabilities(atoms1, elements1, distances1, tol=1e-2)

# Calculate the Warren-Cowley parameters. 'tol' controls the tolerance on distance evaluation
wc_pq_1 = asro_warren_cowley(atoms1, elements1, distances1, tol=1e-2)

print('\nFor the example A1 CrFeCoNi SQS, the conditional pair probabilities on first and second coordination shells are:\n')
print(P_pq_1, '\n')
print('And the Warren-Cowley parameters are:\n')
print(wc_pq_1, '\n')

print('This indicates the considered supercell is fairly disordered, as expected.')
