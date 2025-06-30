from asro_calculator import asro_conditional_probabilities, asro_warren_cowley
from ase.io import read
import numpy as np

#=============#
# A1 CrFeCoNi #
#=============#

elements1 = ['Cr', 'Fe', 'Co', 'Ni']
distances1 = [3.57/np.sqrt(2), 3.57]

atoms1 = read('structures/a1_crfeconi_supercell.xyz')

P_pq_1 = asro_conditional_probabilities(atoms1, elements1, distances1)
wc_pq_1 = asro_warren_cowley(atoms1, elements1, distances1)

print('\nFor the test A1 CrFeCoNi SQS, the conditional pair probabilities on first and second coordination shells are:\n')
print(P_pq_1, '\n')
print('And the Warren-Cowley parameters are:\n')
print(wc_pq_1, '\n')

#=============#
# A1 Fe60Ni40 #
#=============#

elements2 = ['Fe', 'Ni']
distances2 = [3.57/np.sqrt(2), 3.57]

atoms2 = read('structures/a1_fe60ni40_supercell.xyz')

P_pq_2 = asro_conditional_probabilities(atoms2, elements2, distances2)
wc_pq_2 = asro_warren_cowley(atoms2, elements2, distances2)

print('\nFor the test A1 Fe60Ni40 SQS, the conditional pair probabilities on first and second coordination shells are:\n')
print(P_pq_2, '\n')
print('And the Warren-Cowley parameters are:\n')
print(wc_pq_2, '\n')

#=============#
# A2 NbMoTaW  #
#=============#

elements3 = ['Nb', 'Mo', 'Ta', 'W']
distances3 = [3.3*np.sqrt(3)*0.5, 3.3]

atoms3 = read('structures/a2_nbmotaw_supercell.xyz')

P_pq_3 = asro_conditional_probabilities(atoms3, elements3, distances3)
wc_pq_3 = asro_warren_cowley(atoms3, elements3, distances3)

print('\nFor the test A2 NbMoTaW SQS, the conditional pair probabilities on first and second coordination shells are:\n')
print(P_pq_3, '\n')
print('And the Warren-Cowley parameters are:\n')
print(wc_pq_3, '\n')

#=============#
# A2 Fe3Ga    #
#=============#

elements4 = ['Fe', 'Ga']
distances4 = [3.0*np.sqrt(3)*0.5, 3.0]

atoms4 = read('structures/a2_fe3ga_supercell.xyz')

P_pq_4 = asro_conditional_probabilities(atoms4, elements4, distances4)
wc_pq_4 = asro_warren_cowley(atoms4, elements4, distances4)

print('\nFor the test A2 Fe3Ga SQS, the conditional pair probabilities on first and second coordination shells are:\n')
print(P_pq_4, '\n')
print('And the Warren-Cowley parameters are:\n')
print(wc_pq_4, '\n')

