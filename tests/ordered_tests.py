from asro_calculator import asro_conditional_probabilities, asro_warren_cowley
from ase.io import read
import numpy as np

#=============#
# L10 Feni    #
#=============#

elements1 = ['Fe', 'Ni']
distances1 = [3.57/np.sqrt(2), 3.57]

atoms1 = read('structures/fcc_l10_feni.xyz')

P_pq_1 = asro_conditional_probabilities(atoms1, elements1, distances1)
wc_pq_1 = asro_warren_cowley(atoms1, elements1, distances1)

print('\nFor the test FCC L10 FeNi structure, the conditional pair probabilities on first and second coordination shells are:\n')
print(P_pq_1, '\n')
print('And the Warren-Cowley parameters are:\n')
print(wc_pq_1, '\n')

#=============#
# L12 Ni3Al   #
#=============#

elements2 = ['Ni', 'Al']
distances2 = [3.57/np.sqrt(2), 3.57]

atoms2 = read('structures/fcc_l12_ni3al.xyz')

P_pq_2 = asro_conditional_probabilities(atoms2, elements2, distances2)
wc_pq_2 = asro_warren_cowley(atoms2, elements2, distances2)

print('\nFor the test L12 Ni3Al structure, the conditional pair probabilities on first and second coordination shells are:\n')
print(P_pq_2, '\n')
print('And the Warren-Cowley parameters are:\n')
print(wc_pq_2, '\n')

#=============#
# B2 CuZn     #
#=============#

elements3 = ['Cu', 'Zn']
distances3 = [3.0*np.sqrt(3)*0.5, 3.0]

atoms3 = read('structures/bcc_b2_cuzn.xyz')

P_pq_3 = asro_conditional_probabilities(atoms3, elements3, distances3)
wc_pq_3 = asro_warren_cowley(atoms3, elements3, distances3)

print('\nFor the test B2 CuZn structure, the conditional pair probabilities on first and second coordination shells are:\n')
print(P_pq_3, '\n')
print('And the Warren-Cowley parameters are:\n')
print(wc_pq_3, '\n')

#=============#
# D03 Fe3Ga   #
#=============#

elements4 = ['Fe', 'Ga']
distances4 = [3.0*np.sqrt(3)*0.5, 3.0, 3.0*np.sqrt(2)]

atoms4 = read('structures/bcc_d03_fe3ga.xyz')

P_pq_4 = asro_conditional_probabilities(atoms4, elements4, distances4)
wc_pq_4 = asro_warren_cowley(atoms4, elements4, distances4)

print('\nFor the test D03 Fe3Ga structure, the conditional pair probabilities on first, second, and third coordination shells are:\n')
print(P_pq_4, '\n')
print('And the Warren-Cowley parameters are:\n')
print(wc_pq_4, '\n')

#=============#
# L21 Cu2MnAl #
#=============#

elements4 = ['Cu', 'Mn', 'Al']
distances4 = [3.0*np.sqrt(3)*0.5, 3.0, 3.0*np.sqrt(2)]

atoms4 = read('structures/bcc_l21_cu2mnal.xyz')

P_pq_4 = asro_conditional_probabilities(atoms4, elements4, distances4)
wc_pq_4 = asro_warren_cowley(atoms4, elements4, distances4)

print('\nFor the test L21 Cu2MnAl structure, the conditional pair probabilities on first, second, and third coordination shells are:\n')
print(P_pq_4, '\n')
print('And the Warren-Cowley parameters are:\n')
print(wc_pq_4, '\n')

