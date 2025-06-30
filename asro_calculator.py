'''
Routines for calculation of the Warren-Cowley ASRO parameters for a simulation cell described by an ASE Atoms object.

C. D. Woodgate, University of Bristol, 2025
'''

def get_concentrations(atoms, elements):
    '''
    Return the concentrations of elements in a simulation cell described by an atoms object

            Parameters:
                    atoms (ASE atoms object): The atoms to consider
                    elements (list of str): The elements for which to look

            Returns:
                    concentrations (list of floats): List of concentrations of each element in elements
    '''
    
    # Get the labels of the atoms in the Atoms object
    labels = atoms.get_chemical_symbols()
    
    # Define an array of zeros for calculation of the concentrations
    concentrations = np.zeros(len(elements))
    
    # Loop over the labels of atoms in the cell
    for label in labels:
        i = elements.index(label)
        concentrations[i] += 1
    
    # Divide by the number of atoms in the cell
    concentrations = concentrations/float(len(atoms))
    
    return concentrations

def get_coordination_numbers(atoms, shell_distances, tol=1e-3):
    '''
    Return the coordination numbers of the atoms object, where coordination shell distances are specified by shell_distances

            Parameters:
                    atoms (ASE atoms object): The atoms to consider
                    shell_distances (list of floats): The distances for which to look

            Returns:
                    total_pairs (list of floats): List of numbers of atoms on each shell
    '''
    
    # Get the distances between all atoms in the simulation cell, applying the minimum image convention
    atoms_distances = atoms.get_all_distances(mic=True, vector=False)
    
    # Empty array for storing result
    total_pairs = np.zeros(len(shell_distances))
    
    # Loop over all pairs and check if the distance matches
    for i in range(len(atoms)):
        for j in range(len(atoms)):
            for k, dist in enumerate(shell_distances):
                if np.fabs(dist - atoms_distances[i,j]) < tol:
                    total_pairs[k] += 1
    
    # Divide by the number of atoms
    total_pairs = total_pairs/len(atoms)
    
    return total_pairs
