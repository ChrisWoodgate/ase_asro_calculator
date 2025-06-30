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
