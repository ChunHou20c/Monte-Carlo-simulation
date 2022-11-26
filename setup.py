from simulation import constructors
from simulation.DBT1 import DBT1_numerator_matrix
import numpy as np

def setup():

    file_data = constructors._Gro_file_parser('Frames/DBT1-00', 56)
    molecule = file_data['list_of_molecules'][0]
    numerator_matrix = DBT1_numerator_matrix(molecule)
    with open('./cache/DBT1.npy', 'wb') as f:

        np.save(f, numerator_matrix)

    print('numerator cache file setup completed!')

if __name__ == '__main__':

    setup()
