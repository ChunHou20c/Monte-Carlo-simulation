from data_structure.molecule import molecule, atom 
from data_structure.cube import Cube
from simulation.DBT1 import DBT1
from simulation import constructors
from simulation import algorithm
from itertools import combinations

def molecular_data_struct_test():
    atom1 = atom(1,1,1,12)
    print(atom1)

    molecule1 = molecule([atom1 for _ in range(2)])
    print(molecule1)

    dbt = DBT1([atom1 for _ in range(56)])
    print(dbt)

    sulphur = constructors._Create_Atom('S',1,1,1)
    print(sulphur)

    file_data = constructors._Gro_file_parser('Frames/DBT1-00', 56)
    print(file_data.keys())
    print(len(file_data['list_of_molecules']))

def cube_test():
    cube1 = Cube(1)
    print(cube1)
    cube2 = cube1.neighbor_cube(1).neighbor_cube(relative_y=-1)
    print(cube2)

file_data = constructors._Gro_file_parser('Frames/DBT1-00', 56)
molecules = file_data['list_of_molecules']
#do a completion on all the molecules
for index,molecule in enumerate(molecules):
    molecules[index] = algorithm.complete_molecule(molecule)

count = 0
for m1, m2 in combinations(molecules, 2):
    
    if (algorithm.molecular_pair_relation(m1, m2, 1.2, 9.88696) is not None):

        count = count + 1

print(f'number of pair found is {count}')
