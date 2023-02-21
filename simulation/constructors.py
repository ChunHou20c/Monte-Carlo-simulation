"""This module contains all the constructors for the data structures for this project"""

from simulation import DBT
from data_structure.molecule import atom
from typing import Generator, Type

def _Create_Atom(atom_type:str, *args)-> atom:
    """this function return an atom"""
    
    func = {'atom': atom,
            'S': DBT.Sulphur,
            'C': DBT.Carbon,
            'H': DBT.Hydrogen,
            'N': DBT.Nitrogen,
            'O': DBT.Oxygen}
    
    return func[atom_type](*args)

def _gro_atom_line_parser(line:str)-> tuple[str, float, float, float]:
    """this function parse the single line in gromac file to obtain the string and the x, y, z coordinate of the atom"""
    
    data = line.split()
    #this line extract the element from the string, currently it is fine to take the first character only as this will be faster
    #and there is no element with 2 character currently
    element = data[1][0]
    
    #as the index of the atom will get mix up with the element at high number, it is ignored here
    x = float(data[-3])
    y = float(data[-2])
    z = float(data[-1])

    return element, x, y, z

def split(list_to_split:list,chunk_size:int)->Generator:
    """generator to split raw data"""

    for i in range(0,len(list_to_split), chunk_size):
        yield list_to_split[i:i+chunk_size]

def _Gro_file_parser(filepath:str, chunk_size:int, molecule:Type[DBT.DBT])->dict:
    """This function parse the .gro file and store the meta data in a dictionary"""

    with open(filepath, 'r') as f:
        
        title = next(f)
        timestamp = float(title.split()[4])

        num_of_atoms = next(f) 
        
        other_data = f.readlines()
        atomic_data, boundary_data = other_data[:-1], float(other_data[-1].split()[0])
        
        #parse the line and make it into list of atoms
        parsed_atomic_data = [_gro_atom_line_parser(i) for i in atomic_data]
        list_of_atoms = [_Create_Atom(i[0],*i[1:]) for i in parsed_atomic_data]

        list_of_groupped_atoms = [i for i in split(list_of_atoms, chunk_size)]
        
        list_of_molecules = [molecule(atoms, f'{index}DBT') for index,atoms in enumerate(list_of_groupped_atoms)]
        
        return {'title': title, 
                'timestamp': timestamp,
                'num_of_atoms': num_of_atoms,
                'boundary': boundary_data,
                'list_of_molecules': list_of_molecules}
