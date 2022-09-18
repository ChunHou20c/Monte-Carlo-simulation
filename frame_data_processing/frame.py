"""This module is used to process the frame object into data that can be loaded"""
import itertools
import numpy as np
from typing import Generator
from frame_data_processing import molecule
from frame_data_processing import gen_coulomb_matrix
from functools import lru_cache

class frame():
    """The frame class contains the frame number, all the molecules number and also the sides of the box of the frame"""
    def __init__(self, path:str, working_dir:str = '.', molecule_size:int=56, cut_off_distance:float = 1.2) -> None:
        """The constructor to read all the data from the text file, take path of the file as argument"""

        self.molecule_size = molecule_size
        self.distance_cut_off = cut_off_distance

        with open(path) as frame_data:

        #do some processing first to cleanup unused line
        #first line is the frame number

            first_line = next(frame_data)
            step_number = first_line.split()[-1]
            self.dir_name = f'{working_dir}/results/{step_number}.npy' #this last element is the step we want  
            #The step number will be used to  construct the folder name for the molecules and pairs
            
            self.number_of_line = next(frame_data) #this is the number of line to be read as molecular data

            molecular_data = frame_data.readlines()
            self.box_data = molecular_data[-1]

            molecular_data_to_process = molecular_data[:-1]
            list_of_molecule = self.molecule_split(molecular_data_to_process)
            self.molecules=[molecule.molecule(m) for m in list_of_molecule]

    def molecule_split(self, raw_data)->list[list[str]]:
        """This function splits the aggregated molecular raw data into single molecule
        as for now the number of atoms in a molecule is known and fixed"""

        def split(list_to_split:list[str],chunk_size:int)->Generator:
            """generator to split the data"""

            for i in range(0,len(list_to_split), chunk_size):
                yield list_to_split[i:i+chunk_size]

        return list(split(raw_data,self.molecule_size))

    def gen_cm_list(self, numerator_matrix:np.ndarray)->list[np.ndarray]:
        """This function generate the list of matrix template for the coulomb matrix generation"""

        CM_list = []
        matrix_size = self.molecule_size*2
        for molecule1, molecule2 in itertools.combinations(self.molecules,2):
        
            molecular_distance = molecule.distance(molecule1,molecule2)
            if (molecular_distance<= self.distance_cut_off):

                x_list = [i.x for i in molecule1.atoms] + [i.x for i in molecule2.atoms]
                y_list = [i.y for i in molecule1.atoms] + [i.y for i in molecule2.atoms]
                z_list = [i.z for i in molecule1.atoms] + [i.z for i in molecule2.atoms]


                x = np.tile(x_list,(matrix_size,1))
                y = np.tile(y_list,(matrix_size,1))
                z = np.tile(z_list,(matrix_size,1))
                #generate coulomb matrix here
                
                cm = gen_coulomb_matrix.CoulombMatrix(x,y,z,numerator_matrix)
                
                matrix = cm.gen_CM()
                CM_list.append(matrix)
                
        return CM_list

    def print_attributes(self):
        """This method can be used to check the attributes of the data"""
        
        print(f'directory to save = {self.dir_name}')
        print(f'number of lines to save = {self.number_of_line}')
        print(f'box data = {self.box_data}')
        print(f'no of atoms in a single molecule = {self.molecule_size}')
    
@lru_cache(maxsize=1)    
def gen_numerator_matrix(frame:frame)->np.ndarray:
    """this method is used to get product of charge 1 and charge 2, as all the molecules are the same,
    there is only one numerator matrix throughout the program"""
    element_list=[i.element for i in frame.molecules[0].atoms]
    element_matrix = np.tile(element_list*2,(frame.molecule_size*2,1))
    charge_matrix = gen_coulomb_matrix.charge_vect(element_matrix, element_matrix.T, False)
    diagonal_elements = gen_coulomb_matrix.charge_vect(np.array(element_list*2),np.array(element_list*2), True)
    np.fill_diagonal(charge_matrix,diagonal_elements)

    return charge_matrix
