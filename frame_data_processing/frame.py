"""This module is used to process the frame object into data that can be loaded"""
import itertools
import numpy as np
from typing import Callable, Generator
from frame_data_processing import molecule
from frame_data_processing import gen_coulomb_matrix
from frame_data_processing import graph
from functools import lru_cache

class frame():
    """The frame class contains the frame number, all the molecules number and also the sides of the box of the frame"""

    def __init__(self, path:str,
            working_dir:str = '.',
            molecule_size:int=56,
            cut_off_distance:float = 1.2,
            Use_full_CM:bool = True) -> None:
        """The constructor to read all the data from the text file, take path of the file as argument"""
        
        self.first_line = ''
        self.second_line = ''
        self.last_line = ''
        self.working_dir = working_dir
        self.molecule_size = molecule_size
        self.cut_off_distance = cut_off_distance
        self.Use_full_CM = Use_full_CM
        self.CM_save_path = ''
        self.no_of_atoms = 0
        self.x_boundary = 0
        self.y_boundary = 0
        self.z_boundary = 0
        self.molecules = graph.Graph()
        self.cut_by_boundary_keys = []

        self.import_data(path)

    def build_graph(self, molecules: list[molecule.molecule]) ->None:
        """helper method to build the graph
        the molecules that are cut by the boundary should be excluded first"""
            
        #here the molecules are separated into 2 sets, 1 is the normal molecule and another is the molecules that is cut by border
        #these molecules will be treated separately by eventually will be added into the same graph, with different variation of a vertex

        m_cut_by_border = [m for m in molecules if self.is_cut_by_boundary(m)]
    
        #with open('test_files/broken_molecule.txt','w') as f:
        #    for m in m_cut_by_border:

        #        print(m.get_name())
        #        
        #        f.write(m.get_name())
        #        f.write('\n')

        m_not_cut_by_border = [m for m in molecules if m not in m_cut_by_border]
        
        for m in m_cut_by_border:

            m.check_condition()

            m.atom_processing()

            self.molecules.add_vertex(m)

            self.cut_by_boundary_keys.append(m.get_name())
        
        for m1, m2 in itertools.combinations(m_not_cut_by_border, 2):
            
            if (self.molecule_pair_is_close(m1, m2)):
                
                dis = molecule.distance(m1, m2)
                self.molecules.add_edge(m1, m2, dis)

    def import_data(self, path:str):
        """helper function to import the data from file.
        the molecule graph will be build here and update to the object"""

        with open(path) as frame_data:

        #do some processing first to cleanup unused line
        #first line is the frame number

            self.first_line = next(frame_data)
            step_number = self.first_line.split()[-1]
            self.CM_save_path = f'{self.working_dir}/results/{step_number}.npy' #this last element is the step we want  
            #The step number will be used to  construct the folder name for the molecules and pairs
            
            self.second_line = next(frame_data) #this is the number of line to be read as molecular data also the number of atoms in the frame

            molecular_data = frame_data.readlines()
            self.last_line = molecular_data[-1]
            box_data = molecular_data[-1].split()
            
            self.x_boundary = float(box_data[0])
            self.y_boundary = float(box_data[1])
            self.z_boundary = float(box_data[2])
            molecular_data_to_process = molecular_data[:-1]
            
            #list of molecules are list of string that is extracted from the text file
            list_of_molecule = self.molecule_split(molecular_data_to_process)
            molecules=[molecule.molecule(m) for m in list_of_molecule]

            self.build_graph(molecules)

    def molecule_split(self, raw_data: list[str])->list[list[str]]:
        """This function splits the aggregated molecular raw data into single molecule
        as for now the number of atoms in a molecule is known and fixed

        raw_data: list of string

        return: list of list of string"""

        def split(list_to_split:list[str],chunk_size:int)->Generator:
            """generator to split the data"""

            for i in range(0,len(list_to_split), chunk_size):
                yield list_to_split[i:i+chunk_size]

        return list(split(raw_data,self.molecule_size))

    def molecule_pair_is_close(self, m1: molecule.molecule, m2: molecule.molecule)->bool:
        """This function tell whether the two molecule is close enough to be considered as pair
            
            return: bool"""
        
        molecular_distance = molecule.distance(m1, m2)

        if (molecular_distance<= self.cut_off_distance):

            return True

        return False

    def gen_cm_list(self, numerator_matrix:np.ndarray, func:Callable)->list[np.ndarray]:
        """This function generate the list of matrix template for the coulomb matrix generation
        func is a function for filtering molecule
        use self.is_close_to_border of self.is_cut_by_boundary"""
        
        #this method will be changed later to handle the border case molecules

        molecules_to_compare = [m for m in self.molecules if func(m) == False] 

        CM_list = []
        for molecule1, molecule2 in itertools.combinations(molecules_to_compare,2):
            
            if (self.molecule_pair_is_close(molecule1, molecule2)):

                x1 = [i.x for i in molecule1.atoms]
                x2 = [i.x for i in molecule2.atoms]
                y1 = [i.y for i in molecule1.atoms]
                y2 = [i.y for i in molecule2.atoms]
                z1 = [i.z for i in molecule1.atoms]
                z2 = [i.z for i in molecule2.atoms]

                #generate coulomb matrix here
                cm = gen_coulomb_matrix.CoulombMatrix(x1, x2,
                        y1, y2,
                        z1, z2,
                        numerator_matrix,
                        self.molecule_size,
                        self.Use_full_CM)
                
                matrix = cm.gen_CM()
                CM_list.append(matrix)
                
        return CM_list

    def print_attributes(self):
        """This method can be used to check the attributes of the data"""
        
        print(f'directory to save = {self.CM_save_path}')
        print(f'number of lines to save = {self.no_of_atoms}')
        print(f'box data = {self.x_boundary}, {self.y_boundary}, {self.z_boundary}')
        print(f'no of atoms in a single molecule = {self.molecule_size}')

    def is_close_to_boundary(self, molecule:molecule.molecule)->bool:
        """This function check if the molecule is close to the bondary so that it can be excluded"""

        if (molecule.atoms[molecule.N_index].x < 1 or molecule.atoms[molecule.N_index].x > self.x_boundary - 1):
            return True

        if (molecule.atoms[molecule.N_index].y < 1 or molecule.atoms[molecule.N_index].y > self.y_boundary - 1):
            return True

        if (molecule.atoms[molecule.N_index].z < 1 or molecule.atoms[molecule.N_index].z > self.z_boundary - 1):
            return True

        return False

    def is_cut_by_boundary(self, molecule:molecule.molecule)->bool:
        """This method check if the molecule is cut in half by the boundary condition
        
        this method works by first getting the list of xyz coordinate from the molecule,
        which is obtained from get_xyz_list method in the molecule class

        if the range of the data is large (max - min more than limit then the molecule must be cut in border)"""

        limit = 5

        x, y, z = molecule.get_xyz_list()

        if (max(x) - min(x) > limit):

            return True
        
        if (max(y) - min(y) > limit):

            return True

        if (max(z) - min(z) > limit):

            return True

        return False

    def export_molecule(self, key:str, file:str)->None:
        """this method export the gro file with the selected molecule"""

        molecule = self.molecules.get_vertex(key).molecule

        with open(file, 'w') as f:
            
            f.write(self.first_line)
            f.write('{:>5}\n'.format(56))
            
            for line in molecule.export_molecule():

                f.write(line)

            f.write(' 0.00000   0.00000   0.00000\n')

#    def valid_molecules(self, option:str = 'inner')->list[molecule.molecule]:
#        """This function check the valid pair for the comparison
#        this function should return the pairs of molecule that satisfied the conditions
#        molecules that are not close to boundary should be treated differently from the molecules that are close to boundary
#
#        option: inner or outer"""
#        
#        if (option == 'inner'):
#
#            molecular_list = [m for m in self.molecules if self.is_cut_by_boundary(m) == False]
#
#        elif (option == 'outer'):
#
#            molecular_list = [m for m in self.molecules if self.is_cut_by_boundary(m) == True]
#        
#        else:
#            
#            molecular_list = self.molecules
#
#        return molecular_list
#
#    def valid_pairs(self)-> tuple[list[molecule.molecule], list[molecule.molecule], list[molecule.molecule]]:
#        """This method separate the pairs into inner-inner valid pairs, inner-outer pairs and outer-outer pairs
#
#        option: inner or outer"""
#
#        inner_molecules = self.valid_molecules('inner')
#        outer_molecules = self.valid_molecules('outer')
#
#        
#        inner_inner_pair = []
#        inner_outer_pair = []
#        outer_outer_pair = []
#        for molecule1, molecule2 in itertools.combinations(inner_molecules,2):
#            
#            molecular_distance = molecule.distance(molecule1,molecule2)
#            if (molecular_distance<= self.distance_cut_off):
#
#                inner_inner_pair.append([molecule1, molecule2])
#
#        for molecule1, molecule2 in itertools.combinations(inner_molecules,):
#            
#            molecular_distance = molecule.distance(molecule1,molecule2)
#
#            if (molecular_distance<= self.distance_cut_off):
#
#                inner_inner_pair.append([molecule1, molecule2])
#    
    
@lru_cache(maxsize=1)    
def gen_numerator_matrix(frame:frame, use_full_matrix = True)->np.ndarray:
    """this method is used to get product of charge 1 and charge 2, as all the molecules are the same,
    there is only one numerator matrix throughout the program

    use full matrix only when it is necessary (the file will be big)"""
    element_list=[i.element for i in frame.molecules['1DBT'].atoms]

    if(not use_full_matrix):
        element_matrix = np.tile(element_list,(frame.molecule_size,1))
        charge_matrix = gen_coulomb_matrix.charge_vect(element_matrix, element_matrix.T, False)

    else:
        element_matrix = np.tile(element_list*2,(frame.molecule_size*2,1))
        diagonal_elements = gen_coulomb_matrix.charge_vect(np.array(element_list*2),np.array(element_list*2), True)
        charge_matrix = gen_coulomb_matrix.charge_vect(element_matrix, element_matrix.T, False)
        np.fill_diagonal(charge_matrix,diagonal_elements)

    return charge_matrix

