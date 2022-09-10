"""This module is used to process the frame object into data that can be loaded"""
from typing import Generator
from frame_data_processing.molecule import molecule

class frame():
    """The frame class contains the frame number, all the molecules number and also the sides of the box of the frame"""
    def __init__(self, path:str, working_dir:str = '.', molecule_size:int=56) -> None:
        """The constructor to read all the data from the text file, take path of the file as argument"""

        self.molecule_size = molecule_size
        with open(path) as frame_data:

        #do some processing first to cleanup unused line
        #first line is the frame number

            first_line = next(frame_data)
            step_number = first_line.split()[-1]
            self.dir_name = f'{working_dir}/{step_number}/' #this last element is the step we want  
            #The step number will be used to  construct the folder name for the molecules and pairs
            
            self.number_of_line = next(frame_data) #this is the number of line to be read as molecular data

            molecular_data = frame_data.readlines()
            self.box_data = molecular_data[-1]

            molecular_data_to_process = molecular_data[:-1]
            list_of_molecule = self.molecule_split(molecular_data_to_process)
            self.molecules=[molecule(m) for m in list_of_molecule]

    def molecule_split(self, raw_data)->list[list[str]]:
        """This function splits the aggregated molecular raw data into single molecule
        as for now the number of atoms in a molecule is known and fixed"""

        def split(list_to_split:list[str],chunk_size:int)->Generator:
            """generator to split the data"""

            for i in range(0,len(list_to_split), chunk_size):
                yield list_to_split[i:i+chunk_size]

        return list(split(raw_data,self.molecule_size))

    def print_attributes(self):
        """This method can be used to check the attributes of the data"""
        
        print(f'directory to save = {self.dir_name}')
        print(f'number of lines to save = {self.number_of_line}')
        print(f'box data = {self.box_data}')
        print(f'no of atoms in a single molecule = {self.molecule_size}')

