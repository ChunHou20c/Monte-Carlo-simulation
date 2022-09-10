"""This module is used to process the frame object into data that can be loaded"""

class frame():
    """The frame class contains the frame number, all the molecules number and also the sides of the box of the frame"""
    def __init__(self, path:str, working_dir:str = '.') -> None:
        """The constructor to read all the data from the text file, take path of the file as argument"""

        with open(path) as frame_data:

        #do some processing first to cleanup unused line
        #first line is the frame number

            first_line = next(frame_data)
            step_number = first_line.split()[-1]
            self.dir_name = f'{working_dir}/{step_number}/' #this last element is the step we want  
            #The step number will be used to  construct the folder name for the molecules and pairs
            
            self.number_of_line = next(frame_data) #this is the number of line to be read as molecular data
            molecular_data = frame_data.readlines()
            self.molecules=frame.molecule(molecular_data[:-1])
            self.box_data = molecular_data[-1]

    class molecule():
        """This is a inner class to store the molecular data"""

        def __init__(self,file_object)->None:
            """Constructor of the molecule object.
            Take file object as argument and separates the molecule into atoms for processing"""
            
            self.raw_data = file_object

        def print_data(self):
            """Function to check the functionality of the molecular class"""

            print(len(self.raw_data))

    def print_attributes(self):
        """This method can be used to check the attributes of the data"""
        
        print(f'directory to save = {self.dir_name}')
        print(f'number of lines to save = {self.number_of_line}')
        print(f'box data = {self.box_data}')

def main():

    frame0 = frame('Frames/DBT1-00')
    frame0.print_attributes()
    frame0.molecules.print_data()

if(__name__=='__main__'):

    main()
