from functools import lru_cache
from frame_data_processing.atom import atom
from frame_data_processing import atom as Atom

class molecule:
    """This is a class to store the molecular data as this class is used for DBT1 the index for S and N atoms are fixed
    to use this class for other molecule, the atom indices need to be changed"""
    S1_index = 43
    S2_index = 6
    N_index = 21

    def __init__(self,molecular_data:list[str])->None:
        """Constructor of the molecule object.
            Take file object as argument and separates the molecule into atoms for processing"""

        self.raw_data = molecular_data
        self.name = molecular_data[0].split()[0]
        self.atoms = [atom(line) for line in self.raw_data]

    def __str__(self) -> str:
        """print the name of the molecule"""
        
        _str = f"{self.name}"
        return _str


    def print_data(self):
        """Function to check the functionality of the molecular class"""
        
        return_string = ''

        for line in self.raw_data:

            print(line.replace('\n', ''))
            return_string+=line

        return return_string

    def process_data(self, border: float, process_x: bool, process_y: bool, process_z: bool)->None:
        """This method process the molecules that is cut by the boundary and reconstruct the whole raw data of the molecule,
        this function is expected to be used together with the cut by boundary method in the frame object
        
        arguments: 
        border - the border value of the box
        process_x - whether x value is cut by boundary
        process_y - whether y value is cut by boundary
        process_z - whether z value is cut by boundary
        """
        
        tmp_raw_data = ''
        for line in self.raw_data:
            

            datas = line.split()
                
            if (process_x):
                #process x by adding the boundary value to the lower value of the x coordinates
                
                if (datas[-3]<5):
                    
                    modified_value = float(datas[-3])+border
                    datas[-3] = f'{modified_value:.3f}'

            if (process_y):
                #process y by adding the boundary value to the lower value of the x coordinates

                if (datas[-2]<5):
                    
                    modified_value = float(datas[-2])+border
                    datas[-2] = f'{modified_value:.3f}'

            if (process_z):
                #process z by adding the boundary value to the lower value of the x coordinates

                if (datas[-1]<5):
                    
                    modified_value = float(datas[-1])+border
                    datas[-1] = f'{modified_value:.3f}'
        
            new_line = ' '.join(datas)
            tmp_raw_data += new_line
    
        self.raw_data = tmp_raw_data

    def get_name(self):
        """getter function to return the name of the molecule with index for example 1DBT, 2DBT"""

        return self.name

    def get_xyz_list(self)-> tuple[list[float], list[float], list[float]]:

        x_list = [i.x for i in self.atoms]
        y_list = [i.y for i in self.atoms]
        z_list = [i.z for i in self.atoms]

        return x_list, y_list, z_list
    
@lru_cache(maxsize=1)
def distance(m1:molecule,m2:molecule)->float:
    """This method calculates the distance between the two molecule"""
    
    distance_N = Atom.distance(m1.atoms[molecule.N_index], m2.atoms[molecule.S1_index])+\
            Atom.distance(m1.atoms[molecule.N_index], m2.atoms[molecule.S2_index])+\
            Atom.distance(m1.atoms[molecule.N_index], m2.atoms[molecule.N_index])

    distance_S1 = Atom.distance(m1.atoms[molecule.S1_index], m2.atoms[molecule.S1_index])+\
            Atom.distance(m1.atoms[molecule.S1_index], m2.atoms[molecule.S2_index])+\
            Atom.distance(m1.atoms[molecule.S1_index], m2.atoms[molecule.N_index])

    distance_S2 = Atom.distance(m1.atoms[molecule.S2_index], m2.atoms[molecule.S1_index])+\
            Atom.distance(m1.atoms[molecule.S2_index], m2.atoms[molecule.S2_index])+\
            Atom.distance(m1.atoms[molecule.S2_index], m2.atoms[molecule.N_index])
    
    Cut_off_distance = ((distance_N/3) + (distance_S1/3) + (distance_S2/3))/3

    return Cut_off_distance
