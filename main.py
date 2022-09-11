from frame_data_processing.frame import frame
from frame_data_processing import molecule
import itertools

def main():

    frame0 = frame('Frames/DBT1-00', working_dir='/home/chunhou/Documents/FYP' ,molecule_size=56)
    frame0.print_attributes()
    
    distance_cut_off = 1.2
    for molecule1, molecule2 in itertools.combinations(frame0.molecules,2):
        
        molecular_distance = molecule.distance(molecule1,molecule2)
        if (molecular_distance<= distance_cut_off):

            print(f'pair{molecule1.get_name()},{molecule2.get_name()} -  distance : {molecular_distance}')


if(__name__=='__main__'):

    main()
