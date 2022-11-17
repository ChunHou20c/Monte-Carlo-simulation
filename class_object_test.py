"""This script is use for class object testing, atom, molecule"""

from frame_data_processing import atom
from frame_data_processing import molecule
from frame_data_processing import frame

def atom_test():

    with open("test_files/atom1.txt", 'r') as f:
        
        str1 = f.read()
        atom1 = atom.atom(str1)

        atom1.check_string()
        atom1.print_data()

        atom1.modify_coordinate(0.912,'x')
        atom1.modify_coordinate(0.912,'y')
        atom1.modify_coordinate(0.912,'z')

        atom1.print_data()
        atom1.check_string()

def molecule_test():
    
    with open("test_files/3DBT_cut_in_half.txt", 'r') as f:
        
        m_string = f.readlines()

        molecule1 = molecule.molecule(m_string)

        print(molecule1)
        
        molecule1.check_condition()

        molecule1.atom_processing()

        content_to_write = molecule1.raw_data

    with open("test_files/3DBT_processed", 'w') as f:

        for line in content_to_write:

            f.write(line)

def frame_test():

    working_dir = '/home/chunhou/Documents/FYP'
    molecule_size=56

        
    frame0 = frame.frame('Frames/DBT1-00', working_dir=working_dir,
            molecule_size=molecule_size, cut_off_distance= 1.2,
            Use_full_CM = False)
    
    #print(frame0.molecules.get_vertices())
    frame0.export_molecule('3DBT','test_files/test_molecule.gro')


if __name__ == '__main__':
    
    #atom_test()
    #molecule_test()
    frame_test()
