"""This script is use for class object testing, atom, molecule"""

from frame_data_processing import atom
from frame_data_processing import molecule

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
    
    with open("test_files/molecule1.txt", 'r') as f:
        
        m_string = f.readlines()

        molecule1 = molecule.molecule(m_string)

        print(molecule1)

if __name__ == '__main__':

    molecule_test()
