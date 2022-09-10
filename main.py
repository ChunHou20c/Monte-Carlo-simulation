from frame_data_processing.frame import frame

def main():

    frame0 = frame('Frames/DBT1-00')
    frame0.print_attributes()
    frame0.molecules.print_data()
    frame0.molecules.atoms[0].print_data()

if(__name__=='__main__'):

    main()
