from frame_data_processing.frame import frame

def main():

    frame0 = frame('Frames/DBT1-00', working_dir='/home/chunhou/Documents/FYP' ,molecule_size=56)
    frame0.print_attributes()

    for i in range(10):
        frame0.molecules[0].atoms[i].print_data()

if(__name__=='__main__'):

    main()
