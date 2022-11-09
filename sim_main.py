"""This is the main simulation script"""
from frame_data_processing import simulation
from frame_data_processing import frame


working_dir = '/home/chunhou/Documents/FYP'
molecule_size=56

#this part is to prevent recalculation of numerator matrix
frame0 = frame.frame(f'Frames/DBT1-00', working_dir=working_dir,
        molecule_size=molecule_size, cut_off_distance= 1.2,
        Use_full_CM = False)

sim1 = simulation.simulation(frame0)

sim1.check_frame_data()
