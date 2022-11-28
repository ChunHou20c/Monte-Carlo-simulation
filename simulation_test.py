from simulation import simulation
import numpy as np

test_sim = simulation.Simulation('Frames/DBT1-00')
#test_sim.run()

cm_to_save = []

for vertex in test_sim.graph.vert_dict.values():

    for relation in vertex.adjacent.values():

        cm_to_save.append(relation.coulomb_matrix)

cm_list = np.array(cm_to_save)
np.save('cm_from_graph.npy', cm_list)

