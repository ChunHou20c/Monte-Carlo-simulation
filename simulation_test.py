from electron_coupling import marcus_equation
from frame_data_processing.atom import distance
from simulation import simulation
#import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


test_sim = simulation.Simulation('Frames/DBT1-00')
#test_sim.single_jump('2DBT')
#test_sim.run()

#sns.displot(test_sim.electron_coupling_list)
#plt.savefig('testing1.png')

distance_list = []
electron_coupling_list = []
transfer_rate = []

for vertex in test_sim.graph.vert_dict.values():

    for key, relation in vertex.adjacent.items():

        distance_list.append(relation.distance)

        coupling = test_sim.predicted_electron_coupling(vertex.id, key)

        electron_coupling_list.append(coupling)
        
        rate = marcus_equation.transfer_rate(coupling, 0.18, 300)

        transfer_rate.append(rate)

df = pd.DataFrame(list(zip(distance_list, electron_coupling_list, transfer_rate)), columns = ['distance', 'electron_coupling', 'transfer_rate'])

sns.scatterplot(data=df, x='distance', y='electron_coupling')
plt.savefig('coupling vs distance')

plt.clf()

sns.displot(data=df, x = 'electron_coupling')
plt.savefig('coupling_distribution.png')

df.to_csv('data.csv')



#test_sim.electron_coupling_list
#cm_to_save = []
#
#for vertex in test_sim.graph.vert_dict.values():
#
#    for relation in vertex.adjacent.values():
#
#        cm_to_save.append(relation.coulomb_matrix)
#
#cm_list = np.array(cm_to_save)
#np.save('cm_from_graph.npy', cm_list)
#
