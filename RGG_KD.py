import numpy as np
from scipy import spatial
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import rcParams

# make the plots nicer 
rcParams.update({'figure.autolayout': True,'mathtext.fontset': 'stix','font.family' : 'STIXGeneral', 'axes.axisbelow' : 'line'})

nnodes = 1000

for r in arange(.01,.1,10)
	positions =  np.random.rand(nnodes,2)
	kdtree = spatial.KDTree(positions)
	pairs = kdtree.query_pairs(r)
	G = nx.Graph()
	G.add_nodes_from(range(nnodes))
	G.add_edges_from(list(pairs))
	
	# pos = dict(zip(range(nnodes),positions))
	# nx.draw(G,pos,node_size=1)

fig = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('Average Degree $\langle{k}\rangle$',fontsize=18)
ax.set_ylabel('Radius $r$',fontsize=18)

plt.show()