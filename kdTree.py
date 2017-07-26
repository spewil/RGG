from scipy import spatial
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np 
import time 

nnodes = 10**4
print nnodes
r = 0.25
start = time.time()
positions =  np.random.rand(nnodes,2)
kdtree = spatial.KDTree(positions)
pairs = kdtree.query_pairs(r)
print "generation time" + str(time.time() - start)

# print pairs
# G = nx.Graph()
# G.add_nodes_from(range(nnodes))
# G.add_edges_from(list(pairs))
# pos = dict(zip(range(nnodes),positions))
# nx.draw_networkx(G,pos)
# plt.show()