# gen test

import matplotlib
import os
import numpy as np
import network_generation.generation as ng
import networkx as nx
import time
import json

kappa, Num_Nodes, d = 30, 1000, 2

# make an RGG Ensemble with solid boundarys
RGGEs = ng.RGGEnsemble(kappa,Num_Nodes,d)
# generate one sample 
start = time.time()
RGGEs.create_sample()
print "generation time for RGG_" + RGGEs.boundary + " of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)
# export to JSON 
RGGEs.to_disk()

del RGGEs

# ~~~

# import data from JSON
RGGEimport = ng.RGGEnsemble(kappa,Num_Nodes,d)
fn = RGGEimport.get_param_string()
RGGEimport = ng.RGGEnsemble.from_disk('rgg_samples/'+fn+'.pickle')

# find unmatched nodes 
um0 = RGGEimport.samples[0].find_unmatched()
print um0
# plot network 
RGGEimport.samples[0].plot_network(unmatched=um0)
# find the mean degree
print RGGEimport.samples[0].mean_degree()


# # make an RGG with periodic boundarys
# start = time.time()
# RGGp = ng.RGGSample(kappa,Num_Nodes,d,boundary='p')
# print "generation time for RGG_" + RGGp.boundary + " of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)

# # make an ER graph
# start = time.time()
# ER = ng.ERSample(kappa,Num_Nodes)
# print "generation time for ER of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)

# compare with networkx
from scipy import special
import scipy
r = (1.0/((3.141592)**0.5) )*(( ((kappa)/Num_Nodes)*scipy.special.gamma( (d +2.0)/2.0  )   )**(1.0/d ) )

start = time.time()
RGG_nx = nx.random_geometric_graph(Num_Nodes,r)
print "generation time for RGG_nx of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)

# start = time.time()
# ER_nx = nx.gnp_random_graph(Num_Nodes,0.25)
# print "generation time for ER_nx of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)

