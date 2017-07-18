# gen test 

import matplotlib
import os
import numpy as np
import network_generation.generation as gen
import networkx as nx 
import time 
import json

Kappa, Num_Nodes, d = 20, 100, 2

# make an RGG with solid BCs
start = time.time()
RGGs = gen.RGG_Sample(Kappa,Num_Nodes,d)
print "generation time for RGG_" + RGGs.BC + " of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)

# make an RGG with solid BCs
start = time.time()
RGGp = gen.RGG_Sample(Kappa,Num_Nodes,d,Boundary='p')
print "generation time for RGG_" + RGGp.BC + " of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)

# make an ER graph 
start = time.time()
ER = gen.ER_Sample(Kappa,Num_Nodes)
print "generation time for ER of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)

# compare with networkx
from scipy import special
import scipy
r = (1.0/((3.141592)**0.5) )*(( ((Kappa)/Num_Nodes)*scipy.special.gamma( (d +2.0)/2.0  )   )**(1.0/d ) ) 

start = time.time()
RGG_nx = nx.random_geometric_graph(Num_Nodes,r)
print "generation time for RGG_nx of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)

start = time.time()
ER_nx = nx.gnp_random_graph(Num_Nodes,0.25)
print "generation time for ER_nx of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)

