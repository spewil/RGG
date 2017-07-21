# gen test

import matplotlib
import os
import numpy as np
import network_generation.generation as ng
import experiment as ex
import networkx as nx
import time
import json

# RGGE = ex.RGGExperiment( kappa_range, n, d, shortcut_prob=0, boundary='s')

# for each of these N-d-BC ensembles with kappa_range: 
# 	build the experiment object from the data 
# 	store the gamma parameter in an N-d matrix 
#		list of 3-tuples: [(d,N,gamma),(d,N,gamma),...(d,N,gamma)]

# Experiment = ng.RGGExperiment(kappa_range,Num_Nodes,d,boundary='s')
# print Experiment.ensembles

# make an RGG Ensemble with solid boundarys
# RGG = ng.RGGEnsemble(kappa,Num_Nodes,d,boundary='s')
# generate one sample 
# start = time.time()
# RGG.generate_samples(n=1)
# RGG.samples[0].plot_network()
# print "generation time for RGG_" + RGGEs.boundary + " of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)
# export to JSON 
# RGGEs.to_disk()

# del RGGEs

# ~~~

# # import data from JSON
# RGGEimport = ng.RGGEnsemble(kappa,Num_Nodes,d,boundary='s')
# fn = RGGEimport.get_param_string()
# RGGEimport = ng.RGGEnsemble.from_disk('rgg_samples/'+fn+'.pickle')

# # find unmatched nodes 
# um0 = RGGEs.samples[0].find_unmatched()
# print um0
# # plot network 
# RGGEs.samples[0].plot_network(unmatched=um0)
# # find the mean degree
# print RGGEs.samples[0].mean_degree()

# # ~~~

# # make an ER graph
# start = time.time()
# ER = ng.EREnsemble(kappa,Num_Nodes)
# ER.create_sample()
# print "generation time for ER of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)

# ER.samples[0].plot_network()
# # find unmatched nodes 
# um0 = ER.samples[0].find_unmatched()
# print um0
# # plot network 
# ER.samples[0].plot_network(unmatched=um0)
# # find the mean degree
# print ER.samples[0].mean_degree()

# # compare with networkx
# from scipy import special
# import scipy
# r = (1.0/((3.141592)**0.5) )*(( ((kappa)/Num_Nodes)*scipy.special.gamma( (d +2.0)/2.0  )   )**(1.0/d ) )

# start = time.time()
# RGG_nx = nx.random_geometric_graph(Num_Nodes,r)
# print "generation time for RGG_nx of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)

# start = time.time()
# ER_nx = nx.gnp_random_graph(Num_Nodes,0.25)
# print "generation time for ER_nx of " + str(Num_Nodes) + " nodes: " + str(time.time() - start)

