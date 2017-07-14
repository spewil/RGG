#Import the relevant packages:

import matplotlib.pyplot as plt
import os
import cProfile
import numpy as np
from network_generation.HopcroftKarp import HopcroftKarp
import network_generation.generation as NetGen
import network_generation.netobject as NetClass
import time 
from matplotlib import rcParams

# make the plots nicer 
rcParams.update({'figure.autolayout': True,'mathtext.fontset': 'stix','font.family' : 'STIXGeneral', 'axes.axisbelow' : 'line'})

# ~~~ FUNCTIONS 

# find the unmatched nodes 
# must be a better way to do this, but this works for now
def findUnmatched(graph):
	# grab the list of nodes (the lefthand side of the bipartite)
	# turn it into ints 0,1,2,...,N
	unmatched = list(range(len(graph.keys())))
	# grab the matched nodes, they will be doubled (i-->j and j-->i)
	# this is a python list 
	matched = HopcroftKarp(graph.copy()).maximum_matching().values()
	# chuck the matched nodes 
	for node in range(len(graph.keys())):
		for match in matched:
			if node == match:
				unmatched.remove(node)
	return unmatched  

# ~~~~ GENERATE NET

#Pick some parameters:
Number_of_Nodes , Dimension = 1000 , 3

# Experimental Quantities 
ND_avg = []
Avg_deg = []
Avg_var = []
ND_var = []

# 15 steps between 1 and 30 
for Kp in range(1,30,2):
	# ten trials of each param
	# reset the ND list and average degrees
	ND = []
	Mean_Degs = []
	for trial in range(10):

		# Create an instance of the class:
		Sampled_Class = NetGen.Sparse_RGG_Sampler(Kp , Number_of_Nodes , Dimension , Directed = True )

		# Create instance of Network Class:
		Network = NetClass.Network(Sampled_Class.Adjacency) 

		# Compute the mean degree and save:
		Mean_Degs.append(Network.Mean_Degree())

		# Find the number of unmatched nodes 
		ND.append(len(findUnmatched(Sampled_Class.Dictionary)))

	# average the ten trials, get stddev 
	Avg_deg.append(np.mean(Mean_Degs))
	Avg_var.append(np.std(Mean_Degs))
	ND_avg.append(np.mean(ND))
	ND_var.append(np.std(ND))

# scale the ND values by the number of nodes 
ND_avg = np.array(ND_avg)/Number_of_Nodes
ND_var = np.array(ND_var)/Number_of_Nodes

# plot result 
plt.figure()
plt.errorbar(Avg_deg, ND_avg, yerr = ND_var)
plt.xlabel('Mean Degree $\\langle{k}\\rangle$')
plt.ylabel('Fraction of Driver Nodes $n_D$')
plt.title('d = 3, N = 1000, 10 Trials ')
plt.grid(True)
plt.savefig("d3N1000.png")
# plt.show()







