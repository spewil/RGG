#Import the relevant packages:

import matplotlib
import os
import numpy as np
from network_generation.HopcroftKarp import HopcroftKarp
from hopcroftkarp import HopcroftKarp
import network_generation.generation as NetGen
import network_generation.netobject as NetClass
# import networkx as nx
# from networkx.algorithms import bipartite
import time 

# ~~~ FUNCTIONS 

# find the unmatched nodes 
# must be a better way to do this, but this works for now, using the functions at hand
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
Number_of_Nodes , Dimension = 1000 , 2

# 15 steps between 1 and 30 
for Kp in range(1,30,2):
	# ten trials of each param
	for trial in range(10):




#Create an instance of the class:
Sampled_Class = NetGen.Sparse_RGG_Sampler(Mean_Degree_Param , Number_of_Nodes , Dimension , Directed = True )

#Create instance of Network Class:
Network = NetClass.Network(Sampled_Class.Adjacency) 

#Compute the mean degree:
print 'mean degree param: ' + str(Mean_Degree_Param)
Mean_Deg = Network.Mean_Degree()
print("Mean Degree = " +str(Mean_Deg) )

# #We can also call a dictionary of network properties:
# print("Properties:\n"  +str(Network.Properties()) )

# Positions = np.asarray(Sampled_Class.Positions) 
# Network.Plot_Network(File_Name = "Bounded" , Positions = Positions) # , Show_Plot = True)

# # get the dictionary
D = Sampled_Class.Dictionary
# print D
# print len(D.keys())
# print 'dict: ' + str(D)
# for v in D.values():
# 	if len(v) == 1:
# 		print 'one edge'
# 	if v == {}:
# 		print 'empty set'



print 'unmatched nodes: ' + str(findUnmatched(D))

# # ~~~~~~ Periodic BCs ~~~~~~~

# #Create an instance of the class:
# Sampled_Class = NetGen.Sparse_RGG_Sampler(Mean_Degree_Param , Number_of_Nodes , Dimension , Boundaries = 'p')

# #Get the adjacency matrix:
# A = Sampled_Class.Adjacency
# print( type(A) ) 

# #We can also get the normal python matrix out by typing A.todense()

# Network = NetClass.Network(A) 

# #Compute the mean degree:
# Mean_Deg = Network.Mean_Degree()
# print("Mean Degree = " +str(Mean_Deg) )

# #We can also call a dictionary of network properties:
# print("Properties:\n"  +str(Network.Properties()) )

# Positions = np.asarray(Sampled_Class.Positions) 
# Network.Plot_Network(File_Name = "Periodic" , Positions = Positions , Show_Plot = True)


# ~~~~~~ Time Test ~~~~~~~

# from scipy import special
# import scipy

# N = 1000
# d = 2
# Kappa = 10.0
# r = (1.0/((3.141592)**0.5) )*(( ((Kappa)/N)*scipy.special.gamma( (d +2.0)/2.0  )   )**(1.0/d ) ) 

# import networkx as nx
# import time

# start = time.time()

# G=nx.random_geometric_graph(N,r)
# pos=nx.get_node_attributes(G,'pos')
# print( G.number_of_edges() ) 
# print("Time Taken = " + str(time.time() - start ))

# start = time.time()
# Sampled_Class = NetGen.Sparse_RGG_Sampler(Kappa , N , d )
# A = Sampled_Class.Adjacency
# Network = NetClass.Network(A) 
# print( Network.Edges() ) 
# print("Time Taken = " + str(time.time() - start ))

# # about half the time 