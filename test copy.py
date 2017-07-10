#Import the relevant packages:

import matplotlib
import os
import numpy as np
import network_generation.generation as NetGen
import network_generation.netobject as NetClass

#Pick some parameters:
Mean_Degree_Param , Number_of_Nodes , Dimension = 12.0 , 100 , 2

#Create an instance of the class:
Sampled_Class = NetGen.Sparse_RGG_Sampler(Mean_Degree_Param , Number_of_Nodes , Dimension )

#Get the adjacency matrix:
A = Sampled_Class.Adjacency
print( type(A) ) 

#We can also get the normal python matrix out by typing A.todense()

Network = NetClass.Network(A) 

#Compute the mean degree:
Mean_Deg = Network.Mean_Degree()
print("Mean Degree = " +str(Mean_Deg) )

#We can also call a dictionary of network properties:
print("Properties:\n"  +str(Network.Properties()) )

Positions = np.asarray(Sampled_Class.Positions) 
Network.Plot_Network(File_Name = "Bounded" , Positions = Positions , Show_Plot = True)


# ~~~~~~ Periodic BCs ~~~~~~~

#Create an instance of the class:
Sampled_Class = NetGen.Sparse_RGG_Sampler(Mean_Degree_Param , Number_of_Nodes , Dimension , Boundaries = 'p')

#Get the adjacency matrix:
A = Sampled_Class.Adjacency
print( type(A) ) 

#We can also get the normal python matrix out by typing A.todense()

Network = NetClass.Network(A) 

#Compute the mean degree:
Mean_Deg = Network.Mean_Degree()
print("Mean Degree = " +str(Mean_Deg) )

#We can also call a dictionary of network properties:
print("Properties:\n"  +str(Network.Properties()) )

Positions = np.asarray(Sampled_Class.Positions) 
Network.Plot_Network(File_Name = "Periodic" , Positions = Positions , Show_Plot = True)


# ~~~~~~ Time Test ~~~~~~~

from scipy import special
import scipy

N = 1000
d = 2
Kappa = 10.0
r = (1.0/((3.141592)**0.5) )*(( ((Kappa)/N)*scipy.special.gamma( (d +2.0)/2.0  )   )**(1.0/d ) ) 

import networkx as nx
import time

start = time.time()

G=nx.random_geometric_graph(N,r)
pos=nx.get_node_attributes(G,'pos')
print( G.number_of_edges() ) 
print("Time Taken = " + str(time.time() - start ))

start = time.time()
Sampled_Class = NetGen.Sparse_RGG_Sampler(Kappa , N , d )
A = Sampled_Class.Adjacency
Network = NetClass.Network(A) 
print( Network.Edges() ) 
print("Time Taken = " + str(time.time() - start ))


# about half the time 