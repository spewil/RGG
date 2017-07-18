import matplotlib
import os
import numpy as np
import network_generation.generation as netGen
import network_generation.netobject as netObj
from hopcroftkarp import HopcroftKarp
import time 

#Pick some parameters:
Kappa , Num_Nodes , d = 20 , 10 , 2

# Create an instance of the class:
start = time.time()
RGG_Instance = netGen.RGG_Sample(Kappa,Num_Nodes,d,Boundary='s')
print "Network Generation time: " + str(time.time() - start)

# Make a Network Object 
start = time.time()
RGG = netObj.Network(RGG_Instance)
print "Network Object time: " + str(time.time() - start)

# Get the adjacency matrix:
A = RGG.Adjacency
Adense = RGG.Adense
# Check directedness 
print "Is A symmetric? " + str((Adense == Adense.T).all()) 
# print Adense 

# get the bipartite dictionary representation 
D = RGG.BipartiteDict
# print D

# find the unmatched nodes 
start = time.time()
unmatched = RGG.Find_Unmatched()
print "Matching time: " + str(time.time() - start)
print 'umatched nodes: ' + str(unmatched)

# PLOTTING 
RGG.Plot_Network(unmatched = unmatched) 
print np.asarray(RGG.Positions)

#Compute the mean degree:
Mean_Deg = RGG.Mean_Degree()
print "Mean Degree = " + str(Mean_Deg) 

