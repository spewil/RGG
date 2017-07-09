
#Readme:
#Code designed to generate Hard Random Geometric Graphs in a choosen dimension
#Code includes functions for both periodic and unperiodic networks. 

#Inputs: mean degree parameter, Number of nodes and number of dimensions for the embedding space.
#Current functions embedd hard RGGs in the space [0,1]^d 
#Extension: more general domain shapes to be included. 



#Import relevant extensions:
import time
import numpy as np
from scipy.spatial.distance import pdist
from libc.math cimport abs
from scipy import sparse
from scipy.sparse.linalg import eigsh
import scipy
import numpy as np
from scipy import special


#Step function, this is used as the connectivity function for hard random geometric graphs. 
cdef Step(double r,double r_c) :
    
    cdef double P
    
    if r < r_c :
        P = 1.0
    else :
        P = 0.0
        
    return P

#Function to generate Hard RGGs with SOLID BOUNDARY CONDITIONS for 
#choosen values of the number of nodes, means degree and dimension
#extension: alternatively one might want to choose the connection radius for the hard RGG. 

#Input: N, Kappa, d 
#Output: position list of nodes Nxd and edge lists I,J 
def RGG_Generation(int N, float Kappa , int d):

    
	#Variable declarations
    cdef int i 
    cdef int j 
    cdef int Adds = 0 
    cdef double dense
    cdef double r_c
    cdef double[:,:] positions    
    I = [ ]
    J = [ ]
    
    #Define the connection radius for the choosen mean egree number of nodes and dimensions:
    r_c =    (1.0/((3.141592)**0.5) )*(( ((Kappa)/N)*scipy.special.gamma( (d +2.0)/2.0  )   )**(1.0/d ) )


	#Generate N random positions in the domain [0,1]^d
    positions = np.random.uniform(0, 1.0, (N, d))

	#More variable declarations (after the points are generated) 
    cdef int n = positions.shape[0]
    cdef int q = positions.shape[1]
    cdef int p = 2
    cdef double dij
    cdef double dist

    
    #In this loop, for each nodes i and j we:
    #1) Compute the pairwise distance
    #2) Calculate the corresponding connection probability
    #3) Generate an undirected edge between the two nodes with the corresponding probability
    for i in range(N) : 
        for j in range(i+1,N) :        

            dij = 0 

            #Loop over number of dimensions
            for k in range(q):
                # Compute the absolute distance
                dist = abs( positions[i, k] - positions[j, k] )

                # Add to the total distance
                dij = dij + dist**2

            dij = dij**0.5 

            probability = Step(dij , r_c )

            if probability == 1.0 :
                I.append(i)
                J.append(j)
                I.append(j)
                J.append(i)
                Adds = Adds + 1
                    
    return positions ,  I ,J 
    
#Function to generate RGGs with PERIODIC BOUNDARY CONDITIONS
   
