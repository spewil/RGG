# THIS CODE HAS BEEN ADAPTED FROM WORK DONE BY MATT GARROD, 2017 

import numpy as np
from scipy import sparse
import scipy
from scipy import special
from random import random

def Top_Hat(r, r_c, pr) :
	
	"""
	Top hat function with a = 1 , b = pr

	Parameters
	-----------
	r : float 

	Value we wish to evaluate the top-hat function at

	r_c : float

	Width of the top hat

	pr : float

	height of the lower section of the top hat function. This corresponds 
	to the probability of having long-range shortcuts in the RGG class. 


	Returns
	----------
	P : float 

	Value of the top hat function at r

	"""
	
	if r < r_c :
		P = 1.0
	else :
		P = pr
		
	return P


class RGG_Sample : 

	def __init__(self, Kappa, N, d, shortcut_prob = 0, Boundary = 's') :
	

		"""
		
		Generates a directed random geometric graph in a regions [ 0,1]^d with solid or periodic boundary
		conditions.

		Parameters
		------------
		
		Kappa: float

		The mean degree parameter for the network. This controls the expected number of connections possed by
		a node in the bulk of the domain. For Solid Boundary conditions the observed mean degree of the network may be smaller

		N: int

		Number of nodes in the network

		d: int

		The number of dimensions of the hypercube to be generated

		shortcut_prob : float in [0,1]
		
		Controls the probability of adding random shortcuts between nodes. This is the height of the lower part of a 'Top-Hat' 
		Kernel like connection function

		Maximum value 1.0

		Returns
		-------

		S, T -- source and target node lists 
		Positions -- list of node coordinate lists 

		Future Additions
		----------------

		Add the Gaussian boundary condition case 
		Optimize the code... 

		"""		

		#Define the variables used:
		cdef double Dom_Size = 1
		cdef int i 
		cdef int j 
		cdef int Adds = 0 
		cdef double dense
		cdef double r_c
		cdef double[:,:] positions
		cdef int p = 2
		cdef double dij
		cdef double dist

		#Make sure the mean degree is a float so that the next computation works:
		Kappa = float(Kappa)
		# inverted analytic function for r = f(K)
		r_c = (1.0/((3.141592)**0.5) )*(( ((Kappa)/N)*scipy.special.gamma( (d +2.0)/2.0  )   )**(1.0/d ) ) 
		
		#Define arrays to store the integer labels of the connected pairs:
		S = [ ]
		T = [ ]
		
		#Randomly generate the node positions in the hypercube:
		# N-list of d-lists
		positions = np.random.uniform(0, 1.0, (N, d))

		# number of nodes 
		cdef int n = positions.shape[0]
		# number of dimensions 
		cdef int q = positions.shape[1]
				
		for i in range(N) : 
			for j in range(i+1,N) : 

				dij = 0 

				#Loop over number of dimensions
				for k in range(q):
					# Compute the absolute distance
					dist = abs( positions[i, k] - positions[j, k] )
					
					#Extra condition for periodic BCs:
					if Boundary == 'p' :
						if dist>0.5*Dom_Size :
							dist = Dom_Size - dist
					# Add to the total distance
					dij = dij + dist**2

				dij = dij**0.5             
				
				#Compute the connection probability:
				# returns 1 if within radius, shortcut_prob otherwise 
				probability = Top_Hat(dij , r_c , shortcut_prob)
				
				u = np.random.uniform(0,1.0)
			
				# for zero shortcut probability, always true in the radius 
				# for nonzero baseline, might be connecting long-range shortcuts 
				if u < probability :
					# 1/2 probability for direction
					if random() >= 0.5: # i --> j
						S.append(i) 
						T.append(j) 
					else: # j --> i 
						S.append(j)
						T.append(i)

		# source, target, and positions lists 
		self.Source = np.asarray(S)
		# self.Target = T
		# self.Positions = np.asarray(positions)

		# # Params
		# self.Label = 'RGG'
		# self.K = Kappa
		# self.N = N
		# self.d = d
		# self.BC = Boundary

		# # Param_String 
		# if abs(float(shortcut_prob)) == 0.0 : 
		# 	self.Param_String = "RGG_K_" + str(Kappa) + "_N_" + str(N) + "_d_" + str(d) + "_BC_" + str(Boundary)
		# else :
		# 	self.Param_String = "RGG_K_" + str(Kappa) + "_N_" + str(N) + "_d_" + str(d) + "_BC_" + str(Boundary) + '_short_' + str(shortcut_prob)
		

class ER_Sample : 

	def __init__(self,Kappa,N) : 

		"""
		Generates sparse adjacency matrix for an Erdos Renyi Random Graph
		
		 
		Parameters
		------------
		
		Kappa : double
		
		Mean degree parameter:

		N : int 
		
		Number of nodes in the network
		
		
		"""
		# Compute the probability required to get the desired mean degree:
		p = Kappa/float(N)
		
		# Empty array to store edge lists:
		
		S = [ ] 
		T = [ ]
		
		for i in range(N) : 
			for j in range(i+1,N) :
				
				u = np.random.uniform(0,1.0)
				
				if u < p :
					# 1/2 probability for direction
					if random() >= 0.5: # i --> j
						S.append(i) 
						T.append(j) 
					else: # j --> i 
						S.append(j)
						T.append(i)

		# Attributes 
		self.Source = S
		self.Target = T
		# Params 
		self.Label = 'ER'
		self.K = Kappa
		self.N = N
		# Param_String 
		self.Param_String = "ER_K_" + str( N*p ) + "_N_" + str(N)  
		
		

		
		
		
		
		
		
		
		
	   
