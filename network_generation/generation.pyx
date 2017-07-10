# THIS CODE HAS BEEN ADAPTED FROM WORK DONE BY MATT GARROD, 2017 

import numpy as np
from scipy import sparse
import scipy
from scipy import special


def Top_Hat(r, r_c,pr) :
	
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


class Sparse_Erdos_Sample : 

	def __init__(self,N,Kappa) : 

		"""
		Generates sparse adjacency matrix for an Erdos Renyi Random Graph
		
		 
		Parameters
		------------
		
		N : int 
		
		Number of nodes in the network
		
		Kappa : double
		
		Mean degree parameter:
		
		
		"""
		#Compute the probability required to get the desired mean degree:
		p = Kappa/float(N)
		
		#Empty array to stroe edge lists:
		
		I = [ ] 
		J = [ ]
		
		for i in range(N) : 
			for j in range(i+1,N) :
				
				u = np.random.uniform(0,1.0)
				
				if u < p :
					I.append(i)
					J.append(j)
					I.append(j)
					J.append(i)
		#Combine these to obtain a sparse matrix:
		A = sparse.coo_matrix((np.ones_like(I), (I, J)) , shape = (N,N) , dtype=float ).tocsr()
		
		self.Adjacency = A
		self.Label = 'Erdos'
		self.Param_String = "_N_" + str(N) + "_K_" + str( N*p ) 


class Sparse_RGG_Sampler : 

	def __init__(self,Kappa,N,d, shortcut_prob = 0,  Dom_Size = 1.0 , Return_Pos = False , NODE_POS = None , Label = 'RGG', **Optional_Parameters) :
	

		"""
		
		Generates a random geometric graph in a regions [ 0,Dom_Size]^d with solid or periodic boundary
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

		Dom_Size : float 

		The side length of the hypercube. The default value gives the unit hypercube.

		Return_Pos : Bool

		Specify whether it is possible to access positions in this instance of the class. For larger graphs this
		is member intensize. 

		NODE_POS
		The option exists to add in predefined node positions for the network. 
		"""
		#Get the string specifying the boundary conditions from the optional arguments:
		Boundaries = Optional_Parameters.get('Boundaries', None)
		
		#If the shortcut prob is included then we should classify the ensemble name differently:
		

		#Define the variables used:
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
		r_c =    (1.0/((3.141592)**0.5) )*(( ((Kappa)/N)*scipy.special.gamma( (d +2.0)/2.0  )   )**(1.0/d ) ) 
		
		#Define arrays to store the integar labels of the connected pairs:
		I = [ ]
		J = [ ]
		
		#Randomly sample the node positions:
		if NODE_POS == None :
			positions = np.random.uniform(0, 1.0, (N, d))
		else :
			positions = NODE_POS
		
		cdef int n = positions.shape[0]
		cdef int q = positions.shape[1]
		
		n = positions.shape[0]
		q = positions.shape[1]
		
		
		for i in range(N) : 
			for j in range(i+1,N) : 

				dij = 0 

				#Loop over number of dimensions

				for k in range(q):
					# Compute the absolute distance
					dist = abs( positions[i, k] - positions[j, k] )
					
					#Extra condition for periodic BCs:
					if Boundaries == 'P' or Boundaries == 'p' :
						if dist>0.5*Dom_Size :
							dist = Dom_Size - dist
					# Add to the total distance
					dij = dij + dist**2

				dij = dij**0.5             
				
				#Compute the connection probability:
				probability = Top_Hat(dij , r_c , shortcut_prob )


				u = np.random.uniform(0,1.0)
			
				if u < probability :
					I.append(i)
					J.append(j)
					I.append(j)
					J.append(i)
		 


		#combine the edge lists in order to form a sparse adjacency matrix:                               
		A = sparse.coo_matrix((np.ones_like(I), (I, J)) , shape = (N,N) , dtype=float ).tocsr()
		
		self.Adjacency =  A
		self.Positions = positions
		self.Label = Label
		if abs(float(shortcut_prob)) == 0.0 : 
			self.Param_String = "_N_" + str(N) + "_K_" + str(Kappa) + "_d_" + str(d) + "_BC_" + str(Boundaries)
		else :
			self.Param_String = "_N_" + str(N) + "_K_" + str(Kappa) + "_d_" + str(d) + "_BC_" + str(Boundaries) + '_short_' + str(shortcut_prob)
		
		"""
		
		#Return the adjacency matrix and the positions of the nodes if required.
		if Return_Pos == True :
			return A , positions
		else :
			return A	
			
		"""


#General Soft RGG
#Call the connection function of interest as an argument
#Perhaps also call some item which describes the boundaries in some way. 


class Soft_RGG_Sampler : 

	def __init__(self,Kappa,N,d, Connection_Function , shortcut_prob = 0,  Dom_Size = 1.0 , Return_Pos = False , NODE_POS = None , Label = 'RGG', **Optional_Parameters) :
	

		"""
		
		def __init__(self,Kappa,N,d, shortcut_prob = 0, Boundaries = 'P' , Dom_Size = 1.0 , Return_Pos = False , NODE_POS = None , Label = 'RGG', **Optional_Parameters) :
		
		Generates a soft random geometric graph in a regions [ 0,Dom_Size]^d with solid or periodic boundary
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

		Dom_Size : float 

		The side length of the hypercube. The default value gives the unit hypercube.

		Return_Pos : Bool

		Specify whether it is possible to access positions in this instance of the class. For larger graphs this
		is member intensize. 

		NODE_POS
		The option exists to add in predefined node positions for the network. 
		"""
		#Get the string specifying the boundary conditions from the optional arguments:
		Boundaries = Optional_Parameters.get('Boundaries', None)
		
		#If the shortcut prob is included then we should classify the ensemble name differently:
		

		#Define the variables used:
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
		r_c =    (1.0/((3.141592)**0.5) )*(( ((Kappa)/N)*scipy.special.gamma( (d +2.0)/2.0  )   )**(1.0/d ) ) 
		
		#Define arrays to store the integar labels of the connected pairs:
		I = [ ]
		J = [ ]
		#Randomly sample the node positions:
		if NODE_POS == None :
			positions = np.random.uniform(0, 1.0, (N, d))
		else :
			positions = NODE_POS
		
		cdef int n = positions.shape[0]
		cdef int q = positions.shape[1]
		
		n = positions.shape[0]
		q = positions.shape[1]
		
		
		for i in range(N) : 
			for j in range(i+1,N) : 

				dij = 0 

				#Loop over number of dimensions

				for k in range(q):
					# Compute the absolute distance
					dist = abs( positions[i, k] - positions[j, k] )
					
					#Extra condition for periodic BCs:
					if Boundaries == 'P' or Boundaries == 'p' :
						if dist>0.5*Dom_Size :
							dist = Dom_Size - dist
					# Add to the total distance
					dij = dij + dist**2

				dij = dij**0.5             
				
				probability = 0
				"""
				Replace the above line with a the chosen conneciton function!
				
				"""


				u = np.random.uniform(0,1.0)
			
				if u < probability :
					I.append(i)
					J.append(j)
					I.append(j)
					J.append(i)
		 


		#combine the edge lists in order to form a sparse adjacency matrix:                               
		A = sparse.coo_matrix((np.ones_like(I), (I, J)) , shape = (N,N) , dtype=float ).tocsr()
		
		self.Adjacency =  A
		self.Positions = positions
		self.Label = Label
		if abs(float(shortcut_prob)) == 0.0 : 
			self.Param_String = "_N_" + str(N) + "_K_" + str(Kappa) + "_d_" + str(d) + "_BC_" + str(Boundaries)
		else :
			self.Param_String = "_N_" + str(N) + "_K_" + str(Kappa) + "_d_" + str(d) + "_BC_" + str(Boundaries) + '_short_' + str(shortcut_prob)
		
		"""
		
		#Return the adjacency matrix and the positions of the nodes if required.
		if Return_Pos == True :
			return A , positions
		else :
			return A	
			
		"""	
		
def Standard_String_Generator( Net_Type, **Parameters ) :

	"""
	Generates the string to be used at the end of data file names for the 
	specifed network type
	
	Parameters
	------------
	
	Net_Type : str
	
	String specifying the network type
	
	e.g 'RGG' or 'Erdos' 
	
	Parameters : 
	
	list of parameters of the ensemble 
	
	
	"""


	N = Parameters.get('N', None)
	Kappa = Parameters.get('Kappa', None)
	d = Parameters.get('d', None)
	Boundaries = Parameters.get('Boundaries', None)
	Shortcut_Prob = Parameters.get('Shortcut_Prob' , None  ) 
	p = Parameters.get('p', None)
	if Net_Type == 'RGG' :
		
		#Account for two different cases. That with shortcuts and that without:
		if Shortcut_Prob == None :
			Param_String =  "_N_" + str(N) + "_K_" + str(Kappa) + "_d_" + str(d) + "_BC_" + str(Boundaries)
		elif abs(float(Shortcut_Prob)) == 0.0 :
			Param_String =  "_N_" + str(N) + "_K_" + str(Kappa) + "_d_" + str(d) + "_BC_" + str(Boundaries)
				
		else :
			Param_String =  "_N_" + str(N) + "_K_" + str(Kappa) + "_d_" + str(d) + "_BC_" + str(Boundaries) + '_short_' + str(Shortcut_Prob)
			 
	
	elif Net_Type == 'Erdos' : 
		Param_String = "_N_" + str(N) + "_K_" + str( N*p )
	
	else : 
		print("Unrecognised string type") 
		
	
	return Param_String 
		
		
		
		
		
		
		
		
		
		
	   
