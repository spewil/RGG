

import Adjacency_Map as AMap
import time
import numpy as np
from scipy.spatial.distance import pdist
from libc.math cimport abs
from scipy import sparse
from scipy.sparse.linalg import eigsh
import scipy
import numpy as np
from scipy import special

	
def Expander(int x,Edge_List) :

	"""
	
	
	Parameters
	------------
	
	x : int
	
	Label for the node
	
	Edge_List :  List
	
	List of tuples specifying the edges in the network
	Takes the form [ [ 0 1 ]  [ 0 2 ] [ 1 0 ] [ 2 0 ] ] 
	For example. 
	
	
	Returns 
	--------
	
	Component_List : List
	
	List of the labels of nodes within the same connected component as node i
	
	
	"""
	
	stl_adjacency_map = AMap.AdjacencyMap(Edge_List)
	
	Component_List = [ x ] 
	cdef int change = 1 
	cdef int j 
	cdef int i 
	
	
	while change >= 1   :
		len_1 = len(Component_List)
		
		Current_array_store = Component_List
		
		for j in range( len(Component_List) ) : 
			
			Vec = stl_adjacency_map.get(Component_List[j]  )
			#print("\ncurrents")
			#print(Vec)
			#print(Current_array)
			#print(" ")
			
					
			if Vec != None : 
				
				To_adds = [ ]
				
				for i in range(len(Vec)) :
					
					#print(Vec[i] in Current_array)
					T_F = Vec[i] in Component_List
					if T_F == False :
						
						To_adds.append( Vec[i] )
						
						
				#print("To add " + str(To_adds))       
				Component_List = Component_List + To_adds 
			
		#Terminate if one of the know values if already in it: 
				
		#change = 0 
		
		#print( Current_array )
		
		len_2 = len(Component_List)
		change = len_2 - len_1
		
		
	return Component_List
	
	

def Components(A) : 

	"""
	
	Parameters
	-----------
	
	A : Sparse adjacecny matrix
	
	Returns
	----------
	
	comp :  List
	
	List of Lists containing the indices of the nodes contained
	in the different components in the network
	
	"""
	
	N = A.shape[0]
	
	#Begin by turning sparse matrix into tuples:
	#Each edge will be represented as a tuple
	Edge_List = np.transpose(np.nonzero(A)).astype(np.int32)
	#print("Edge List: "  + str(Edge_List ) ) 
	
	T_F = False
	
	#Empty array containing a list of the nodes which are already specified within a component:
	Done_List = [  ]
    
    #Empty array containing the list of components:
	Comps = [ ] 
	
	cdef int i 
	
	#go through nodes and do the expander
	for i in range(N) : 
		
		T_F = i in Done_List
		
		#print(Done_List)
		if T_F == False :
			#Given a node label and list of edges we 
			#can identify the component which that node belongs to
			Current_Comp = Expander(i,Edge_List) 
			      
			Done_List = Done_List + Current_Comp 
			Comps.append( Current_Comp )
		
	return Comps

