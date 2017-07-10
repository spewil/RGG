# THIS CODE HAS BEEN ADAPTED FROM WORK DONE BY MATT GARROD, 2017 

#Import required packages:
import numpy as np
from scipy import sparse
import networkx as nx
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt


class Network :
	
		
	def __init__( self , A , Get_Dense_Adj = False ) :

		"""
		Class to store information about a network given the data about it.

		Unless otherwise required we work with sparse matrices since this allows
		computations to be carried out for larger networks with more ease


		Parameters
		-------------

		A: Scipy Sparse Matrix (csr)
		
		Get_Dense_Adj : bool
		
		
		Return
		------------
		
		


		"""
		
		self.Adjacency = A
		#Does size need to be imported as a seperate parameter?
		self.Network_Size = A.shape[0]

		if Get_Dense_Adj == True : 
			self.Dense_Adjacency = A.todense()



		#After initializing we have all the functions to compute properties of the specific network
		#Some of these may be taken from networkx for example. 

		#Functions below give us various properties of the network:

	def Mean_Degree(self) :
		Degree_array = np.array( self.Adjacency.sum(axis=0) ) 
		return np.mean( Degree_array[0] )


	def Min_Degree(self) :
		Degree_array = np.array( self.Adjacency.sum(axis=0) ) 
		return min( Degree_array[0]  ) 
		
	def Max_Degree(self) : 
		Degree_array = np.array( self.Adjacency.sum(axis=0) ) 
		return max( Degree_array[0]  ) 
		
	def Edges(self) :
		"""
		Returns the number of edges in the network
		
		"""
		return len( self.Adjacency.data)/2
		
	def Degrees(self) :
	
		"""
		Degrees:
		
		Returns
		----------
		
		Degrees[0] :  list
		
		A list containing the degrees of the network.
		
		
		
		"""
		degrees = np.array( self.Adjacency.sum(axis=0) )
		return degrees[0]

	def Laplacian(self) :
	
		"""
		
		Constructs the laplacian matrix of the network
		(note is the laplacian defined with the convention L = D - A,
		where D is the diagonal matrix of degrees and A is the adjacency matrix of the network)
		
		Returns
		----------
		
		L : scipy sparse matrix
		
		Laplacian matrix for the network
		
		
		
		"""
	
	
		degrees = np.array( self.Adjacency.sum(axis=0) )
		Is = np.arange(0,self.Network_Size,1)
		D = sparse.coo_matrix((degrees[0], (Is, Is)) , shape = (self.Network_Size,self.Network_Size) , dtype=float).tocsr()
		return  - self.Adjacency + D


	def Mean_Path_Length(self) :
		
		"""
		
		Uses the built in function from networkx to compute the average shortest path
		length of the network.
		
		This method only works if the graph is connected. 
		
		Returns
		---------
		
		Path_Len : float 
		
		Average shortest path length between nodes in the network
		
		#Average shortst path length of the largest connected component of the network
		
		
		"""
		
		#NB: must first compute the fiant comp?? 
		
		A_Giant = self.Largest_Connected_Component()
		G=nx.from_numpy_matrix(A_Giant.todense())
		Path_Len = nx.average_shortest_path_length(G)
		
		return Path_Len
		
	def Av_Clustering_Coefficient(self) :
	
		"""
		
		Calculate the average clustering coefficient of the network
		using the method from networkx
		
		Returns :
		
		Clust : float
		
		Average Clustering Coefficient of the network
		
		
		"""
		
		G=nx.from_numpy_matrix(self.Adjacency.todense())
		Clust = nx.average_clustering(G)
		return Clust
		
	def Degree_Assortativity(self) : 
		#print("Type = " + str(  self.Adjacency ) ) 
		#print( self.Adjacency ) 
		G=nx.from_numpy_matrix(self.Adjacency.todense())
		Assort = nx.degree_assortativity_coefficient(G) 
		return Assort
		
	def Edge_Betweenness_Centrality_Distribution(self) :
		
		G=nx.from_numpy_matrix(self.Adjacency.todense())
		E_B_C_D = nx.edge_betweenness_centrality(G)
		
		return E_B_C_D.values()
		
	def Diameter(self) : 
		G=nx.from_numpy_matrix(self.Adjacency.todense())
		Diam = nx.diameter(G)
		return Diam
		
	def Components(self) :
	
		"""
		Return a list of the connected components in the network
		
		
		"""
		Components = Giant.Components(self.Adjacency  )
		return Components
		
	def Largest_Connected_Component(self) :
	
		"""
		
		Returns: 
		
		Adjacency matrix corresponding to the largest connected component of the network
		
		Optional? positions of the nodes in the largest connected component
		
		"""
	
		index, Giant_Comp_Vec = max(enumerate(self.Components() ), key = lambda tup: len(tup[1]))
		A_Giant = self.Adjacency[Giant_Comp_Vec, :].tocsc()[:, Giant_Comp_Vec]
		A_Giant = A_Giant.tocsr()
		N_Giant = len(Giant_Comp_Vec)
		return A_Giant
		
	def Largest_Connected_Min_Degree(self) :
		A_Giant = self.Largest_Connected_Component()
		Degree_array = np.array( A_Giant.sum(axis=0) ) 
		return min( Degree_array[0]  ) 
		
	def Largest_Comp_Fraction(self) :
		return float(self.Largest_Connected_Component().shape[0])/float(self.Adjacency.shape[0])
		
	def Plot_Network(self, File_Name = "Network" , Positions = None , Show_Plot = False , Colours = None , Size = 20) :
		
		"""
		Plots a network using 
		
		Returns 
		---------
		figure handle??
		
		"""
		
		#(optional argument? directory to save in?? ) 
		
		#Create a figure and clear the current 
		#figure just in case
		
		#If no node positions specified just use standard layout:
		if Positions.all() == None: 
			plt.figure()
			plt.clf
			Gfig = open( File_Name +  '.png' , 'w' )
			graph = nx.from_numpy_matrix(self.Adjacency.todense())
			if Colours == None :
				nx.draw_networkx(graph, node_size= Size, with_labels=False)
			else :
				nx.draw_networkx(graph , node_size = 100 , with_labels=False , node_color = Colours ) 
			
			plt.savefig( Gfig , format = 'png' )
			if Show_Plot == True :
				plt.show()
			# Gfig.close()
			
		else : 
			plt.figure(2)
			plt.clf()
			Gfig = open( File_Name +  '.png' , 'w' )
			graph = nx.from_numpy_matrix(self.Adjacency.todense())
			#nx.draw_spectral(graph, node_size=15, with_labels=False)
			#nx.draw(graph,Positions, node_size=15, with_labels=False)
			if Colours == None :
				nx.draw_networkx(graph, Positions , node_size= Size, with_labels=False)
			else :
				nx.draw_networkx(graph , Positions ,  node_size = 100 , with_labels=False , node_color = Colours ) 
			plt.savefig( Gfig , format = 'png' )
			if Show_Plot == True :
				plt.show()
			# Gfig.close()
			


	def Properties(self, *Additional_Properties) : 
	
		"""
		
		Parameters
		------------
		
		
		Additional_Properties : 
		
		Tuple of strings containing names of the properties of interest
		
		Possible strings are:
		- "Algebraic Connectivity"
		- "Largest Laplace Eigenvalue"
		- "Largest Adjacency Eigenvalue"
		-  "Mean Path Length"
		- "Average Clustering Coefficient"
		
		
		Returns
		-----------
		
		Properties : dict
		
		A dictionary containing the different properties of the ensemble
		Each property is labelled according to its name
		
		
		"""
		
		Property_Dict = { "Mean Degree" : self.Mean_Degree() , "Minimum Degree" : self.Min_Degree() , "Maximum Degree" : self.Max_Degree() }
		Property_Dict["Number of Nodes"] = self.Network_Size
		
		#For any additional functions added we add the property of interest to the array with a name based on the function:
		
		
		"""
		Functional approach - for now just have to hard code things...
		
		
		for i in range( len(Additional_Properties_Functions) ) :
			#if Required_Methods
			Property_Dict[Additional_Properties_Functions[i].__name__] = Additional_Properties_Functions[i]( Required_Methods[i] )
			
		"""
		#Spectral Properties:
		if "Algebraic Connectivity" in Additional_Properties : 
			Property_Dict["Algebraic Connectivity"] = SEV.Smallest_Non_Zero_EV( self.Laplacian() ) 
			
		if "Largest Laplace Eigenvalue" in Additional_Properties :
			Property_Dict["Largest Laplace Eigenvalue"] = SEV.Largest_EV( self.Laplacian() )
		
		if "Largest Adjacency Eigenvalue" in Additional_Properties :
			Property_Dict["Largest Adjacency Eigenvalue"] = SEV.Largest_EV( self.Adjacency ) 
		
		#Other Properties:
		if "Mean Path Length" in Additional_Properties :
			Property_Dict["Mean Path Length"] = self.Mean_Path_Length()
			
		if "Average Clustering Coefficient" in Additional_Properties : 
			Property_Dict["Average Clustering Coefficient"] = self.Av_Clustering_Coefficient()
			
		if "Largest Connected Component Fraction" in Additional_Properties : 
			Property_Dict["Largest Connected Component Fraction"] = self.Largest_Comp_Fraction()
		
		if "Degree Distribution" in Additional_Properties :
			Property_Dict["Degree Distribution"] = self.Degrees()
			
		if "Adjacency Spectrum" in Additional_Properties :
			Property_Dict["Adjacency Spectrum"] = SEV.Sparse_Spectrum( self.Adjacency ) 
		
		if "Laplacian Spectrum" in Additional_Properties :
			Property_Dict["Laplacian Spectrum"] = SEV.Sparse_Spectrum( self.Laplacian() ) 
			print("Sample Laplace" ) 		
		
		if "Degree Assortativity" in Additional_Properties :
			Property_Dict["Degree Assortativity"] = self.Degree_Assortativity() 
			
		if "Largest Connected Minimum Degree" in Additional_Properties :
			Property_Dict["Largest Connected Minimum Degree"] = self.Largest_Connected_Min_Degree()
			
		if "Max Edge Betweenness Centraility" in Additional_Properties :
			Property_Dict["Max Edge Betweenness Centraility"] = max( self.Edge_Betweenness_Centrality_Distribution() ) 
		
		return Property_Dict
		
		


   
