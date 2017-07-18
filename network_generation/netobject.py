
#Import required packages:
import numpy as np
from scipy import sparse
import networkx as nx
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from hopcroftkarp import HopcroftKarp
import pydot

class Network :
	
	###############	
	# CONSTRUCTOR #
	###############

	def __init__( self , Sample) :

		"""
		Class to store information about a network given the data about it.

		Unless otherwise required we work with sparse matrices since this allows
		computations to be carried out for larger networks with more ease


		Parameters
		-------------

		Sample
			A network generation object 
				Source
				Target
				Positions
		"""
		
		# make this work OO version... this is stupid 
		# Params
		self.Label = Sample.Label
		self.K = Sample.K
		self.N = Sample.N
		self.d = Sample.d
		self.BC = Sample.BC
		self.Param_String = Sample.Param_String

		# Data 
		self.S = Sample.Source 
		self.T = Sample.Target 
		# if it's an RGG, store the positions 
		if self.Label == 'RGG':
			self.Positions = Sample.Positions

		### Create A ### 

		#combine the edge lists in order to form a sparse adjacency matrix:       
		# scipy.sparse.coo_matrix(arg1, shape=None, dtype=None, copy=False)
		# arg1 is (data, row, column)                        
		A = sparse.coo_matrix((np.ones_like(self.S), (self.S, self.T)) , shape = (self.N,self.N) , dtype=float ).tocsr()
		self.Adjacency = A
		self.Adense = A.todense()

		### Create D ###

		# bipartite dictionary representation of directed graph 

		# build dictionary from the edge list (python)
		nodelist = range(self.N)
		# ['1s', '2s', '3s', ...] for the lefthand side of the bipartite 
		nodelist = [str(num) + 's' for num in nodelist]
		
		# generate an empty dict of lists for targets 
		D = {k: [] for k in nodelist}
		
		# go through and add targets to source key -- target sets
		for idx, val in enumerate(self.S):
			D[str(val)+'s'].append(self.T[idx])

		self.BipartiteDict = D

	#

	def Find_Unmatched(self):

		graph = self.BipartiteDict.copy() 

		# list of node indices 0,1,2,...,N 
		unmatched = range(self.N)

		# find matched nodes
		# (they will be doubled (i-->j and j-->i) ) 
		matched = HopcroftKarp(graph).maximum_matching().values()
		# remove the matched nodes to leave the unmatched 
		for node in range(self.N):
			for match in matched:
				# if node exists in matching
				if node == match:
					# remove it 
					unmatched.remove(node)
		return unmatched  

	############
	# PLOTTING #
	############

	def Plot_Network(self, unmatched = None, Size = 20, label=False) :            

		if self.d != 2:
			raise ValueError('The graph is not 2D!')

		else: 
			File_Name = str(self.Param_String) + '_plot'
			
			plt.figure()
			plt.clf
			Gfig = open( File_Name +  '.png' , 'w' )
			
			# make a networkx graph object 
			graph = nx.DiGraph(self.Adjacency.todense())

			# make position and label dicts 
			posDict = {}
			for i in range(len(self.Positions)):
				posDict[i] = self.Positions[i] 

			# color the nodes 
			if unmatched != None: 
				colorList = []
				for i in range(self.N):
					if i in unmatched:
						colorList.append('green')
					else:
						colorList.append('red')

			# draw the network 
			if label: 
				nx.draw_networkx(graph, pos=posDict, with_labels=True, node_color=colorList, node_size=Size, width = 0.3)
			else: 
				nx.draw_networkx(graph, pos=posDict, with_labels=False, node_color=colorList, node_size=Size, width = 0.3)

			plt.savefig( Gfig , format = 'png' )

			Gdot = nx.nx_pydot.to_pydot(graph)   #converts G to Gp in pydot format
 
			# ... some code that makes the Gp pydot graph prettier
			# the default for a DiGraph will have arrows
			# i will often do the following
			# edges = Gp.get_edge_list()
			# for e in edges:
			#      source = e.get_source()
			#      target = e.get_destination()
			#      if source == 'some string that is the name of the source with edges that you want to make look different' and target == 'some string2':
			#        e.set_style('bold')
			#        e.set_color('red')
			#        e.set_arrowhead('tee')  # use tee instead of arrow
			#        e.set_label('label for edge')
			# can do similar stuff with nodes
			 
			# now output your graph to a file and display it
			File_Name = str(self.Param_String) + '_plot'
			Gdot.write_png(File_Name + '_dot.png', prog='dot')  
			# writes Gp to png file #use prog='neato' or 'fdp' for undirected graph (no default arrows for this)
			# the next 2 lines open and display the png file
			# im = Image.open(File_Name + '_dot.png') 
			# im.show()
			 

	####################
	# NETWORK MEASURES #
	####################

	def Mean_Degree(self) :
		Degree_array = np.array( self.Adjacency.sum(axis=0) ) 
		# 2x because of directedness  
		return 2*np.mean( Degree_array[0] )

	def Properties(self, *Additional_Properties) : 
	
		"""
		
		Parameters
		------------
		
		
		Additional_Properties : 
		
		Tuple of strings containing names of the properties of interest
		
		Possible strings are:
		- "Algebraic Connectivity"
		
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
			
		return Property_Dict
		
		


   
