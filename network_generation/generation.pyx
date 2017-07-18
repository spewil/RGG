# coding: utf-8
from random import random
from scipy import sparse
from scipy import special
import networkx as nx
from hopcroftkarp import HopcroftKarp 
import numpy as np
import pickle
import matplotlib.pyplot as plt
import scipy


class BaseEnsemble(object):

	def __init__(self, *params):
		raise Exception('need to override this')

	def create_sample(self):
		raise Exception('need to override this')

	def generate_samples(self, n=100):
		self.samples = [self.create_sample() for i in xrange(n)]

	def to_dict(self):
		raise Exception('need to override this')

	@classmethod
	def from_dict(self, data):
		raise Exception('need to override this')

	def to_disk(self, folder='rgg_samples'):
		raise Exception('need to override this')

	@classmethod
	def from_disk(self, path):
		raise Exception('need to override this')


class RGGEnsemble(BaseEnsemble):
	def __init__(self, kappa, n, d, shortcut_prob=0, boundary='s', samples=None):
		# Params 
		self.kappa = kappa
		self.n = n
		self.d = d
		self.shortcut_prob = shortcut_prob
		self.boundary = boundary
		# a list of sample objects 
		self.samples = samples or []
		self.num_samples = len(self.samples)
		# labels the type of ensemble 
		self.label = 'RGG'

		# param_string 
		if abs(float(shortcut_prob)) == 0.0 : 
			self.param_string = "RGG_K_" + str(kappa) + "_N_" + str(n) + "_d_" + str(d) + "_boundary_" + str(boundary)
		else :
			self.param_string = "RGG_K_" + str(kappa) + "_N_" + str(n) + "_d_" + str(d) + "_boundary_" + str(boundary) + '_short_' + str(shortcut_prob)

	def create_sample(self):
		kappa = self.kappa
		n = self.n
		d = self.d
		shortcut_prob = self.shortcut_prob
		boundary = self.boundary

		#Define the variables used:
		cdef double dom_size = 1
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
		kappa = float(kappa)
		# inverted analytic function for r = f(K)
		r_c = (1.0/((3.141592)**0.5) )*(( ((kappa)/n)*scipy.special.gamma( (d +2.0)/2.0  )   )**(1.0/d ) )

		#Define arrays to store the integer labels of the connected pairs:
		source = [ ]
		target = [ ]

		#Randomly generate the node positions in the hypercube:
		# N-list of d-lists
		positions = np.random.uniform(0, 1.0, (n, d))

		# number of nodes -- nn for cython 
		cdef int nn = positions.shape[0]
		# number of dimensions
		cdef int q = positions.shape[1]

		for i in range(nn) :
			for j in range(i+1,nn) :

				dij = 0

				#Loop over number of dimensions
				for k in range(q):
					# Compute the absolute distance
					dist = abs( positions[i, k] - positions[j, k] )

					#Extra condition for periodic boundarys:
					if boundary == 'p' :
						if dist>0.5*dom_size :
							dist = dom_size - dist
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
						source.append(i)
						target.append(j)
					else: # j --> i
						source.append(j)
						target.append(i)

		# store the generated sample as a sample object with these S,T,P and params 
		self.samples.append(RGGSample(source, target, positions, self.kappa, self.n, self.d, self.boundary, self.shortcut_prob))

		# param_string 
		if abs(float(shortcut_prob)) == 0.0 : 
			self.param_string = "RGG_K_" + str(kappa) + "_N_" + str(n) + "_d_" + str(d) + "_boundary_" + str(boundary)
		else :
			self.param_string = "RGG_K_" + str(kappa) + "_N_" + str(n) + "_d_" + str(d) + "_boundary_" + str(boundary) + '_short_' + str(shortcut_prob)

	## SAVING AND LOADING DATA 

	def to_disk(self, folder='rgg_samples'):
		filename = './%s/%s.pickle' % (folder, self.param_string())
		with open(filename, 'w') as fout:
			pickle.dump(self.to_dict(), fout)

	@classmethod
	def from_disk(self, path):
		with open(path, 'r') as fin:
			data = pickle.load(fin)
		return EREnsemble.from_dict(data)

class RGGSample(object):
	def __init__(self, source, target, positions, kappa, n, d, boundary, shortcut_prob):
		self.source = source
		self.target = target
		self.positions = positions
		self.kappa = kappa
		self.n = n
		self.d = d
		self.boundary = boundary
		self.shortcut_prob = shortcut_prob

		self.adjacency = None
		self.adjacency_dense = None
		self.bipartite_dict = None

		self.calculate_adjacency_representation()
		self.calculate_dictionary_representation()

	def calculate_adjacency_representation(self):
		### Create A ###
		#combine the edge lists in order to form a sparse adjacency matrix:
		# scipy.sparse.coo_matrix(arg1, shape=None, dtype=None, copy=False)
		# arg1 is (data, row, column)
		A = sparse.coo_matrix(
				(np.ones_like(self.source), (self.source, self.target)),
				shape=(self.n, self.n),
				dtype=float)
		self.adjacency = A.tocsr()
		self.adjacency_dense = self.adjacency.todense()

	def calculate_dictionary_representation(self):
		
		### Create D ###
		# bipartite dictionary representation of directed graph

		# build dictionary from the edge list (python)
		nums = xrange(self.n)
		# ['1s', '2s', '3s', ...] for the lefthand side of the bipartite
		nodelist = (str(num) + 's' for num in nums)

		# generate an empty dict of lists for targets
		D = {k: [] for k in nodelist}

		# go through and add targets to source key -- target sets
		for idx, val in enumerate(self.source):
			D[str(val)+'s'].append(self.target[idx])

		self.bipartite_dict = D

	def find_unmatched(self):

		self.calculate_dictionary_representation 
		graph = self.bipartite_dict

		# list of node indices 0,1,2,...,N 
		unmatched = range(self.n)

		# find matched nodes
		# (they will be doubled (i-->j and j-->i) ) 
		matched = HopcroftKarp(graph).maximum_matching().values()
		# remove the matched nodes to leave the unmatched 
		for node in range(self.n):
			for match in matched:
				# if node exists in matching
				if node == match:
					# remove it 
					unmatched.remove(node)
		
		return unmatched  

	def plot_network(self, unmatched=None, size=20, node_label=False):
		
		if self.d != 2:
			raise ValueError('The graph is not 2D!')

		else: 
			File_Name = str(self.param_string) + '_plot'
			
			plt.figure()
			plt.clf
			Gfig = open( File_Name +  '.eps' , 'w' )
			
			# make a networkx graph object 
			graph = nx.DiGraph(self.Adjacency.todense())

			# make position and label dicts 
			posDict = {}
			for i in range(len(self.positions)):
				posDict[i] = self.positions[i] 

			# color the nodes 
			if unmatched != None: 
				colorList = []
				for i in range(self.n):
					if i in unmatched:
						colorList.append('green')
					else:
						colorList.append('red')

			# draw the network 
			if node_label: 
				nx.draw_networkx(graph, pos=posDict, with_labels=True, node_color=colorList, node_size=size, width = 0.3)
			else: 
				nx.draw_networkx(graph, pos=posDict, with_labels=False, node_color=colorList, node_size=size, width = 0.3)

			plt.savefig(Gfig, format='eps', dpi=1000)
			# plt.savefig( Gfig , format = 'png' )
		
	def mean_degree(self) :

		degree_array = np.array( self.adjacency_dense.sum(axis=0) ) 

		# 2x because of directedness  
		return 2*np.mean( degree_array )

	def properties(self, *additional_properties):
		raise NotImplementedError()

	def to_dict(self):
		return {
			'kappa': self.kappa,
			'n': self.n,
			'd': self.d,
			'boundary': self.boundary,
			'shortcut_prob': self.shortcut_prob,
			'source': self.source,
			'target': self.target,
			'positions': self.positions,
		}

	@classmethod
	def from_dict(cls, data):
		return RGGSample(
			kappa=data['kappa'],
			n=data['n'],
			d=data['d'],
			boundary=data['boundary'],
			shortcut_prob=data['shortcut_prob'],
			source=data['source'],
			target=data['target'],
			positions=data['positions'],
		)

class EREnsemble(BaseEnsemble):
	
	def __init__(self, kappa, n, samples=None):
		self.kappa = kappa
		self.n = n
		self.samples = samples or []

	def create_sample(self) :

		"""
		Generates sparse adjacency matrix for a directed Erdos Renyi Random Graph


		Parameters
		------------

		kappa : double

		Mean degree parameter:

		N : int

		Number of nodes in the network
		"""
		kappa = self.kappa
		n = self.n

		# Compute the probability required to get the desired mean degree:
		p = kappa/float(n)

		# Empty array to store edge lists:

		source = []
		target = []

		for i in range(n) :
			for j in range(i+1,n) :

				u = np.random.uniform(0,1.0)

				if u < p :
					# 1/2 probability for direction
					if random() >= 0.5: # i --> j
						source.append(i)
						target.append(j)
					else: # j --> i
						source.append(j)
						target.append(i)

		# Attributes
		self.source = source
		self.target = target
		# Params
		self.label = 'ER'
		self.kappa = kappa
		self.n = n
		# param_string
		self.param_string = "ER_K_" + str( n*p ) + "_N_" + str(n)

class ERSample : 

	def __init__(self,source,target,kappa,n) : 
		self.source = source
		self.target = target
		self.kappa = kappa
		self.n = n

		self.adjacency = None
		self.adjacency_dense = None
		self.bipartite_dict = None

		self.calculate_adjacency_representation()
		self.calculate_dictionary_representation()

	def param_string(self):
		# param_string
		if abs(float(self.shortcut_prob)) == 0.0 :
			return "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_boundary_" + str(self.boundary)

		return "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_boundary_" + str(self.boundary) + '_short_' + str(self.shortcut_prob)

	def to_dict(self):
		return {
			'kappa': self.kappa,
			'n': self.n,
			'd': self.d,
			'boundary': self.boundary,
			'shortcut_prob': self.shortcut_prob,
			'samples': [s.to_dict() for s in self.samples]
		}

	@classmethod
	def from_dict(self, data):
		return RGGEnsemble(
			kappa=data['kappa'],
			n=data['n'],
			d=data['d'],
			boundary=data['boundary'],
			shortcut_prob=data['shortcut_prob'],
			samples=[ERSample.from_dict(s) for s in data['samples']],
		)

	def to_disk(self, folder='er_samples'):
		filename = './%s/%s.pickle' % (folder, self.param_string())
		with open(filename, 'w') as fout:
			pickle.dump(self.to_dict(), fout)

	@classmethod
	def from_disk(self, path):
		with open(path, 'r') as fin:
			data = pickle.load(fin)
		return EREnsemble.from_dict(data)

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
