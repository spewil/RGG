# coding: utf-8
from random import random
from scipy import sparse
from scipy import special
from scipy import stats
import scipy.spatial as sp
import numpy as np
import pickle
import matplotlib.pyplot as plt
import networkx as nx
from hopcroftkarp import HopcroftKarp
import Giant_Extract as Giant

class BaseEnsemble(object):
	def __init__(self, *params):
		raise Exception('need to override this')

	def create_sample(self):
		raise Exception('need to override this')

	def generate_samples(self, n=100):
		self.samples = [self.create_sample() for i in xrange(n)]

	def generate_samples_kd(self, n=100):
		self.samples = [self.create_sample_kd() for i in xrange(n)]

	def to_dict(self):
		raise Exception('need to override this')

	@classmethod
	def from_dict(self, data):
		raise Exception('need to override this')

	def to_disk(self, folder=''):
		raise Exception('need to override this')

	@classmethod
	def from_disk(self, path):
		raise Exception('need to override this')

class RGGEnsemble(BaseEnsemble):
	def __init__(self, kappa, n, d, shortcut_prob=0, boundary='s', samples=[], num_radii=1):
		self.kappa = kappa
		self.n = n
		self.d = d
		self.shortcut_prob = shortcut_prob
		self.boundary = boundary
		self.samples = samples
		self.num_samples = self.get_num_samples()

		# radius computation 
		self.r_list = [self.compute_radius() for r in xrange(num_radii)]
		self.radius = np.mean(self.r_list)
		self.radius_std = stats.sem(self.r_list)

	def get_num_samples(self):
		return len(self.samples)

	def find_nD_list(self):
		sample_nD_list = []
		for sample in self.samples:
			sample_nD_list.append(sample.find_num_unmatched()/float(sample.n))	
		return sample_nD_list	

	def find_nD_list_LCC(self):
		sample_nD_list = []
		for sample in self.samples:
			sample_nD_list.append(sample.find_num_unmatched()/float(sample.LCC_size))	
		return sample_nD_list		

	def compute_radius(self):
		# generate a single sample (positions)
		# compute the distance list 
		# sort and find the r_critical, store to ensemble 

		Kappa = self.kappa
		N = self.n
		d = self.d
		shortcut_prob = self.shortcut_prob
		boundary = self.boundary

		#Define the variables used:
		cdef double dom_size = 1
		cdef int i
		cdef int j
		cdef int k
		cdef double[:,:] positions
		cdef double dij
		cdef double dist

		#Define arrays to store the integer labels of the connected pairs:
		source= [ ]
		target= [ ]

		# N-list of d-lists
		# Randomly generate positions in a d-dimensional standard Normal distribution
		if boundary == 'g':
			positions = np.random.multivariate_normal(np.zeros(d), np.eye(d), size=N)
		#Randomly generate the node positions in the hypercube:
		else:
			positions = np.random.uniform(0, 1.0, (N, d))

		# number of nodes
		cdef int n = positions.shape[0]
		# number of dimensions
		cdef int q = positions.shape[1]

		# compute distance list 
		d_list = [ ]

		# boundary == 's' or 'g' 
		if boundary != 'p':
			d_list = sp.distance.pdist(positions)
		# boundary == 'p'
		else: 
			for i in range(n) :
				for j in range(i+1,n) :
					dij = 0
					#loop over number of dimensions
					for k in range(q):
						# compute the absolute distance
						dist = abs( positions[i, k] - positions[j, k] )
						if dist>0.5*dom_size :
							dist = dom_size - dist
						# add to the total distance
						dij = dij + dist**2
					dij = dij**0.5
					d_list.append(dij)

		d_list_sorted = np.sort(d_list)
		
		# find the value required to get 
		# the required number of edges 
		# above the required threshold:
		# <k> = 2E/N
		num_edges = int((Kappa*N)/2.0)

		# return the critical radius 
		return d_list_sorted[num_edges]		

	def create_sample(self):
		Kappa = self.kappa
		N = self.n
		d = self.d
		shortcut_prob = self.shortcut_prob
		boundary = self.boundary

		#Define the variables used:
		cdef double dom_size = 1
		cdef int i
		cdef int j
		cdef int k
		cdef double r_c
		cdef double[:,:] positions
		cdef double dij
		cdef double dist

		#Define arrays to store the integer labels of the connected pairs:
		source= [ ]
		target= [ ]

		# N-list of d-lists
		# Randomly generate positions in a d-dimensional standard Normal distribution
		if boundary == 'g':
			positions = np.random.multivariate_normal(np.zeros(d), np.eye(d), size=N)
		#Randomly generate the node positions in the hypercube:
		else:
			positions = np.random.uniform(0, 1.0, (N, d))

		# number of nodes
		cdef int n = positions.shape[0]
		# number of dimensions
		cdef int q = positions.shape[1]

		r_c = self.radius

		for i in range(n) :
			for j in range(i+1,n) :
				dij = 0
				#Loop over number of dimensions
				for k in range(q):
					# Compute the absolute distance
					dist = abs( positions[i, k] - positions[j, k] )
					if boundary == 'p':
						if dist>0.5*dom_size :
							dist = dom_size - dist
					# Add to the total distance
					dij = dij + dist**2
				dij = dij**0.5

				# compute the connection probability:
				# returns 1 if within radius, shortcut prob otherwise
				if dij < r_c :
					probability = 1.0
				else :
					probability = shortcut_prob

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

		return RGGSample(
			    source=source,
				target=target,
				positions=np.asarray(positions),
				kappa=Kappa,
				n=N,
				d=d,
				boundary=boundary,
				shortcut_prob=shortcut_prob,
				)

	def create_sample_kd(self):
		Kappa = self.kappa
		N = self.n
		d = self.d
		shortcut_prob = self.shortcut_prob
		boundary = self.boundary

		#Define arrays to store the integer labels of the connected pairs:
		source= [ ]
		target= [ ]

		# N-list of d-lists
		# Randomly generate positions in a d-dimensional standard Normal distribution
		if boundary == 'g':
			positions = np.random.multivariate_normal(np.zeros(d), np.eye(d), size=N)
		#Randomly generate the node positions in the hypercube:
		else:
			positions = np.random.uniform(0, 1.0, (N, d))

		r_c = self.radius

		# make k-d tree 
		if boundary == 'p':
			tree = sp.cKDTree(positions, copy_data = True, boxsize=1.)
		else:
			tree = sp.cKDTree(positions, copy_data = True)

		# find the list of nearest-neighbor pairs 
		n_pairs = tree.query_pairs(r_c,output_type='ndarray')

		# for each 2-vec in the array of neighbor pairs 
		for t in n_pairs:
			# flip coin 
			if np.random.random() < 0.5:
				source.append(t[0])
				target.append(t[1])
			else:
				source.append(t[1])
				target.append(t[0])

		return RGGSample(
			    source=source,
				target=target,
				positions=np.asarray(positions),
				kappa=Kappa,
				n=N,
				d=d,
				boundary=boundary,
				shortcut_prob=shortcut_prob,
				)

	def to_LCC(self):
		for sample in self.samples:
			sample.to_LCC()

	def get_param_string(self):
		# param_string
		if abs(float(self.shortcut_prob)) == 0.0 :
			return "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + str(self.boundary) + "_num_r_" + str(len(self.r_list))
		else:	
			return "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + str(self.boundary) + "_num_r_" + str(len(self.r_list)) + '_short_' + str(self.shortcut_prob)

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
			samples=[RGGSample.from_dict(s) for s in data['samples']],
		)

	def to_disk(self, folder='/Users/spencerw/Documents/DATA/rgg_samples'):
		filename = '%s/%s.pickle' % (folder, self.get_param_string())
		with open(filename, 'w') as fout:
			pickle.dump(self.to_dict(), fout)

	@classmethod
	def from_disk(self, path):
		with open(path, 'r') as fin:
			data = pickle.load(fin)
		return RGGEnsemble.from_dict(data)

class RGGSample(object):
	def __init__(self, source, target, positions, kappa, n, d, boundary, shortcut_prob):
		self.source = source
		self.target = target
		self.positions = np.asarray(positions)
		self.kappa = kappa
		self.n = n
		self.d = d
		self.boundary = boundary
		self.shortcut_prob = shortcut_prob

		# # param_string
		# if abs(float(self.shortcut_prob)) == 0.0 :
		# 	self.param_string = "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + str(self.boundary)
		# else:
		# 	self.param_string = "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + str(self.boundary) + '_short_' + str(self.shortcut_prob)

		self.adjacency = None
		self.adjacency_dense = None
		self.bipartite_dict = None

		self.calculate_adjacency_representation()
		self.calculate_dictionary_representation()

	def get_param_string(self):
		# param_string
		if abs(float(self.shortcut_prob)) == 0.0 :
			return "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + str(self.boundary)
		else:	
			return "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + str(self.boundary) + '_short_' + str(self.shortcut_prob)

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

	def to_LCC(self):

		#### update the adjacency 
		# make it symmetric so the calculation works 
		symmetric = self.adjacency + self.adjacency.transpose(copy=True)

		# Till's code  
		index, LCC_vec = max(enumerate(Giant.Components(symmetric)), key = lambda tup: len(tup[1]))

		# take only nodes that are in the LCC_vec
		LCC_adjacency_dense = self.adjacency_dense[np.ix_(LCC_vec,LCC_vec)]
		LCC_adjacency = sparse.csr_matrix(LCC_adjacency_dense)

		# add sample object attributes 
		self.adjacency = LCC_adjacency
		self.adjacency_dense = LCC_adjacency_dense
		self.LCC_size = len(LCC_vec)
		self.LCC_vec = LCC_vec
		self.positions = self.positions[LCC_vec]

		#### update the dictionary 
		# Note: if node is in the LCC, then all its neighbors are, too!
		# pop all values that AREN'T in the LCC
		for node in range(self.n):
			if node not in LCC_vec:
				self.bipartite_dict.pop(str(node)+'s',None)

	def find_unmatched(self):

		graph = self.bipartite_dict.copy() 
		num_nodes = len(graph.keys())

		# list of node indices 0,1,2,...,N 
		unmatched = range(num_nodes)

		# find matched nodes
		# (they will be doubled (i-->j and j-->i) ) 
		matched = HopcroftKarp(graph).maximum_matching().values()
		# remove the matched nodes to leave the unmatched 
		for node in range(num_nodes):
			for match in matched:
				# if node exists in matching
				if node == match:
					# remove it 
					unmatched.remove(node)
		return unmatched  

	def find_num_unmatched(self):
		unmatched = self.find_unmatched()

		return len(unmatched)

	def plot_network(self, unmatched=None, size=20, label_nodes=False, save=True, LCC=False):
		
		if self.d != 2:
			raise ValueError('The graph is not 2D!')
		else:		
			if LCC == True: 
				File_Name = './plots/' + str(self.get_param_string()) + '_LCC'+ '_plot' + '.eps'
			else:
				File_Name = './plots/' + str(self.get_param_string()) + '_plot' + '.eps'
				
			plt.figure()
			plt.clf
			Gfig = open( File_Name , 'w' )
			
			# make a networkx graph object 
			graph = nx.DiGraph(self.adjacency_dense)

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
				if label_nodes: 
					nx.draw_networkx(graph, pos=posDict, with_labels=True, node_color=colorList, node_size=size, width = 0.3)
				else: 
					nx.draw_networkx(graph, pos=posDict, with_labels=False, node_color=colorList, node_size=size, width = 0.3)
			else:
				if label_nodes: 
					nx.draw_networkx(graph, pos=posDict, with_labels=True, node_size=20, width = 0.3)
				else: 
					nx.draw_networkx(graph, pos=posDict, with_labels=False, node_size=20, width = 0.3)
	
			if self.boundary != 'g':	
				plt.xlim([0,1])
				plt.ylim([0,1])

			if save == True: 
				plt.savefig( File_Name , format='eps', dpi=1000 )

	def mean_degree(self) :
		# sum down columns (in-degree)
		Degree_array = np.array( self.adjacency_dense.sum(axis=0) ) 
		# 2x because of directedness  
		# [0] index because A is a np.matrix type!! 
		return 2*np.mean( Degree_array[0] )

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

	def to_disk(self, folder='/Users/spencerw/Documents/DATA/rgg_samples'):
		filename = '%s/%s.pickle' % (folder, self.get_param_string())
		with open(filename, 'w') as fout:
			pickle.dump(self.to_dict(), fout)

	@classmethod
	def from_disk(self, path):
		with open(path, 'r') as fin:
			data = pickle.load(fin)
		return RGGSample.from_dict(data)

class EREnsemble(BaseEnsemble):
	def __init__(self, kappa, n, samples=None):
		self.kappa = kappa
		self.n = n
		self.samples = samples or []

	def get_num_samples(self):
		return len(self.samples)

	def find_nD_list(self):
		sample_nD_list = []
		for sample in self.samples:
			sample_nD_list.append(sample.find_num_unmatched()/float(sample.n))	
		return sample_nD_list	

	def find_nD_list_LCC(self):
		sample_nD_list = []
		for sample in self.samples:
			sample_nD_list.append(sample.find_num_unmatched()/float(sample.LCC_size))	
		return sample_nD_list	

	def create_sample(self) :
		Kappa = self.kappa
		N = self.n
		# Params
		self.Label = 'ER'

		# Compute the probability required to get the desired mean degree:
		p = Kappa/float(N)

		# Empty array to store edge lists:

		source = [ ]
		target = [ ]

		for i in range(N) :
			for j in range(i+1,N) :

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
		
		return ERSample(
			    source=source,
				target=target,
				kappa=Kappa,
				n=N,
				)

	def to_LCC(self):
		for sample in self.samples:
			sample.to_LCC()

	def get_param_string(self):
		# param_string
		return "ER_K_" + str(self.kappa) + "_N_" + str(self.n) 

	def to_dict(self):
		return {
			'kappa': self.kappa,
			'n': self.n,
			'samples': [s.to_dict() for s in self.samples]
		}

	@classmethod
	def from_dict(self, data):
		return EREnsemble(
			kappa=data['kappa'],
			n=data['n'],
			samples=[ERSample.from_dict(s) for s in data['samples']],
		)

	def to_disk(self, folder='/Users/spencerw/Documents/DATA/er_samples'):
		filename = '%s/%s.pickle' % (folder, self.get_param_string())
		with open(filename, 'w') as fout:
			pickle.dump(self.to_dict(), fout)

	@classmethod
	def from_disk(self, path):
		with open(path, 'r') as fin:
			data = pickle.load(fin)
		return EREnsemble.from_dict(data)

class ERSample(object):
	def __init__(self, source, target, kappa, n):
		self.source = source
		self.target = target
		self.kappa = kappa
		self.n = n

		self.adjacency = None
		self.adjacency_dense = None
		self.bipartite_dict = None

		self.calculate_adjacency_representation()
		self.calculate_dictionary_representation()

	def get_param_string(self):
		# param_string
		return "ER_K_" + str(self.kappa) + "_N_" + str(self.n) 

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

	def to_LCC(self):

		#### update the adjacency 
		# make it symmetric so the calculation works 
		symmetric = self.adjacency + self.adjacency.transpose(copy=True)

		# Till's code  
		index, LCC_vec = max(enumerate(Giant.Components(symmetric)), key = lambda tup: len(tup[1]))

		# take only nodes that are in the LCC_vec
		LCC_adjacency_dense = self.adjacency_dense[np.ix_(LCC_vec,LCC_vec)]
		LCC_adjacency = sparse.csr_matrix(LCC_adjacency_dense)

		# add sample object attributes 
		self.adjacency = LCC_adjacency
		self.adjacency_dense = LCC_adjacency_dense
		self.LCC_size = len(LCC_vec)
		self.LCC_vec = LCC_vec

		#### update the dictionary 
		# convert the LCC adjacency into a dictionary 
		nums = xrange(self.LCC_size)
		nodelist = (str(num) + 's' for num in nums)
		D = {k: [] for k in nodelist}
		
		# row=out, col=in 
		# nonzero returns tuple of row,col indices of nonzeros 
		outgoing, incoming = (np.nonzero(LCC_adjacency_dense))
		for i in xrange(len(outgoing)):
			D[str(outgoing[i])+'s'].append(incoming[i])

		self.bipartite_dict = D

	def find_unmatched(self):

		graph = self.bipartite_dict.copy() 
		num_nodes = len(graph.keys())

		# list of node indices 0,1,2,...,N 
		unmatched = range(num_nodes)

		# find matched nodes
		# (they will be doubled (i-->j and j-->i) ) 
		matched = HopcroftKarp(graph).maximum_matching().values()
		# remove the matched nodes to leave the unmatched 
		for node in range(num_nodes):
			for match in matched:
				# if node exists in matching
				if node == match:
					# remove it 
					unmatched.remove(node)
		return unmatched  

	def find_num_unmatched(self):
		unmatched = self.find_unmatched()

		return len(unmatched)

	def plot_network(self, unmatched=None, size=20, label_nodes=False, show=False):

		File_Name = 'plots/' + str(self.get_param_string()) + '_plot' + '.eps'
		
		plt.figure()
		plt.clf
		Gfig = open( File_Name , 'w' )
		
		# make a networkx graph object 
		graph = nx.DiGraph(self.adjacency_dense)

		# color the nodes 
		if unmatched != None: 
			colorList = []
			for i in range(self.n):
				if i in unmatched:
					colorList.append('green')
				else:
					colorList.append('red')
			# draw the network 
			if label_nodes: 
				nx.draw_random(graph, with_labels=True, node_color=colorList, node_size=size, width = 0.3)
			else: 
				nx.draw_random(graph, with_labels=False, node_color=colorList, node_size=size, width = 0.3)
		else:
			if label_nodes: 		
				nx.draw_random(graph, with_labels=True, node_size=20, width = 0.3)
			else: 
				nx.draw_random(graph, with_labels=False, node_size=20, width = 0.3)
		
		if show:
			plt.show()
		else:
			plt.savefig( File_Name , format='eps', dpi=1000 )

	def mean_degree(self) :
		Degree_array = np.array( self.adjacency_dense.sum(axis=0) ) 
		# 2x because of directedness  
		return 2*np.mean( Degree_array[0] )

	def properties(self, *additional_properties):
		raise NotImplementedError()

	def to_dict(self):
		return {
			'kappa': self.kappa,
			'n': self.n,
			'source': self.source,
			'target': self.target,
		}

	@classmethod
	def from_dict(cls, data):
		return ERSample(
			kappa=data['kappa'],
			n=data['n'],
			source=data['source'],
			target=data['target'],
		)


######### Gaussian ########## 

class GaussianEnsemble(BaseEnsemble):
	def __init__(self, kappa, n, d, shortcut_prob=0, samples=[], num_radii=1):
		self.kappa = kappa
		self.n = n
		self.d = d
		self.shortcut_prob = shortcut_prob
		self.samples = samples
		self.num_samples = self.get_num_samples()

		# we compute the radius num_radii times for every sample 
		# radius computation 
		self.r_list = [self.compute_radius() for r in xrange(num_radii)]
		self.radius = np.mean(self.r_list)
		self.radius_std = stats.sem(self.r_list)

	def get_num_samples(self):
		return len(self.samples)

	def find_nD_list(self):
		sample_nD_list = []
		for sample in self.samples:
			sample_nD_list.append(sample.find_num_unmatched())	
		return sample_nD_list	

	def compute_radius(self):
		# generate a single sample (positions)
		# compute the distance list 
		# sort and find the r_critical, store to ensemble 

		Kappa = self.kappa
		N = self.n
		d = self.d
		shortcut_prob = self.shortcut_prob

		#Define the variables used:
		cdef double dom_size = 1
		cdef int i
		cdef int j
		cdef int k
		cdef double[:,:] positions
		cdef double dij
		cdef double dist

		#Define arrays to store the integer labels of the connected pairs:
		source= [ ]
		target= [ ]

		# N-list of d-lists
		# Randomly generate positions in a d-dimensional standard Normal distribution
		positions = np.random.multivariate_normal(np.zeros(d), np.eye(d), size=N)

		### FIRST RADIUS CALCULATION ###

		# compute distance list 
		d_list = sp.distance.pdist(positions)
		d_list_sorted = np.sort(d_list)
		
		# find the value required to get 
		# the required number of edges 
		# above the required threshold:
		# <k> = 2E/N
		num_edges = int((Kappa*N)/2.0)

		# return the critical radius 
		r_c = d_list_sorted[num_edges]	

		##################################	

		### MAKE A SAMPLE ###

		# number of nodes
		cdef int n = positions.shape[0]
		# number of dimensions
		cdef int q = positions.shape[1]

		# Define arrays to store the integer labels of the connected pairs:
		source= [ ]
		target= [ ]

		# N-list of d-lists
		# Randomly generate positions in a d-dimensional standard Normal distribution
		positions = np.random.multivariate_normal(np.zeros(d), np.eye(d), size=N)

		# form the sample with this radius 
		# make a distance matrix to save computation
		d_mat = np.zeros((N,N))
		for i in range(n) :
			d_mat[i][i] = 0
			for j in range(i+1,n) :
				dij = 0
				#Loop over number of dimensions
				for k in range(q):
					# Compute the absolute distance
					dist = abs( positions[i, k] - positions[j, k] )
					# Add to the total distance
					dij = dij + dist**2
				dij = dij**0.5
				d_mat[i][j] = dij**0.5

				# compute the connection probability:
				# returns 1 if within radius, shortcut_prob otherwise
				if dij < r_c :
					probability = 1.0
				else :
					probability = shortcut_prob

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
		# symmetric
		d_mat = d_mat + d_mat.T 

		sample = GaussianSample(
			    source=source,
				target=target,
				positions=np.asarray(positions),
				kappa=Kappa,
				n=N,
				d=d,
				shortcut_prob=shortcut_prob,
				)

		# change sample to LCC to recompute radius 
		sample.to_LCC()
		LCC_nodes = sample.LCC_vec

		### SECOND RADIUS CALCULATION ###

		# compute distance list 
		d_list = sp.distance.pdist(sample.positions)
		d_list_sorted = np.sort(d_list)
		
		# find the value required to get 
		# the required number of edges 
		# above the required threshold:
		# <k> = 2E/N
		num_edges = int((Kappa*sample.LCC_size)/2.0)

		# return the NEW critical radius 
		r_c = d_list_sorted[num_edges]	

		return r_c

	def create_sample(self):
		
		Kappa = self.kappa
		N = self.n
		d = self.d
		shortcut_prob = self.shortcut_prob

		#Define the variables used:
		cdef double dom_size = 1
		cdef int i
		cdef int j
		cdef int k
		cdef double r_c
		cdef double[:,:] positions
		cdef double dij
		cdef double dist

		#Define arrays to store the integer labels of the connected pairs:
		source= [ ]
		target= [ ]

		# N-list of d-lists
		# Randomly generate positions in a d-dimensional standard Normal distribution
		positions = np.random.multivariate_normal(np.zeros(d), np.eye(d), size=N)

		r_c = self.radius 

		# number of nodes
		cdef int n = positions.shape[0]
		# number of dimensions
		cdef int q = positions.shape[1]

		# form the sample with this radius 
		# make a distance matrix
		for i in range(n) :
			for j in range(i+1,n) :
				dij = 0
				#Loop over number of dimensions
				for k in range(q):
					# Compute the absolute distance
					dist = abs( positions[i, k] - positions[j, k] )
					# Add to the total distance
					dij = dij + dist**2
				dij = dij**0.5

				# compute the connection probability:
				# returns 1 if within radius, shortcut_prob otherwise
				if dij < r_c :
					probability = 1.0
				else :
					probability = shortcut_prob

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

		sample = GaussianSample(
			    source=source,
				target=target,
				positions=np.asarray(positions),
				kappa=Kappa,
				n=N,
				d=d,
				shortcut_prob=shortcut_prob,
				)

		sample.to_LCC()

		return sample 

	def to_LCC(self):
		for sample in self.samples:
			sample.to_LCC()

	def get_param_string(self):
		# param_string
		if abs(float(self.shortcut_prob)) == 0.0 :
			return "Gauss_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + "_num_r_" + str(len(self.r_list))
		else:
			return "Gauss_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + "_num_r_" + str(len(self.r_list)) + '_short_' + str(self.shortcut_prob)

	def to_dict(self):
		return {
			'kappa': self.kappa,
			'n': self.n,
			'd': self.d,
			'shortcut_prob': self.shortcut_prob,
			'samples': [s.to_dict() for s in self.samples]
		}

	@classmethod
	def from_dict(self, data):
		return GaussianEnsemble(
			kappa=data['kappa'],
			n=data['n'],
			d=data['d'],
			shortcut_prob=data['shortcut_prob'],
			samples=[GaussianSample.from_dict(s) for s in data['samples']],
		)

	def to_disk(self, folder='/Users/spencerw/Documents/DATA/rgg_samples'):
		filename = '%s/%s.pickle' % (folder, self.get_param_string())
		with open(filename, 'w') as fout:
			pickle.dump(self.to_dict(), fout)

	@classmethod
	def from_disk(self, path):
		with open(path, 'r') as fin:
			data = pickle.load(fin)
		return GaussianEnsemble.from_dict(data)

class GaussianSample(object):
	def __init__(self, source, target, positions, kappa, n, d, shortcut_prob):
		self.source = source
		self.target = target
		self.positions = np.asarray(positions)
		self.kappa = kappa
		self.n = n
		self.d = d
		self.shortcut_prob = shortcut_prob

		# # param_string
		# if abs(float(self.shortcut_prob)) == 0.0 :
		# 	self.param_string = "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + str(self.boundary)
		# else:
		# 	self.param_string = "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + str(self.boundary) + '_short_' + str(self.shortcut_prob)

		self.adjacency = None
		self.adjacency_dense = None
		self.bipartite_dict = None

		self.calculate_adjacency_representation()
		self.calculate_dictionary_representation()

	def get_param_string(self):
		# param_string
		if abs(float(self.shortcut_prob)) == 0.0 :
			return "Gauss_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) 
		else:	
			return "Gauss_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + '_short_' + str(self.shortcut_prob)

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

	def to_LCC(self):

		#### update the adjacency 
		# make it symmetric so the calculation works 
		symmetric = self.adjacency + self.adjacency.transpose(copy=True)

		# Till's code  
		index, LCC_vec = max(enumerate(Giant.Components(symmetric)), key = lambda tup: len(tup[1]))

		# take only nodes that are in the LCC_vec
		LCC_adjacency_dense = self.adjacency_dense[np.ix_(LCC_vec,LCC_vec)]
		LCC_adjacency = sparse.csr_matrix(LCC_adjacency_dense)

		# add sample object attributes 
		self.adjacency = LCC_adjacency
		self.adjacency_dense = LCC_adjacency_dense
		self.LCC_size = len(LCC_vec)
		self.LCC_vec = LCC_vec
		self.LCC_positions = self.positions[LCC_vec]

		#### update the dictionary 
		# convert the LCC adjacency into a dictionary 
		nums = xrange(self.LCC_size)
		nodelist = (str(num) + 's' for num in nums)
		D = {k: [] for k in nodelist}
		
		# row=out, col=in 
		# nonzero returns tuple of row,col indices of nonzeros 
		outgoing, incoming = (np.nonzero(LCC_adjacency_dense))
		for i in xrange(len(outgoing)):
			D[str(outgoing[i])+'s'].append(incoming[i])

		self.bipartite_dict = D

	def find_unmatched(self):

		graph = self.bipartite_dict.copy() 
		num_nodes = len(graph.keys())

		# list of node indices 0,1,2,...,N 
		unmatched = range(num_nodes)

		# find matched nodes
		# (they will be doubled (i-->j and j-->i) ) 
		matched = HopcroftKarp(graph).maximum_matching().values()
		# remove the matched nodes to leave the unmatched 
		for node in range(num_nodes):
			for match in matched:
				# if node exists in matching
				if node == match:
					# remove it 
					unmatched.remove(node)
		return unmatched  

	def find_num_unmatched(self):
		unmatched = self.find_unmatched()

		return len(unmatched)

	def plot_network(self, unmatched=None, size=20, label_nodes=False, save=True, LCC=False):
		
		if self.d != 2:
			raise ValueError('The graph is not 2D!')
		else:		
			if LCC == True: 
				File_Name = './plots/' + str(self.get_param_string()) + '_LCC'+ '_plot' + '.eps'
			else:
				File_Name = './plots/' + str(self.get_param_string()) + '_plot' + '.eps'
				
			plt.figure()
			plt.clf
			Gfig = open( File_Name , 'w' )
			
			# make a networkx graph object 
			graph = nx.DiGraph(self.adjacency_dense)

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
				if label_nodes: 
					nx.draw_networkx(graph, pos=posDict, with_labels=True, node_color=colorList, node_size=size, width = 0.3)
				else: 
					nx.draw_networkx(graph, pos=posDict, with_labels=False, node_color=colorList, node_size=size, width = 0.3)
			else:
				if label_nodes: 
					nx.draw_networkx(graph, pos=posDict, with_labels=True, node_size=20, width = 0.3)
				else: 
					nx.draw_networkx(graph, pos=posDict, with_labels=False, node_size=20, width = 0.3)
			
			if save == True: 
				plt.savefig( File_Name , format='eps', dpi=1000 )

	def mean_degree(self) :
		Degree_array = np.array( self.adjacency_dense.sum(axis=0) ) 
		# 2x because of directedness  
		return 2*np.mean( Degree_array[0] )

	def properties(self, *additional_properties):
		raise NotImplementedError()

	def to_dict(self):
		return {
			'kappa': self.kappa,
			'n': self.n,
			'd': self.d,
			'shortcut_prob': self.shortcut_prob,
			'source': self.source,
			'target': self.target,
			'positions': self.positions,
		}

	@classmethod
	def from_dict(cls, data):
		return GaussianSample(
			kappa=data['kappa'],
			n=data['n'],
			d=data['d'],
			shortcut_prob=data['shortcut_prob'],
			source=data['source'],
			target=data['target'],
			positions=data['positions'],
		)

	def to_disk(self, folder='/Users/spencerw/Documents/DATA/rgg_samples'):
		filename = '%s/%s.pickle' % (folder, self.get_param_string())
		with open(filename, 'w') as fout:
			pickle.dump(self.to_dict(), fout)

	@classmethod
	def from_disk(self, path):
		with open(path, 'r') as fin:
			data = pickle.load(fin)
		return GaussianSample.from_dict(data)

