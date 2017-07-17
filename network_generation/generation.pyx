# coding: utf-8
from random import random
from scipy import sparse
from scipy import special
import numpy as np
import pickle
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
		self.kappa = kappa
		self.n = n
		self.d = d
		self.shortcut_prob = shortcut_prob
		self.boundary = boundary
		self.samples = samples or []

	def create_sample(self):
		Kappa = self.kappa
		N = self.n
		d = self.d
		shortcut_prob = self.shortcut_prob
		Boundary = self.boundary

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

		return RGGSample(
			source=np.asarray(S),
			target=T,
			positions=np.asarray(positions),
			kappa=Kappa,
			n=N,
			d=d,
			boundary=Boundary,
			shortcut_prob=shortcut_prob,
		)

	def param_string(self):
		# Param_String
		if abs(float(self.shortcut_prob)) == 0.0 :
			return "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + str(self.boundary)

		return "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + str(self.boundary) + '_short_' + str(self.shortcut_prob)

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

	def to_disk(self, folder='rgg_samples'):
		filename = './%s/%s.pickle' % (folder, self.param_string())
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
		raise NotImplementedError()

	def plot_network(self, unmatched=None, sizea=20, label=False):
		raise NotImplementedError()

	def mean_degree(self) :
		raise NotImplementedError()

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
		Generates sparse adjacency matrix for an Erdos Renyi Random Graph


		Parameters
		------------

		Kappa : double

		Mean degree parameter:

		N : int

		Number of nodes in the network
		"""
		Kappa = self.kappa
		N = self.n

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
