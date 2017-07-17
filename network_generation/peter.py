import pickle
import numpy as np
from scipy import sparse

class PeterRGGSample(object):
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
		nodelist = xrange(self.n)
		# ['1s', '2s', '3s', ...] for the lefthand side of the bipartite
		nodelist = (str(num) + 's' for num in nodelist)

		# generate an empty dict of lists for targets
		D = {k: [] for k in nodelist}

		# go through and add targets to source key -- target sets
		for idx, val in enumerate(self.source):
			D[str(val)+'s'].append(self.target[idx])

		self.bipartite_dict = D

	def param_string(self):
		# Param_String
		if abs(float(self.shortcut_prob)) == 0.0 :
			return "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + str(self.boundary)

		return "RGG_K_" + str(self.kappa) + "_N_" + str(self.n) + "_d_" + str(self.d) + "_BC_" + str(self.boundary) + '_short_' + str(self.shortcut_prob)

	def find_unmatched(self):
		raise NotImplementedError()

	def plot_network(self, unmatched=None, sizea=20, label=False):
		raise NotImplementedError()

	def mean_degree(self) :
		raise NotImplementedError()

	def properties(self, *additional_properties):
		raise NotImplementedError()

        # XXX These are the new methods



        def to_dict(self):
            return {
                'kappa': self.kappa,
                'n': self.n,
                'd': self.d,
                'boundary': self.boundary,
                'shortcut_prob': self.shortcut_prob,
                # TODO: change the numpy arrays into lists
                'source': self.source,
                'target': self.target,
                'positions': self.positions,
            }

        @classmethod
        def from_dict(cls, data):
            return PeterRGGSample(
                kappa=data['kappa'],
                n=data['n'],
                d=data['d'],
                boundary=data['boundary'],
                shortcut_prob=data['shortcut_prob'],
                # TODO: change the following back into numpy arrays
                source=data['source'],
                target=data['target'],
                positions=data['positions'],
            )
