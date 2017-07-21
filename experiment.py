# a run of ensembles 

import scipy
import numpy as np
import network_generation.generation as ng

class RGGExperiment(object):

	def __init__(self, kappa_range, n, d, shortcut_prob=0, boundary='s', ensembles=[], find_scaling = True):
		
		# params 
		self.kappa_range = kappa_range
		self.n = n
		self.d = d
		self.shortcut_prob = shortcut_prob
		self.boundary = boundary	

		# populate ensembles 
		self.ensembles = self.add_ensembles()

		# find mean nD list, mean k list
		self.mean_nD_list = self.find_nD_stats()[0]
		self.std_nD_list = self.find_nD_stats()[1]
		self.mean_degree_list = self.find_degree_stats()[0]
		self.std_degree_list = self.find_degree_stats()[1]

		if find_scaling:
			# find single parameter defining the Experiment 
			self.scaling_params = self.find_scaling_params()
			self.gamma = self.find_gamma()


	def add_ensembles(self):

		ensembles = []

		for kappa in self.kappa_range: 
			
			# import data 
			RGG = ng.RGGEnsemble(kappa,self.n,self.d,boundary=self.boundary)
			try: 
				RGG = ng.RGGEnsemble.from_disk('rgg_samples/'+RGG.get_param_string()+'.pickle')
			except:
				raise ValueError('Data does not exist')
			
			ensembles.append(RGG)

		return ensembles

	def find_nD_stats(self):
		mean_nD_list = []
		std_nD_list = []

		# find unmatched nodes 
		for ensemble in self.ensembles:
			sample_nD_list = []
			for sample in ensemble.samples:
				sample_nD_list.append(sample.find_num_unmatched())
			mean_nD_list.append(np.mean(sample_nD_list)/self.n)
			std_nD_list.append(np.std(sample_nD_list)/self.n)

		return mean_nD_list, std_nD_list

	def find_degree_stats(self):

		mean_degree_list = []
		std_degree_list = []
		# find unmatched nodes 
		for ensemble in self.ensembles:
			sample_degree_list = []
			for sample in ensemble.samples:
				# mean degree of the graph 
				sample_degree_list.append(sample.mean_degree())
			# mean of the list of mean degree  
			mean_degree_list.append(np.mean(sample_degree_list))
			std_degree_list.append(np.std(sample_degree_list))

		return mean_degree_list, std_degree_list

	def find_scaling_params(self):
		# curve fit the <k> vs <nD> exponential relationship 
		# f(t) = a*np.exp(b*t)
		# returns [[a,b],[covariance matrix]] 
		scaling_params = scipy.optimize.curve_fit(lambda t,a,b: a*np.exp(b*t),  self.mean_degree_list,  self.mean_nD_list,  p0=(1, -0.5))
		return scaling_params

	def find_gamma(self):
		return self.scaling_params[0][1]

class ERExperiment(object):

	def __init__(self, kappa_range, n, ensembles=[]):
		
		# params 
		self.kappa_range = kappa_range
		self.n = n

		# populate ensembles 
		self.ensembles = self.add_ensembles()

		# find mean nD list, mean k list
		self.mean_nD_list = self.find_nD_stats()[0]
		self.std_nD_list = self.find_nD_stats()[1]
		self.mean_degree_list = self.find_degree_stats()[0]
		self.std_nD_list = self.find_degree_stats()[1]

		# find single parameter defining the Experiment 
		self.scaling_params = self.find_scaling_params()
		self.gamma = self.find_gamma()


	def add_ensembles(self):

		ensembles = []

		for kappa in self.kappa_range: 
			
			# import data 
			ER = ng.EREnsemble(kappa,self.n)
			try: 
				ER = ng.EREnsemble.from_disk('er_samples/'+ER.get_param_string()+'.pickle')
			except:
				raise ValueError('Data does not exist')
			
			ensembles.append(ER)

		return ensembles

	def find_nD_stats(self):
		mean_nD_list = []
		std_nD_list = []

		# find unmatched nodes 
		for ensemble in self.ensembles:
			sample_nD_list = []
			for sample in ensemble.samples:
				sample_nD_list.append(sample.find_num_unmatched())
			mean_nD_list.append(np.mean(sample_nD_list)/self.n)
			std_nD_list.append(np.std(sample_nD_list)/self.n)

		return mean_nD_list, std_nD_list

	def find_degree_stats(self):

		mean_degree_list = []
		std_degree_list = []
		# find unmatched nodes 
		for ensemble in self.ensembles:
			sample_degree_list = []
			for sample in ensemble.samples:
				# mean degree of the graph 
				sample_degree_list.append(sample.mean_degree())
			# mean of the list of mean degree  
			mean_degree_list.append(np.mean(sample_degree_list))
			std_degree_list.append(np.std(sample_degree_list))

		return mean_degree_list, std_degree_list

	def find_scaling_params(self):
		# curve fit the <k> vs <nD> exponential relationship 
		# f(t) = a*np.exp(b*t)
		# returns [[a,b],[covariance matrix]] 
		scaling_params = scipy.optimize.curve_fit(lambda t,a,b: a*np.exp(b*t),  self.mean_degree_list,  self.mean_nD_list,  p0=(1, -0.5))
		return scaling_params

	def find_gamma(self):
		return self.scaling_params[0][1]


