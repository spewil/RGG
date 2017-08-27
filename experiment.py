# a run of ensembles 

import scipy
import numpy as np
import network_generation.generation as ng
import os 
import networkx as nx 
from scipy import stats

class RGGExperiment(object):

	def __init__(self, kappa_range, n, d, shortcut_prob=0, boundary='s', ensembles=[],num_radii=1):
		
		# params 
		self.kappa_range = kappa_range
		self.n = n
		self.d = d
		self.shortcut_prob = shortcut_prob
		self.boundary = boundary	
		self.num_radii = num_radii

		# populate ensembles 
		self.ensembles = self.add_ensembles()

	def add_ensembles(self):
		ensembles = []
		for kappa in self.kappa_range: 
			# import data 
			RGG = ng.RGGEnsemble(kappa,self.n,self.d,boundary=self.boundary, num_radii=self.num_radii)
			try: 
				RGG = ng.RGGEnsemble.from_disk('/Users/spencerw/Documents/DATA/rgg_samples/' + RGG.get_param_string()+'.pickle')
			except:
				raise ValueError('Data does not exist for: ' + str(RGG.get_param_string()))
			ensembles.append(RGG)
		return ensembles	

	def find_nD_stats(self):
		mean_nD_list = []
		std_nD_list = []
		stddev_nD_list = []
		# find unmatched nodes 
		for ensemble in self.ensembles:
			sample_nD_list = []	
			for sample in ensemble.samples:
				sample_nD_list.append(sample.find_num_unmatched())
			mean_nD_list.append(np.mean(sample_nD_list)/self.n)
			std_nD_list.append(stats.sem([sample/self.n for sample in sample_nD_list]))
			stddev_nD_list.append(np.std(sample_nD_list)/self.n)
		return mean_nD_list, std_nD_list, stddev_nD_list

	def find_degree_stats(self):
		mean_degree_list = []
		std_degree_list = []
		stddev_degree_list = []
		# find unmatched nodes 
		for ensemble in self.ensembles:
			sample_degree_list = []
			for sample in ensemble.samples:
				# mean degree of the graph 
				sample_degree_list.append(sample.mean_degree())
			# mean of the list of mean degree  
			mean_degree_list.append(np.mean(sample_degree_list))
			std_degree_list.append(stats.sem(sample_degree_list))
			stddev_degree_list.append(np.std(sample_degree_list))
		return mean_degree_list, std_degree_list, stddev_degree_list

	def find_n_stats(self):
		mean_n_list = []
		std_n_list = []
		stddev_n_list = []
		# find unmatched nodes 
		for ensemble in self.ensembles:
			sample_n_list = []
			for sample in ensemble.samples:
				# mean degree of the graph 
				try:
					sample_n_list.append(sample.LCC_size)
				except AttributeError:
					sample_n_list.append(sample.n)
			# mean of the list of mean degree  
			mean_n_list.append(np.mean(sample_n_list))
			std_n_list.append(stats.sem(sample_n_list))
			stddev_n_list.append(np.std(sample_n_list))
		return mean_n_list, std_n_list, stddev_n_list

	def find_scaling_params(self):
		# curve fit the <k> vs <nD> exponential relationship 
		# f(t) = a*np.exp(b*t)
		# returns [[a,b],[covariance matrix]] 
		scaling_params = scipy.optimize.curve_fit(lambda t,a,b: a*np.exp(b*t),  self.mean_degree_list,  self.mean_nD_list,  p0=(1, -0.5))
		return scaling_params

	def to_LCC(self):
		for ensemble in self.ensembles:
			ensemble.to_LCC()

		# update mean nD list, mean k list, mean n list 
		# self.mean_nD_list = self.find_nD_stats()[0]
		# self.std_nD_list = self.find_nD_stats()[1]
		# self.mean_degree_list = self.find_degree_stats()[0]
		# self.std_degree_list = self.find_degree_stats()[1]
		# self.mean_n_list = self.find_n_stats()[0]
		# self.std_n_list = self.find_n_stats()[1]

class ERExperiment(object):

	def __init__(self, kappa_range, n, ensembles=[]):
		
		# params 
		self.kappa_range = kappa_range
		self.n = n

		# populate ensembles 
		self.ensembles = self.add_ensembles()

	def add_ensembles(self):
		ensembles = []
		for kappa in self.kappa_range: 	
			# import data 
			ER = ng.EREnsemble(kappa,self.n)
			try: 
				ER = ng.EREnsemble.from_disk('/Users/spencerw/Documents/DATA/er_samples/'+ER.get_param_string()+'.pickle')
			except:
				raise ValueError('Data does not exist')
			ensembles.append(ER)
		return ensembles

	def find_nD_stats(self):
		mean_nD_list = []
		std_nD_list = []
		stddev_nD_list = []
		# find unmatched nodes 
		for ensemble in self.ensembles:
			sample_nD_list = []	
			for sample in ensemble.samples:
				sample_nD_list.append(sample.find_num_unmatched())
			mean_nD_list.append(np.mean(sample_nD_list)/self.n)
			std_nD_list.append(stats.sem([sample/self.n for sample in sample_nD_list]))
			stddev_nD_list.append(np.std(sample_nD_list)/self.n)
		return mean_nD_list, std_nD_list, stddev_nD_list

	def find_degree_stats(self):
		mean_degree_list = []
		std_degree_list = []
		stddev_degree_list = []
		# find unmatched nodes 
		for ensemble in self.ensembles:
			sample_degree_list = []
			for sample in ensemble.samples:
				# mean degree of the graph 
				sample_degree_list.append(sample.mean_degree())
			# mean of the list of mean degree  
			mean_degree_list.append(np.mean(sample_degree_list))
			std_degree_list.append(stats.sem(sample_degree_list))
			stddev_degree_list.append(np.std(sample_degree_list))
		return mean_degree_list, std_degree_list, stddev_degree_list

	def find_scaling_params(self):
		# curve fit the <k> vs <nD> exponential relationship 
		# f(t) = a*np.exp(b*t)
		# returns [[a,b],[covariance matrix]] 
		scaling_params = scipy.optimize.curve_fit(lambda t,a,b: a*np.exp(b*t),  self.mean_degree_list,  self.mean_nD_list,  p0=(1, -0.5))
		return scaling_params

	def find_n_stats(self):
		mean_n_list = []
		std_n_list = []
		stddev_n_list = []
		# find unmatched nodes 
		for ensemble in self.ensembles:
			sample_n_list = []
			for sample in ensemble.samples:
				# mean degree of the graph 
				try:
					sample_n_list.append(sample.LCC_size)
				except AttributeError:
					sample_n_list.append(sample.n)
			# mean of the list of mean degree  
			mean_n_list.append(np.mean(sample_n_list))
			std_n_list.append(stats.sem(sample_n_list))
			stddev_n_list.append(np.std(sample_n_list))
		return mean_n_list, std_n_list, stddev_n_list
	
	def to_LCC(self):
		for ensemble in self.ensembles:
			ensemble.to_LCC()

	def find_gamma(self):
		return self.scaling_params[0][1]


