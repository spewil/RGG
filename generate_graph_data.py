# generate data  

import network_generation.generation as ng

kappa_range = kappa_range = [1,3,5,10,15,20]
n_range = range(100,1001,100)
d_range = range(2,11,1)	
boundaries = ['s','p']

num_samples = 10

print 'running ' + str(len(n_range)*len(d_range)*len(boundaries)*len(kappa_range)) + ' generations'
i=1

## RGGs ##
for kappa in kappa_range:
	for n in n_range:
		for d in d_range:
			for boundary in boundaries: 
				# make an RGG Ensemble with solid boundarys
				RGG = ng.RGGEnsemble(kappa,n,d,boundary=boundary)

				# generate samples 
				# sample data is stored in the object 
				RGG.generate_samples(n=num_samples)

				# data to disk 
				RGG.to_disk()

				print 'saved run ' + str(i)
				i+=1 

print 'running ' + str(len(n_range)*len(kappa_range)) + ' generations'
i=1				

## ERs ##
for kappa in kappa_range:
	for n in n_range:
		# make an ER Ensemble
		ER = ng.EREnsemble(kappa,n)

		# generate samples 
		# sample data is stored in the object 
		ER.generate_samples(n=num_samples)

		# data to disk 
		ER.to_disk()

		print 'saved run ' + str(i)
		i+=1 

