# generate data  

import network_generation.generation as ng

# kappa_range = kappa_range = [1,3,5,10,15,20]
# n_range = [100,1000]
# d_range =  [4,5,6,7,8,10] #[2,3,6,9]	
# boundaries = ['s','p','g']
# num_samples = 100 

kappa_range = kappa_range = [1,3,5,10,15,20]
n_range = [10000]
d_range = [2,9,25] #[2,3,6,9,15,25,40]
boundaries = ['s'] #['s','p','g']
num_samples = 10

print 'running ' + str(len(n_range)*len(d_range)*len(boundaries)*len(kappa_range)) + ' generations'
i=1

## Generate RGG Data ##
for kappa in kappa_range:
    for n in n_range:
        for d in d_range:
            for boundary in boundaries: 
                RGG = ng.RGGEnsemble(kappa, n, d, boundary=boundary, num_radii=3)
                RGG.generate_samples(n=num_samples)
                # data to disk 
                RGG.to_disk()
                print 'saved run ' + str(i)
                i+=1 

# GENERATE ER DATA 
print 'running ' + str(len(n_range)*len(kappa_range)) + ' generations'
i=1
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