# save gamma data 

import numpy as np 
import experiment as ex
import pickle

kappa_range = kappa_range = [1,3,5,10,15,20]
n_range = range(100,1001,100)
d_range = range(2,11,1)

## SOLID 

gamma_matrix = []
gamma_std_matrix = []
print 'number of imports: ' + str(len(n_range)*len(d_range))
i = 1 
for n in n_range:
    gamma_d_list = []
    for d in d_range: 
        RGGE = ex.RGGExperiment(kappa_range, n, d, shortcut_prob=0, boundary='s')
        gamma_d_list.append(RGGE.gamma)
        gamma_std = np.sqrt(RGGE.scaling_params[1][1][1])
        gamma_d_std_list.append(gamma_std)
        print 'import ' + str(i) + ' done'
        i += 1
    gamma_matrix.append(gamma_d_list)
    gamma_std_matrix.append(gamma_d_std_list)
    
gamma_matrix = np.array(gamma_matrix)
gamma_std_matrix = np.array(gamma_std_matrix)

folder = 'scaling_data'
filename = './%s/%s.pickle' % (folder, 'RGG_gamma_matrix_s')
with open(filename, 'w') as fout:
	pickle.dump(gamma_matrix, fout)
filename = './%s/%s.pickle' % (folder, 'RGG_gamma__std_matrix_s')
with open(filename, 'w') as fout:
	pickle.dump(gamma_matrix, fout)

## PERIODIC

gamma_matrix = []
gamma_std_matrix = []

print 'number of imports: ' + str(len(n_range)*len(d_range))
i = 1 
for n in n_range:
    gamma_d_list = []
    for d in d_range: 
        RGGE = ex.RGGExperiment(kappa_range, n, d, shortcut_prob=0, boundary='p')
        gamma_d_list.append(RGGE.gamma)
        gamma_std = np.sqrt(RGGE.scaling_params[1][1][1])
        gamma_d_std_list.append(gamma_std)
        print 'import ' + str(i) + ' done'
        i += 1
    gamma_matrix.append(gamma_d_list)
    gamma_std_matrix.append(gamma_d_std_list)
    
gamma_matrix = np.array(gamma_matrix)
gamma_std_matrix = np.array(gamma_std_matrix)

filename = './%s/%s.pickle' % (folder, 'RGG_gamma_matrix_p')
with open(filename, 'w') as fout:
	pickle.dump(gamma_matrix, fout)

filename = './%s/%s.pickle' % (folder, 'RGG_gamma__std_matrix_p')
with open(filename, 'w') as fout:
	pickle.dump(gamma_std_matrix, fout)






