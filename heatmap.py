import matplotlib.pyplot as plt
import numpy as np
import pickle 

# # TESTING 
# # increasing y coordinate is higher index  
# test_matrix = np.array([[(d*5)+n for d in range(5)] for n in range(5)])
# print test_matrix
# plt.clf()
# plt.imshow(test_matrix, origin='lower')
# plt.xlabel('d')
# plt.ylabel('n')
# plt.colorbar()
# plt.show()

kappa_range = kappa_range = [1,3,5,10,15,20]
n_range = range(100,1001,100)
d_range = range(2,11,1)

## SOLID 

path = 'RGG_gamma_matrix_s.pickle'
with open(path, 'r') as fin:
	gamma_matrix = pickle.load(fin)

# gamma matrix --> 

extent = [d_range[0], d_range[-1], n_range[0]/100, n_range[-1]/100]
plt.clf()
plt.imshow(gamma_matrix, extent=extent, interpolation='nearest', origin='lower')
plt.colorbar()
plt.xlabel('Dimension $d$')
plt.ylabel('Graph Size $N$ [Hundreds]')
plt.savefig('./plots/heatmap_RGG_solid.png')

## PERIODIC 

path = 'RGG_gamma_matrix_p.pickle'
with open(path, 'r') as fin:
	gamma_matrix = pickle.load(fin)

extent = [d_range[0], d_range[-1], n_range[0]/100, n_range[-1]/100]
plt.figure()
plt.imshow(gamma_matrix, extent=extent, interpolation='nearest', origin='lower')
plt.colorbar()
plt.xlabel('Dimension $d$')
plt.ylabel('Graph Size $N$ [Hundreds]')
plt.savefig('./plots/heatmap_RGG_periodic.png')
