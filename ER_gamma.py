# save gamma data 

import numpy as np 
import experiment as ex
import matplotlib.pyplot as plt
import pickle

kappa_range = kappa_range = [1,3,5,10,15,20]
n_range = range(100,1001,100)

## SOLID 

gamma_list = []
gamma_std_list = []

print 'number of imports: ' + str(len(n_range))
i = 1 
for n in n_range:
    ERE = ex.ERExperiment(kappa_range, n)
    print 'import ' + str(i) + ' done'
    i += 1
    gamma_list.append(ERE.gamma)
    gamma_std_list.append(ERE.scaling_params[1][1][1])
    
gamma_list = np.array(gamma_list)
gamma_std_list = np.array(gamma_std_list)

plt.clf()
plt.errorbar(n_range,gamma_list,yerr=gamma_std_list)
plt.xlabel('Graph Size $N$	')
plt.ylabel('Scaling $\\gamma$')
plt.tight_layout()
plt.savefig('./plots/ER_gamma.png')