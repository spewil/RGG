import numpy as np
import matplotlib.pyplot as plt 
import experiment as ex
from network_generation import generation as ng

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make an RGG Sample and draw the network
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print 'Generating RGG Data'

kappa = 5
n = 100
d = 2
boundary = 's'
num_samples = 1

RGG = ng.RGGEnsemble(kappa, n, d, shortcut_prob=0, boundary=boundary, num_radii=3)
RGG.generate_samples(n=num_samples)
RGG.to_LCC()
sample = RGG.samples[0]

sample.plot_network(unmatched=sample.find_unmatched(),label_nodes=True,show=True)

print 'Mean Degree: ',sample.mean_degree()
print 'Size of LCC: ',sample.LCC_size
print 'Driver Nodes: ',sample.find_unmatched()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make an RGG Experiment and plot the data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print 'Creating RGG Experiment'

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

kappa_range = [1,3,5,10,15,20]
n = 100 # change to 1000 for interpretable results
num_radii = 10
boundary = 's' #(or 'p' or 'g')
d_range = [2,3,6,9,15,25]

colors = ['ro','bo','go','ko','mo','co']
i=0

# if the following returns "Data Not Found", 
# Generate the data using the script "generate_RGG_nD_k_data.py"
# using the required parameters.

handles = []
labels = []
for d in d_range:
    RGG = ex.RGGExperiment(kappa_range, n, d, shortcut_prob=0, boundary=boundary, num_radii=num_radii)
    RGG.to_LCC()
    for ensemble in RGG.ensembles:
        for sample in ensemble.samples:
            data, = ax.plot(sample.mean_degree(),float(sample.find_num_unmatched())/sample.LCC_size,colors[i],markersize=2)
    handles.append(data)
    labels.append('d = '+str(d))
    print 'dimension ',str(d),' is done'
    
    i += 1

ax.set_xlabel('Mean Sample Degree $\\langle{k}\\rangle$')
ax.set_ylabel('Fraction of Driver Nodes $n_D$')
ax.legend(handles,labels)
fig.tight_layout()
# fig.savefig('./plots/scatter_N_' + str(n) + '.eps', dpi=1000)

plt.show()
