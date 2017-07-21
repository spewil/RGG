import numpy as np 
import experiment as ex
import matplotlib.pyplot as plt
import pickle

kappa_range = kappa_range = [1,3,5,10,15,20]
n_range = range(100,1001,100)
d_range = range(2,11,1)

## SOLID 

n = 1000
plt.figure()
for d in d_range:
	RGGE = ex.RGGExperiment(kappa_range, n, d, shortcut_prob=0, boundary='s',find_scaling=False)
	plt.errorbar(RGGE.mean_degree_list, RGGE.mean_nD_list, yerr=RGGE.std_nD_list, xerr=RGGE.std_degree_list)

x = np.arange(1,20,.1)
plt.plot(x,[np.exp(-xx/2) for xx in x],'k--')

plt.xlabel('$\\langle{k}\\rangle$')
plt.ylabel('$n_D$')
plt.title('$n_D$ vs $\\langle{k}\\rangle$ over $d$ for $N=1000$')
legend_string = ['d='+str(d) for d in d_range]
legend_string.reverse()
legend_string.append('ER Theory')
legend_string.reverse()
plt.legend(legend_string)
plt.yscale('log')
plt.savefig('./plots/nD_vs_k_n1000_s_log.eps')

d = 7
plt.figure()
for n in n_range:
	RGGE = ex.RGGExperiment(kappa_range, n, d, shortcut_prob=0, boundary='s',find_scaling=False)
	plt.errorbar(RGGE.mean_degree_list, RGGE.mean_nD_list, yerr=RGGE.std_nD_list, xerr=RGGE.std_degree_list)

x = np.arange(1,10,.1)
plt.plot(x,[np.exp(-xx/2) for xx in x],'k--')

plt.xlabel('$\\langle{k}\\rangle$')
plt.ylabel('$n_D$')
plt.title('$n_D$ vs $\\langle{k}\\rangle$ over $N$ for $d=7$')
legend_string = ['N='+str(n) for n in n_range]
legend_string.reverse()
legend_string.append('ER Theory')
legend_string.reverse()
plt.legend(legend_string)
plt.yscale('log')
plt.savefig('./plots/nD_vs_k_d7_s_log.eps')

## PERIODIC 

n = 1000
plt.figure()
for d in d_range:
	RGGE = ex.RGGExperiment(kappa_range, n, d, shortcut_prob=0, boundary='p',find_scaling=False)
	plt.errorbar(RGGE.mean_degree_list, RGGE.mean_nD_list, yerr=RGGE.std_nD_list, xerr=RGGE.std_degree_list)

x = np.arange(1,20,.1)
plt.plot(x,[np.exp(-xx/2) for xx in x],'k--')

plt.xlabel('$\\langle{k}\\rangle$')
plt.ylabel('$n_D$')
plt.title('$n_D$ vs $\\langle{k}\\rangle$ over $d$ for $N=1000$')
legend_string = ['d='+str(d) for d in d_range]
legend_string.reverse()
legend_string.append('ER Theory')
legend_string.reverse()
plt.legend(legend_string)
plt.yscale('log')
plt.savefig('./plots/nD_vs_k_n1000_p_log.eps')

d = 7
plt.figure()
for n in n_range:
	RGGE = ex.RGGExperiment(kappa_range, n, d, shortcut_prob=0, boundary='p',find_scaling=False)
	plt.errorbar(RGGE.mean_degree_list, RGGE.mean_nD_list, yerr=RGGE.std_nD_list, xerr=RGGE.std_degree_list)

x = np.arange(1,20,.1)
plt.plot(x,[np.exp(-xx/2) for xx in x],'k--')

plt.xlabel('$\\langle{k}\\rangle$')
plt.ylabel('$n_D$')
plt.title('$n_D$ vs $\\langle{k}\\rangle$ over $N$ for $d=7$')
legend_string = ['N='+str(n) for n in n_range]
legend_string.reverse()
legend_string.append('ER Theory')
legend_string.reverse()
plt.legend(legend_string)
plt.yscale('log')
plt.savefig('./plots/nD_vs_k_d7_p_log.eps')
