import numpy as np 
import experiment as ex
import matplotlib.pyplot as plt
import pickle

def find_ER_nD(kappa, plot=False):
        
    # make sure this is a float!
    z0 = float(kappa/2)
    w1irange = np.arange(0,1,.01)
    w1o = lambda x: np.exp(-z0*np.exp(-z0*x))
    w1i = lambda y: -(1/z0)*np.log((1/z0)*np.log(1/y))

    if plot:
        plt.figure()
        plt.plot(w1irange,[w1o(i) for i in w1irange])
        plt.plot(w1irange,w1irange,'k-')
        
        w1o = lambda x: np.exp(-z0*np.exp(-z0*x))
    w1i = lambda y: -(1/z0)*np.log((1/z0)*np.log(1/y))
    
    # BABY NEWTON METHOD 
    xold = 0.5
    yold = 0.5
    check = 100
    while abs(check)>0.001:

        # move up
        ynew = w1o(xold)
        # move right 
        xnew = ynew
        
        if plot: 
            plt.plot((xold,xold),(yold,ynew),'r-')
            plt.plot((xold,xnew),(ynew,ynew),'r-')

        # compute check 
        check = xold - xnew 

        # update for next iteration
        xold = xnew 
        yold = ynew

    w1 = xnew 
    w2 = 1 - np.exp(-z0*w1)
    nD = w1 - w2 + z0*w1*(1-w2)

    return nD

###########################################################

kappa_range = kappa_range = [1,3,5,10,20]
d_range = [2,3,6,9]
n = 1000

## SOLID 

plt.figure()
for d in d_range:
	print d 
	RGGE = ex.RGGExperiment(kappa_range, n, d, shortcut_prob=0, boundary='s',find_scaling=False)
	plt.errorbar(RGGE.mean_degree_list, RGGE.mean_nD_list, fmt='o', markersize=5, yerr=RGGE.std_nD_list, xerr=RGGE.std_degree_list)

plt.xlabel('$\\langle{k}\\rangle$')
plt.ylabel('$n_D$')
plt.title('$n_D$ vs $\\langle{k}\\rangle$ over $d$ for $N=1000$')
legend_string = ['d='+str(d) for d in d_range]
legend_string.reverse()

# ER THEORY
nD_list = []
for kappa in kappa_range:
    nD_list.append(find_ER_nD(kappa))
plt.plot(kappa_range,nD_list,'b-')   

legend_string.append('ER Theory')
plt.legend(legend_string)
filename = './plots/nD_vs_k_n' + str(n) + '_b_s.eps' 
plt.savefig(filename)
plt.show()


# ## PERIODIC 

# plt.figure(2)
# for d in d_range:
# 	RGGE = ex.RGGExperiment(kappa_range, n, d, shortcut_prob=0, boundary='p',find_scaling=False)
# 	plt.errorbar(RGGE.mean_degree_list, RGGE.mean_nD_list, yerr=RGGE.std_nD_list, xerr=RGGE.std_degree_list)

# x = np.arange(1,20,.1)
# plt.plot(x,[np.exp(-xx/2) for xx in x],'k--')

# plt.xlabel('$\\langle{k}\\rangle$')
# plt.ylabel('$n_D$')
# plt.title('$n_D$ vs $\\langle{k}\\rangle$ over $d$ for $N=1000$')
# legend_string = ['d='+str(d) for d in d_range]
# legend_string.reverse()
# legend_string.append('ER Theory')
# legend_string.reverse()
# plt.legend(legend_string)
# plt.savefig('./plots/nD_vs_k_n100_p.eps')


##################################
