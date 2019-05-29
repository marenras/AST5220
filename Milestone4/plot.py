import numpy as np
import matplotlib.pyplot as plt
import sys as sys

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 13}

plt.rc('font', **font)


save = False
if len(sys.argv) > 1:
    if sys.argv[1] in ['save', 'Save']:
        save = True
    else:
	    print("Please write 'save' in command line to save the plots.")
	    sys.exit()

# Reading in l-values for the tranfer function and the integrand
l_values = np.loadtxt('data/l_values.dat', unpack=True)

# Reading in power spectra values with the best parameters
l, C_l = np.loadtxt('data/powerspectras/CMB_spectrum_best_fit.dat', unpack=True)

# Reading in power spectra observed data
l_real, CMB_real, ndC, pdC = np.loadtxt('data/COM_PowerSpect_CMB-TT-full_R3.01.txt', unpack=True)


# Reading in tranfer function and integrand for 6 values of k
theta = []
integrand = []
kcH = []

for k in range(1,7):
    kcH_, theta_, integrand_ = np.loadtxt('data/theta_integrand_%d.dat' %(k), unpack=True)
    kcH.append(kcH_)
    theta.append(theta_)
    integrand.append(integrand_)


fig_CMB = plt.figure()
plt.errorbar(l_real, CMB_real, yerr=[ndC,pdC], color='lightblue', label='Observed', zorder=1)
plt.plot(l, C_l*5775/max(C_l), label='Best fit')
plt.ylabel(r'$l(l+1) C_l / 2 \pi$ [$\mu$K$^2$]')
plt.xlabel(r'$l$')
plt.legend(loc='best')
plt.xlim(l[0], l[-1])
plt.grid()
plt.show()

fig_theta = plt.figure()
for i in range(6):
   plt.plot(kcH[i], theta[i], label=r'l = %d' %(l_values[i]))
plt.ylabel(r'$\Theta_l$')
plt.xlabel(r'$ck/H_0$')
plt.legend(loc='best')
plt.grid()
plt.show()

fig_integrand = plt.figure()
for i in range(6):
   plt.plot(kcH[i], integrand[i], label=r'l = %d' %(l_values[i]))
plt.ylabel(r'$\Theta_l^2/k$')
plt.xlabel(r'$ck/H_0$')
plt.legend(loc='best')
plt.grid()
plt.show()

if save:
    fig_CMB.savefig("figures/CMB_best_fit.pdf", bbox_inches='tight')
