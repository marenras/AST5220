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
	    print("Please write 'save' in command line to save the plots as tex-files.")
	    sys.exit()



theta = []
integrand = []
kcH = []

l_values = np.loadtxt('data/l_values.dat', unpack=True)
x, test = np.loadtxt('data/Sj_l.dat', unpack=True)
l, C_l = np.loadtxt('data/powerspectras/CMB_spectrum_original.dat', unpack=True)
l1, C_l1 = np.loadtxt('data/powerspectras/CMB_spectrum_best_fit.dat', unpack=True)
l2, C_l2 = np.loadtxt('data/powerspectras/CMB_spectrum_h2.dat', unpack=True)
l_real, CMB_real, ndC, pdC = np.loadtxt('data/COM_PowerSpect_CMB-TT-full_R3.01.txt', unpack=True)

for k in range(1,7):
    kcH_, theta_, integrand_ = np.loadtxt('data/theta_integrand_%d.dat' %(k), unpack=True)
    kcH.append(kcH_)
    theta.append(theta_)
    integrand.append(integrand_)


fig_CMB = plt.figure()
plt.errorbar(l_real, CMB_real, yerr=[ndC,pdC], color='lightblue', label='Observed', zorder=1)
plt.plot(l, C_l*5775/max(C_l), label='Default')
plt.plot(l1, C_l1*5775/max(C_l1), label='Computed1')
plt.plot(l2,C_l2*5775/max(C_l2), label='Computed2')
plt.ylabel(r'$C_l$')
plt.xlabel(r'$l$')
plt.legend(loc='best')
plt.xlim(l[0], l[-1])
plt.grid()
plt.show()




#plt.plot(x, test)
#plt.show()


#fig_theta = plt.figure()
#for i in range(6):
#    plt.plot(kcH[i], theta[i], label=r'l = %d' %(l_values[i]))
#plt.ylabel(r'$\Theta$')
#plt.xlabel(r'$ckH$')
#plt.legend(loc='best')
#plt.grid()
#plt.show()
