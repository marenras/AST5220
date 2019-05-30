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


# Reading in power spectra values
l, C_l = np.loadtxt('data/CMB_spectrum_original.dat', unpack=True)

# Changing n
l_n1, C_l_n1 = np.loadtxt('data/CMB_spectrum_n1.dat', unpack=True)
l_n2, C_l_n2 = np.loadtxt('data/CMB_spectrum_n2.dat', unpack=True)

# Changing omega_b
l_b1, C_l_b1 = np.loadtxt('data/CMB_spectrum_b1.dat', unpack=True)
l_b2, C_l_b2 = np.loadtxt('data/CMB_spectrum_b2.dat', unpack=True)

# Changing omega_m
l_m1, C_l_m1 = np.loadtxt('data/CMB_spectrum_m1.dat', unpack=True)
l_m2, C_l_m2 = np.loadtxt('data/CMB_spectrum_m2.dat', unpack=True)

# Changing h
l_h1, C_l_h1 = np.loadtxt('data/CMB_spectrum_h1.dat', unpack=True)
l_h2, C_l_h2 = np.loadtxt('data/CMB_spectrum_h2.dat', unpack=True)

# Reading in power spectra observed data
l_real, CMB_real, ndC, pdC = np.loadtxt('data/COM_PowerSpect_CMB-TT-full_R3.01.txt', unpack=True)

fig_CMB_default = plt.figure()
plt.errorbar(l_real, CMB_real, yerr=[ndC,pdC], color='lightblue', label='Observed', zorder=1)
plt.plot(l, C_l*5775/max(C_l), label='Default')
plt.ylabel(r'$l(l+1) C_l / 2 \pi$ [$\mu$K$^2$]')
plt.xlabel(r'$l$')
plt.legend(loc='best')
plt.xlim(l[0], l[-1])
plt.grid()
plt.show()


fig_CMB_n = plt.figure()
plt.errorbar(l_real, CMB_real, yerr=[ndC,pdC], color='lightblue', label='Observed', zorder=1)
plt.plot(l, C_l*5775/max(C_l), label='Default')
plt.plot(l_n1, C_l_n1*5775/max(C_l_n1), label=r'$n_s = 1.1$')
plt.plot(l_n2, C_l_n2*5775/max(C_l_n2), label=r'$n_s = 0.8$')
plt.ylabel(r'$l(l+1) C_l / 2 \pi$ [$\mu$K$^2$]')
plt.xlabel(r'$l$')
plt.legend(loc='best')
plt.xlim(l[0], l[-1])
plt.grid()
plt.show()

fig_CMB_b = plt.figure()
plt.errorbar(l_real, CMB_real, yerr=[ndC,pdC], color='lightblue', label='Observed', zorder=1)
plt.plot(l, C_l*5775/max(C_l), label='Default')
plt.plot(l_b1, C_l_b1*5775/max(C_l_b1), label=r'$\Omega_b = 0.066$')
plt.plot(l_b2, C_l_b2*5775/max(C_l_b2), label=r'$\Omega_b = 0.026$')
plt.ylabel(r'$l(l+1) C_l / 2 \pi$ [$\mu$K$^2$]')
plt.xlabel(r'$l$')
plt.legend(loc='best')
plt.xlim(l[0], l[-1])
plt.grid()
plt.show()

fig_CMB_m = plt.figure()
plt.errorbar(l_real, CMB_real, yerr=[ndC,pdC], color='lightblue', label='Observed', zorder=1)
plt.plot(l, C_l*5775/max(C_l), label='Default')
plt.plot(l_m1, C_l_m1*5775/max(C_l_m1), label=r'$\Omega_m = 0.324$')
plt.plot(l_m2, C_l_m2*5775/max(C_l_m2), label=r'$\Omega_m = 0.124$')
plt.ylabel(r'$l(l+1) C_l / 2 \pi$ [$\mu$K$^2$]')
plt.xlabel(r'$l$')
plt.legend(loc='best')
plt.xlim(l[0], l[-1])
plt.grid()
plt.show()

fig_CMB_h = plt.figure()
plt.errorbar(l_real, CMB_real, yerr=[ndC,pdC], color='lightblue', label='Observed', zorder=1)
plt.plot(l, C_l*5775/max(C_l), label='Default')
plt.plot(l_h1, C_l_h1*5775/max(C_l_h1), label=r'$h = 0.8$')
plt.plot(l_h2, C_l_h2*5775/max(C_l_h2), label=r'$h = 0.6$')
plt.ylabel(r'$l(l+1) C_l / 2 \pi$ [$\mu$K$^2$]')
plt.xlabel(r'$l$')
plt.legend(loc='best')
plt.xlim(l[0], l[-1])
plt.grid()
plt.show()


if save:
    fig_CMB_default.savefig("figures/CMB_default.pdf", bbox_inches='tight')
    fig_CMB_n.savefig("figures/CMB_n.pdf", bbox_inches='tight')
    fig_CMB_b.savefig("figures/CMB_b.pdf", bbox_inches='tight')
    fig_CMB_m.savefig("figures/CMB_m.pdf", bbox_inches='tight')
    fig_CMB_h.savefig("figures/CMB_h.pdf", bbox_inches='tight')
