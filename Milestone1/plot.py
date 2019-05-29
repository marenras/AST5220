import numpy as np
import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save
import sys as sys


save = False
if len(sys.argv) > 1:
    if sys.argv[1] in ['save', 'Save']:
        save = True
    else:
	    print("Please write 'save' in command line to save the plots as tex-files.")
	    sys.exit()

x, H_x = np.loadtxt('data/hubble_x.dat', unpack=True)
z, H_z = np.loadtxt('data/hubble_z.dat', unpack=True)
x_eta, eta = np.loadtxt('data/eta.dat', unpack=True)
x_eta, omega_m, omega_b, omega_r, omega_lambda = np.loadtxt('data/omegas.dat', unpack=True)
x_t, splined_eta = np.loadtxt('data/splined_eta.dat', unpack=True)


# Converting from km to Mpc:
unit = 1e-3 * 3.08568025e22

plt.figure()
plt.semilogy(x, H_x * unit)
plt.title("Hubble parameter with respect to x")
plt.xlabel(r'$x = \ln a$')
plt.ylabel('H(x) [km s$^{-1}$Mpc$^{-1}$]')
plt.grid()
if save:
    tikz_save("figures/hubble_x.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")

plt.figure()
plt.loglog(z, H_z * unit)
plt.title("Hubble parameter with respect to z")
plt.gca().invert_xaxis()
plt.xlabel(r'$z$')
plt.ylabel('H(z) [km s$^{-1}$Mpc$^{-1}$]')
plt.grid()
if save:
    tikz_save("figures/hubble_z.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")

plt.figure()
plt.semilogy(x_eta, eta, label=r'$\eta$')
plt.semilogy(x_t, splined_eta, '--', color='red', label=r'Splined $\eta$')
plt.legend()
plt.title("Conformal time")
plt.xlabel(r'$x = \ln a$')
plt.ylabel(r'$\eta$')
plt.grid()
if save:
    tikz_save("figures/eta.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")


plt.figure()
plt.plot(x_eta, omega_m, label=r'$\Omega_m$')
plt.plot(x_eta, omega_b, label=r'$\Omega_b$')
plt.plot(x_eta, omega_r, label=r'$\Omega_r$')
plt.plot(x_eta, omega_lambda, label=r'$\Omega_\Lambda$')
plt.legend()
plt.title('Density parameters')
plt.xlabel(r'$x = \ln a$')
plt.ylabel(r'$\Omega$')
plt.grid()
if save:
    tikz_save("figures/omegas.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")
plt.show()
