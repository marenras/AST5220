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

z, X_e = np.loadtxt('data/X_e.dat', unpack=True)
x, tau, dtau, ddtau = np.loadtxt('data/tau.dat', unpack=True)
x, g, dg, ddg = np.loadtxt('data/g.dat', unpack=True)

plt.figure()
plt.plot(z, X_e)
plt.title("Fractional electron density")
plt.xlabel(r'Redshift $z$')
plt.ylabel(r'$X_e$')
plt.xlim(1800, 0)
plt.grid()
if save:
    tikz_save("figures/X_e.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")

plt.figure()
plt.semilogy(x, tau, label=r'$\tau$')
plt.semilogy(x[:-5], abs(dtau[:-5]), label=r'$|\dot{tau}|$')
plt.semilogy(x[:-8], ddtau[:-8], label=r'$\ddot{tau}$')
plt.title('Optical depth')
plt.xlabel(r'$x = \ln a$')
plt.ylabel(r'$\tau$')
plt.legend()
plt.grid()
if save:
    tikz_save("figures/tau.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")


plt.figure()
plt.plot(x, g, label=r'$\tilde{g}$')
plt.plot(x, dg/10, label=r'$\dot{\tilde{g}}/10$')
plt.plot(x, ddg/300, label=r'$\ddot{\tilde{g}}/300$')
plt.title('Visibility function')
plt.xlabel(r'$x = \ln a$')
plt.ylabel(r'$\tilde{g}$')
plt.legend()
plt.xlim(-8, -6)
plt.grid()
if save:
    tikz_save("figures/g.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")


plt.show()
