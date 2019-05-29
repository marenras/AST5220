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


phi = []
psi = []
delta = []
delta_b = []
v = []
v_b = []
theta0 = []
theta1 = []

k_values = np.loadtxt('data/k_values.dat', unpack=True)

for k in range(1,7):
    x, phi_, psi_, delta_, delta_b_,  v_, v_b_, theta0_, theta1_ = np.loadtxt('data/perturbations_%d.dat' %(k), unpack=True)
    phi.append(phi_)
    psi.append(psi_)
    delta.append(delta_)
    delta_b.append(delta_b_)
    v.append(v_)
    v_b.append(v_b_)
    theta0.append(theta0_)
    theta1.append(theta1_)

fig_phi = plt.figure()
for i in range(6):
    plt.plot(x, phi[i], label=r'$kc/H_0$ = %.2f' %(k_values[i]))
plt.ylabel(r'$\Phi$')
plt.xlabel(r'$x = \ln a$')
plt.legend(loc='best')
plt.grid()
plt.show()

fig_psi = plt.figure()
for i in range(6):
    plt.plot(x, psi[i], label=r'$kc/H_0$ = %.2f' %(k_values[i]))
plt.ylabel(r'$\Psi$')
plt.xlabel(r'$x = \ln a$')
plt.legend(loc='best')
plt.grid()
plt.show()

fig_delta = plt.figure()
for i in range(6):
    plt.semilogy(x, delta[i], label=r'$kc/H_0$ = %.2f' %(k_values[i]))
plt.ylabel(r'$\delta$')
plt.xlabel(r'$x = \ln a$')
plt.legend(loc='best')
plt.grid()
plt.show()

fig_delta_b = plt.figure()
for i in range(6):
    plt.semilogy(x, delta_b[i], label=r'$kc/H_0$ = %.2f' %(k_values[i]))
plt.ylabel(r'$\delta_b$')
plt.xlabel(r'$x = \ln a$')
plt.legend(loc='best')
plt.grid()
plt.show()

fig_v = plt.figure()
for i in range(6):
    plt.plot(x, v[i], label=r'$kc/H_0$ = %.2f' %(k_values[i]))
plt.ylabel(r'$v$')
plt.xlabel(r'$x = \ln a$')
plt.legend(loc='best')
plt.grid()
plt.show()

fig_v_b = plt.figure()
for i in range(6):
    plt.plot(x, v_b[i], label=r'$kc/H_0$ = %.2f' %(k_values[i]))
plt.ylabel(r'$v_b$')
plt.xlabel(r'$x = \ln a$')
plt.legend(loc='best')
plt.grid()
plt.show()

fig_theta0 = plt.figure()
for i in range(6):
    plt.plot(x, theta0[i], label=r'$kc/H_0$ = %.2f' %(k_values[i]))
plt.ylabel(r'$\Theta_0$')
plt.xlabel(r'$x = \ln a$')
plt.legend(loc='best')
plt.grid()
plt.show()

if save:
    fig_phi.savefig("figures/phi.pdf", bbox_inches='tight')
    fig_psi.savefig('figures/psi.pdf', bbox_inches='tight')
    fig_delta.savefig("figures/delta.pdf", bbox_inches='tight')
    fig_delta_b.savefig("figures/delta_b.pdf", bbox_inches='tight')
    fig_v.savefig("figures/v.pdf", bbox_inches='tight')
    fig_v_b.savefig("figures/v_b.pdf", bbox_inches='tight')
    fig_theta0.savefig("figures/theta0.pdf", bbox_inches='tight')
