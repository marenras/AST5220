import numpy as np
import matplotlib.pyplot as plt



j, x = np.loadtxt('j_l1.dat', unpack=True)
j2, x = np.loadtxt('j_l2.dat', unpack=True)
S, x = np.loadtxt('S.dat', unpack=True)

plt.figure()
plt.plot(x,j)
plt.ylabel('j')

plt.figure()
plt.plot(x,j2)
plt.ylabel('j2')

plt.figure()
plt.plot(x,S)
plt.ylabel('S')

plt.show()
