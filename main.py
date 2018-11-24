from PricingLib import BSprice
from PricingLib import MCsimulations
from PricingLib import MCLEuler
from mpl_toolkits . mplot3d import Axes3D
from greeks import GreekBS, GreekFD, GreekMCL, GreekMCLlikelihoodRatio
from random import randint
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np


# Option parameters
sigma = 0.12  # Flat volatility
S = 100.0  # Stock price
K = 105.0  # Fixed strike
epsilon = 0.4  # The % on the left / right of Strike .
# Asset prices are centered around Spot (" ATM Spot ")
T1 = 30/260  # Shortest expiry in days
T2 = 720/260  # Longest expiry in days
r = 0.05  # Continuous risk free rate
d = 0.02  # Continuous div rate
e = 0.01  # small value used to calculate price sensitivities
N = 10000  # number of simulations
#
#A = [BSprice(1, 0.05, 0.00, 90, i, 0.20).digitcall() for i in range(50, 150, 5)]
#B = [MCsimulations(1, 0.05, 0.00, 90, i, 0.20, 10000).digitcall() for i in range(50, 150, 5)]
#C = [0.01*BSprice(1, 0.05, 0.00, 90, i, 0.20).callprice() for i in range(50, 150, 5)]

#A = [0.01*MCsimulations(1, 0.05, 0.00, 100, 100, i*0.01, 10000).callprice() for i in range(1, 50, 1)]
# A = [MCsimulations(1, 0.05, 0.00, 90, 100, i*0.01, 10000).digitcall() for i in range(1, 50, 1)]
# B = [MCsimulations(1, 0.05, 0.00, 110, 100, i*0.01, 10000).digitcall() for i in range(1, 50, 1)]
# C = [MCsimulations(1, 0.05, 0.00, 100, 100, i*0.01, 10000).digitcall() for i in range(1, 50, 1)]
#C = [MCLEuler(1, 0.05, 0.00, 100, i, 0.25, 1000, 10).digitcall() for i in range(1, 200, 5)]

A = [GreekBS(T2, r, d, K, i, sigma).delta() for i in range(50, 150, 1)]
B = [GreekFD(T2, r, d, K, i, sigma, e).delta() for i in range(50, 150, 1)]
C = [GreekMCL(T2, r, d, K, i, sigma, N, e, 1).delta() for i in range(50, 150, 1)]
D = [GreekMCLlikelihoodRatio(T2, r, d, K, i, sigma, N).delta() for i in range(50, 150, 1)]



# C = [GreekMCL(T2, r, d, K, i, sigma, N, e, randint(1, 10000)).gamma() for i in range(50, 150, 1)]

plt.plot(A, label='BS delta')
plt.plot(B, label='FD delta')
plt.plot(C, label='MCL delta')
plt.plot(D, label='LR delta')
plt.legend()
plt.show()

# # Grid definition
# dx, dy = 40, 40  # Steps throughout asset price and expiries axis
#
# # xx: Asset price axis , yy: expiry axis , zz: greek axis
# xx, yy = np.meshgrid(np.linspace(K*(1 - epsilon), (1 + epsilon) * K, dx), np.linspace(T1, T2, dy))
# zz = np.array([GreekBS(y, r, d, K, x, sigma).rho() for x, y in zip(np.ravel(xx), np.ravel(yy))])
# zz = zz. reshape(xx.shape)
#
# # Plot greek surface
# print(" Plotting surface ... ")
# fig = plt.figure()
# fig.suptitle('Call delta', fontsize=20)
# ax = fig.gca(projection='3d')
# surf = ax. plot_surface(xx, yy, zz, rstride=1, cstride=1, alpha=0.75, cmap=cm.RdYlBu)
# ax.set_xlabel('Asset price')
# ax.set_ylabel('Expiry')
# ax.set_zlabel('Delta')
# # Plot 3D contour
# zzlevels = np.linspace(zz.min(), zz. max(), num=8, endpoint=True)
# xxlevels = np.linspace(xx.min(), xx. max(), num=8, endpoint=True)
# yylevels = np.linspace(yy.min(), yy. max(), num=8, endpoint=True)
#
# ax.set_xlim(xx.min(), xx.max())
# ax.set_ylim(yy. min(), yy.max())
# ax.set_zlim(zz.min(), zz.max())
# plt.show()