from scipy.optimize import root_scalar
import numpy as np

eps = np.finfo(float).eps


def Lx(x, mu):
    return -(1 - mu) / ((x + mu) * abs(x + mu)) - mu / ((x - (1 - mu)) * abs(x - (1 - mu))) + x


def L1(mu):
    return root_scalar(Lx, args=(mu), bracket=[-mu + eps, (1 - mu) - eps]).root


def L2(mu):
    return root_scalar(Lx, args=(mu), bracket=[(1-mu) + eps, 2]).root


def L3(mu):
    return root_scalar(Lx, args=(mu), bracket=[-2, -mu - eps]).root


def fyy(mu, Lx):
    return 1 - (1 - mu) / abs(Lx + mu)**3 - mu / abs(Lx - (1 - mu))**3


fyy1 = []
fyy2 = []
fyy3 = []
mu_range = np.linspace(eps, eps*1000, 1000)

for mu in mu_range:
    fyy1.append(fyy(mu, L1(mu)))
    fyy2.append(fyy(mu, L2(mu)))
    fyy3.append(fyy(mu, L3(mu)))

import matplotlib.pyplot as plt

plt.plot(mu_range, fyy1, 'b')
plt.plot(mu_range, fyy2, 'r')
plt.plot(mu_range, fyy3, 'g')
plt.grid()