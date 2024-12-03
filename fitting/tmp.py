import pandas as pd
from os.path import join
import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt
from scipy.integrate import odeint
import scipy.optimize as optimize
import plotly.graph_objects as go
from style import *


def model(y, t):
    N0 = y[0]
    dN0 = 0.8 * N0 * (1 - N0 / 1)
    return dN0


J = [0.8 * (1 - N / 1) for N in np.linspace(0.01, 1, 20)]
Cm = np.linspace(0, 16, 10)
J_grid, Cm_grid = np.meshgrid(np.array(J), np.array(Cm))
JD = J_grid * 1 / (1 + Cm_grid / 0.00030537131107306284)
plt.contourf(J_grid, Cm_grid, JD, levels=100)
plt.colorbar(label="μ [1/h]")
plt.show()

J = [0.8 * (1 - N / 1) for N in np.linspace(0.01, 1, 20)]
Cf = np.linspace(0, 0.25, 10)
J_grid, Cf_grid = np.meshgrid(np.array(J), np.array(Cf))
JD = J_grid - Cf_grid / (Cf_grid + 0.00030537131107306284)
plt.contourf(J_grid, Cf_grid, JD, levels=100)
plt.colorbar(label="μ [1/h]")
plt.show()
