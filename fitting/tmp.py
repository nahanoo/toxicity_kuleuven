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


J = np.linspace(0, 0.9, 100)
Cf = np.linspace(0, 0.25, 100)
J_grid, Cf_grid = np.meshgrid(np.array(J), np.array(Cf))
JD = J_grid - 1.5 * Cf_grid / (Cf_grid + 0.029)
plt.contourf(J_grid, Cf_grid, JD, levels=20, cmap="coolwarm")
plt.colorbar(label="Netto growth rate [1/h]")
plt.xlabel("Per capita growth rate µ [1/h]")
plt.ylabel("Cefotaxime concentration μg/mL")
plt.title("J = -1.5 1/h")
plt.tight_layout()
plt.savefig("../figures/fig5_cefotaxime_isoclines_1_5.svg")
plt.close()

J = np.linspace(0, 0.7, 100)
Cf = np.linspace(0, 16, 100)
J_grid, Cf_grid = np.meshgrid(np.array(J), np.array(Cf))
JD = J_grid * 1 / (1 + Cf_grid / 0.32)
plt.contourf(J_grid, Cf_grid, JD, levels=20, cmap="coolwarm")
plt.colorbar(label="Netto growth rate [1/h]")
plt.xlabel("Per capita growth rate µ [1/h]")
plt.ylabel("Chloramphenicol concentration μg/mL")
plt.tight_layout()
plt.savefig("../figures/fig5_chloramphenicol.svg")
plt.close()
