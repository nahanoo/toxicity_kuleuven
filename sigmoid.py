import numpy as np
from scipy.integrate import odeint
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import pandas as pd

r = 0.05
K = 50
M = 10
Q = 0.025
N0 = 0.075

def cr(y,t):
    N,R = y
    dN = r * R / (R + K) * N
    dR = -r * R / (R + K) * N / Q
    return dN, dR


def four_param_logistic(x, A, Km, B, C):
    return A / (1 + (x / Km)**B) + C


df = pd.read_csv('/home/eric/ChiBioFlow/data/at_oa/growth_assays/carbon_sources.csv')
df = df[(df['well'] == 'G10') & (df['strain'] == 'ct')]

"""xs = np.linspace(0,72,200)
Y = odeint(cr,[N0,M],xs)
N, R = Y[:,0], Y[:,1]"""
xs = np.linspace(0.5,75.5,len(df)) [10:]
N = np.array(df['OD'].to_list())[10:] /2

params,_ = curve_fit(four_param_logistic,xs,N)
A, Km, B, C = params

fitted_population_data = four_param_logistic(xs, A, Km, B, C)
plt.figure(figsize=(10, 5))
plt.plot(xs, N, 'bo', label='Original Data')
plt.plot(xs, fitted_population_data, 'r-', label='Fitted Curve')
plt.xlabel('Time')
plt.ylabel('Bacterial Population')
plt.title('4PL Model Fitting to Bacterial Growth Data')
plt.legend()
plt.grid(True)
plt.show()