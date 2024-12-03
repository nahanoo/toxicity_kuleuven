import pandas as pd
from os.path import join
import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt
from scipy.integrate import odeint
import scipy.optimize as optimize
import plotly.graph_objects as go
from style import *

paramters = ''
fig_total = go.Figure()

f = join("~", "curves", "metadata", "pooled_df_joint_metadata.csv")
pooled_meta = pd.read_csv(f)
f = join("~", "curves", "export", "240903_fran", "measurement_data.csv")
raw = pd.read_csv(f)
strain = 'Salmonella Typhimurium mcherry:pGDPI CTX-M-15'
conc = 0.25
filter = (pooled_meta["project"]
          == "240903_fran") & (pooled_meta["Experiment description"].isin([
              "Repeat 1 of TSB gradient experiment with all strains.",
              "Repeat 2 of TSB gradient experiment with all strains.",
              "Repeat 3 of TSB gradient experiment with all strains.",
          ]) & (pooled_meta['cs_conc'] == conc) &
                               (pooled_meta['species'] == strain))

meta = pooled_meta[filter]
xs = [
    raw[raw["linegroup"] == lg]["time"].to_numpy() for lg in meta['linegroup']
]
ys = [
    raw[raw["linegroup"] == lg]["measurement"].to_numpy()
    for lg in meta['linegroup']
]
show_legend = True
for x, y in zip(xs, ys):
    fig_total.add_trace(
        go.Scatter(x=x,
                   y=y,
                   name='Salmonella Typhimurium mcherry:pGDPI CTX-M-15',
                   showlegend=show_legend,
                   line=dict(color=colors['st'])), )
    show_legend = False

x = np.average(xs, axis=0)
y = np.average(ys, axis=0)
x_window = [0, 10.2]
for i, q in enumerate(x):
    if x_window[0] < q:
        j = i
        break
for i, q in enumerate(x):
    if x_window[1] < q:
        k = i
        break
x, y = x[j:k], y[j:k]
q = 0.52 / conc


def max_growth_rate(x_window, x, y):
    for i, q in enumerate(x):
        if x_window[0] < q:
            j = i
            break
    for i, q in enumerate(x):
        if x_window[1] < q:
            k = i
            break
    x, y = x[j:k], y[j:k]
    return linregress(x, np.log(y))[0]


def monod(y, t, v, Km, q):
    n, c = y
    dndt = v * c / (Km + c) * n
    dcdt = -v * c / (Km + c) * n / q
    return np.array([dndt, dcdt])


def simulate_monod(Km, v, t, q, n, c0, n0):
    y = odeint(monod, [n0, c0], t, args=(v, Km[0], q))
    return np.sum((n - y[:, 0])**2)


v = max_growth_rate([2.5, 8], x, y)
c0 = 0.25
n0 = 0.01
Km = optimize.minimize(
    simulate_monod,
    [0.0001],
    args=(
        v,
        x,
        q,
        y,
        c0,
        n0,
    ),
    bounds=((0, np.inf), ),
).x[0]
sim = odeint(monod, [n0, c0], x, args=(v, Km, q))
fig_fit_st = go.Figure()
show_legend = True
for x, y in zip(xs, ys):
    fig_fit_st.add_trace(
        go.Scatter(x=x,
                   y=y,
                   name='Salmonella Typhimurium mcherry:pGDPI CTX-M-15',
                   showlegend=show_legend,
                   line=dict(color=colors['st'])))
    show_legend = False
fig_fit_st.add_trace(
    go.Scatter(x=x,
               y=sim[:, 0],
               name='Model',
               line=dict(color=colors['st'], dash='dash')))
fig_fit_st = style_plot(fig_fit_st)
fig_fit_st.update_layout(width=width * 2, height=height)
fig_fit_st.update_xaxes(title='Time [h]'), fig_fit_st.update_yaxes(
    title='Abundance [OD]')
fig_fit_st.write_image('../figures/fig1_st_fit.svg')
paramters += 'ST v:' + str(v) + ',ST Km:' + str(Km) + '\n'

f = join("~", "curves", "metadata", "pooled_df_joint_metadata.csv")
pooled_meta = pd.read_csv(f)
f = join("~", "curves", "export", "240903_fran", "measurement_data.csv")
raw = pd.read_csv(f)
strain = 'Escherichia coli  sfGFP:pGDPI CAT'
conc = 0.25
filter = (pooled_meta["project"]
          == "240903_fran") & (pooled_meta["Experiment description"].isin([
              "Repeat 1 of TSB gradient experiment with all strains.",
              "Repeat 2 of TSB gradient experiment with all strains.",
              "Repeat 3 of TSB gradient experiment with all strains.",
          ]) & (pooled_meta['cs_conc'] == conc) &
                               (pooled_meta['species'] == strain))

meta = pooled_meta[filter]
xs = [
    raw[raw["linegroup"] == lg]["time"].to_numpy() for lg in meta['linegroup']
]
ys = [
    raw[raw["linegroup"] == lg]["measurement"].to_numpy()
    for lg in meta['linegroup']
]
show_legend = True
for x, y in zip(xs, ys):
    fig_total.add_trace(
        go.Scatter(x=x,
                   y=y,
                   name='Escherichia coli  sfGFP:pGDPI CAT',
                   showlegend=show_legend,
                   line=dict(color=colors['ecoli'])), )
    show_legend = False

x = np.average(xs, axis=0)
y = np.average(ys, axis=0)
x_window = [0, 15]
for i, q in enumerate(x):
    if x_window[0] < q:
        j = i
        break
for i, q in enumerate(x):
    if x_window[1] < q:
        k = i
        break
x, y = x[j:k], y[j:k]
q = 0.52 / conc

v = max_growth_rate([2.5, 8], x, y)
c0 = 0.25
n0 = 0.01
Km = optimize.minimize(
    simulate_monod,
    [0.0001],
    args=(
        v,
        x,
        q,
        y,
        c0,
        n0,
    ),
    bounds=((0, np.inf), ),
).x[0]
sim = odeint(monod, [n0, c0], x, args=(v, Km, q))
fig_fit_e_coli = go.Figure()
show_legend = True
for x, y in zip(xs, ys):
    fig_fit_e_coli.add_trace(
        go.Scatter(x=x,
                   y=y,
                   name='Escherichia coli  sfGFP:pGDPI CAT',
                   showlegend=show_legend,
                   line=dict(color=colors['ecoli'])))
    show_legend = False
fig_fit_e_coli.add_trace(
    go.Scatter(x=x,
               y=sim[:, 0],
               name='Model',
               line=dict(color=colors['ecoli'], dash='dash')))
fig_fit_e_coli = style_plot(fig_fit_e_coli)
fig_fit_e_coli.update_layout(width=width * 2, height=height)
fig_fit_e_coli.update_xaxes(title='Time [h]'), fig_fit_e_coli.update_yaxes(
    title='Abundance [OD]')
fig_fit_e_coli.write_image('../figures/fig1_e_coli_fit.svg')

fig_total = style_plot(fig_total)
fig_total.update_layout(width=width * 2, height=height)
fig_total.update_xaxes(title='Time [h]'), fig_total.update_yaxes(
    title='Abundance [OD]')
fig_total.write_image('../figures/fig1_e_coli_st.svg')
paramters += 'Ecoli v:' + str(v) + ',Ecoli Km:' + str(Km) + '\n'
with open('parameters.txt', 'w') as handle:
    handle.write(paramters)
