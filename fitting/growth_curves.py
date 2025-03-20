import pandas as pd
from os.path import join
import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt
from scipy.integrate import odeint
import scipy.optimize as optimize
import plotly.graph_objects as go
from style import *
from plotly.subplots import make_subplots

paramters = ""
fig_total = go.Figure()

f = join("~", "curves", "metadata", "pooled_df_joint_metadata.csv")
pooled_meta = pd.read_csv(f)
f = join("~", "curves", "export", "240903_fran", "measurement_data.csv")
raw = pd.read_csv(f)
strain = "Salmonella Typhimurium mcherry:pGDPI CTX-M-15"
conc = 1
filter = (pooled_meta["project"] == "240903_fran") & (
    pooled_meta["Experiment description"].isin(
        [
            "Repeat 1 of TSB gradient experiment with all strains.",
            "Repeat 2 of TSB gradient experiment with all strains.",
            "Repeat 3 of TSB gradient experiment with all strains.",
        ]
    )
    & (pooled_meta["cs_conc"] == conc)
    & (pooled_meta["species"] == strain)
)

meta = pooled_meta[filter]
xs = [raw[raw["linegroup"] == lg]["time"].to_numpy() for lg in meta["linegroup"]]
ys = [raw[raw["linegroup"] == lg]["measurement"].to_numpy() for lg in meta["linegroup"]]
show_legend = True
for x, y in zip(xs, ys):
    fig_total.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name="Salmonella Typhimurium <br>mcherry:pGDPI CTX-M-15",
            showlegend=show_legend,
            line=dict(color=colors["st"]),
        ),
    )
    show_legend = False


def logistic(y, t, v, K):
    N = y[0]
    dN = v * (1 - N / K) * N
    return dN


ts = np.linspace(0, x[-1], 500)
sim = odeint(logistic, [0.001], ts, args=(0.7, 1.08))
fig_fit_st = go.Figure()
show_legend = True
for x, y in zip(xs, ys):
    fig_fit_st.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name="Salmonella Typhimurium <br>mcherry:pGDPI CTX-M-15",
            showlegend=show_legend,
            line=dict(color=colors["st"]),
        )
    )
    show_legend = False
fig_fit_st.add_trace(
    go.Scatter(x=ts, y=sim[:, 0], name="Model", line=dict(color="black", dash="dash"))
)
fig_fit_st = style_plot(fig_fit_st)
fig_fit_st.update_layout(width=width, height=height * 0.8)
fig_fit_st.update_xaxes(title="Time [h]"), fig_fit_st.update_yaxes(
    title="Abundance [OD]"
)
fig_fit_st.write_image("../plots/experiments/st_fit.svg")

f = join("~", "curves", "metadata", "pooled_df_joint_metadata.csv")
pooled_meta = pd.read_csv(f)
f = join("~", "curves", "export", "240903_fran", "measurement_data.csv")
raw = pd.read_csv(f)
strain = "Escherichia coli  sfGFP:pGDPI CAT"
conc = 1
filter = (pooled_meta["project"] == "240903_fran") & (
    pooled_meta["Experiment description"].isin(
        [
            "Repeat 1 of TSB gradient experiment with all strains.",
            "Repeat 2 of TSB gradient experiment with all strains.",
            "Repeat 3 of TSB gradient experiment with all strains.",
        ]
    )
    & (pooled_meta["cs_conc"] == conc)
    & (pooled_meta["species"] == strain)
)

meta = pooled_meta[filter]
xs = [raw[raw["linegroup"] == lg]["time"].to_numpy() for lg in meta["linegroup"]]
ys = [raw[raw["linegroup"] == lg]["measurement"].to_numpy() for lg in meta["linegroup"]]
show_legend = True
for x, y in zip(xs, ys):
    fig_total.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name="Escherichia coli  <br>sfGFP:pGDPI CAT",
            showlegend=show_legend,
            line=dict(color=colors["ecoli"]),
        ),
    )
    show_legend = False

ts = np.linspace(0, x[-1], 500)
sim = odeint(logistic, [0.001], ts, args=(0.9, 1.08))
fig_fit_e_coli = go.Figure()
show_legend = True
for x, y in zip(xs, ys):
    fig_fit_e_coli.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name="Escherichia coli  <br>sfGFP:pGDPI CAT",
            showlegend=show_legend,
            line=dict(color=colors["ecoli"]),
        )
    )
    show_legend = False
fig_fit_e_coli.add_trace(
    go.Scatter(x=ts, y=sim[:, 0], name="Model", line=dict(color="black", dash="dash"))
)
fig_fit_e_coli = style_plot(fig_fit_e_coli)
fig_fit_e_coli.update_layout(width=width * 2, height=height)
fig_fit_e_coli.update_xaxes(title="Time [h]"), fig_fit_e_coli.update_yaxes(
    title="Abundance [OD]"
)
fig_fit_e_coli.write_image("../plots/experiments/e_coli_fit.svg")
fig_total.update_layout(title="Mono-culture growth curves 100% CFA")
fig_total = style_plot(fig_total, marker_size=1.5, line_thickness=2)
fig_total.update_layout(width=width, height=height * 0.7)
fig_total.update_xaxes(title="Time [h]"), fig_total.update_yaxes(title="Abundance [OD]")
fig_total.write_image("../plots/experiments/e_coli_st.svg")
