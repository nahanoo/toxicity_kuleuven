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


def logistic(y, t, v, K):
    N = y[0]
    dN = v * (1 - N / K) * N
    return dN


fig_st_e = go.Figure()


def st():
    fig = go.Figure()
    f = join("~", "curves", "metadata", "pooled_df_joint_metadata.csv")
    pooled_meta = pd.read_csv(f)
    f = join("~", "curves", "export", "240903_fran", "measurement_data.csv")
    raw = pd.read_csv(f)
    strain = "Salmonella Typhimurium mcherry:pGDPI CTX-M-15"
    conc = 0.25
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
    ys = [
        raw[raw["linegroup"] == lg]["measurement"].to_numpy()
        for lg in meta["linegroup"]
    ]
    sim = odeint(logistic, [0.001], xs[0], args=(1, 0.55))[:, 0]
    for i, (x, y) in enumerate(zip(xs, ys)):
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name="<i>Salmonella Typhimurium</i> <br>mcherry:pGDPI CTX-M-15",
                line=dict(color=colors["st"]),
                showlegend=i == 0,
            )
        )
    fig.add_trace(
        go.Scatter(
            x=xs[0],
            y=sim,
            name="Logistic fit",
            line=dict(color=colors["st"], dash="dash"),
            showlegend=True,
        )
    )
    fig.update_layout(
        xaxis=dict(title="Time [h]", range=[0, 24], dtick=8),
        yaxis=dict(title="OD", range=[0, 0.6], dtick=0.0),
        width=width,
        height=height,
    )
    fig = style_plot(
        fig,
        line_thickness=1,
        buttom_margin=0,
        left_margin=40,
        top_margin=0,
        right_margin=0,
    )
    fig.write_image("../plots/experiments/st_fit_25.svg")
    fig_st_e.add_trace(fig["data"][0])
    fig_st_e.add_trace(fig["data"][1])
    fig_st_e.add_trace(fig["data"][2])


def ecoli():
    fig = go.Figure()
    f = join("~", "curves", "metadata", "pooled_df_joint_metadata.csv")
    pooled_meta = pd.read_csv(f)
    f = join("~", "curves", "export", "240903_fran", "measurement_data.csv")
    raw = pd.read_csv(f)
    strain = "Escherichia coli  sfGFP:pGDPI CAT"
    conc = 0.25
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
    ys = [
        raw[raw["linegroup"] == lg]["measurement"].to_numpy()
        for lg in meta["linegroup"]
    ]
    sim = odeint(logistic, [0.0019], xs[0], args=(0.8, 0.55))[:, 0]
    for i, (x, y) in enumerate(zip(xs, ys)):
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name="<i>Escherichia coli</i><br>sfGFP:pGDPI CAT",
                line=dict(color=colors["ecoli"]),
                showlegend=i == 0,
            )
        )
    fig.add_trace(
        go.Scatter(
            x=xs[0],
            y=sim,
            name="Logistic fit",
            line=dict(color=colors["ecoli"], dash="dash"),
            showlegend=True,
        )
    )
    fig.update_layout(
        xaxis=dict(title="Time [h]", range=[0, 24], dtick=8),
        yaxis=dict(title="OD", range=[0, 0.6], dtick=0.0),
        width=width,
        height=height,
    )
    fig = style_plot(
        fig,
        line_thickness=1,
        buttom_margin=0,
        left_margin=40,
        top_margin=0,
        right_margin=0,
    )
    fig.write_image("../plots/experiments/e_fit_25.svg")
    fig_st_e.add_trace(fig["data"][0])
    fig_st_e.add_trace(fig["data"][1])
    fig_st_e.add_trace(fig["data"][2])


fig_st_e.update_layout(
    xaxis=dict(title="Time [h]", range=[0, 24], dtick=8),
    yaxis=dict(title="OD", range=[0, 0.6], dtick=0.0),
    width=width,
    height=height,
)
fig_st_e = style_plot(
    fig_st_e,
    line_thickness=1,
    buttom_margin=0,
    left_margin=40,
    top_margin=0,
    right_margin=0,
)
st()
ecoli()
fig_st_e.write_image("../plots/experiments/e_st_fit_25.svg")
