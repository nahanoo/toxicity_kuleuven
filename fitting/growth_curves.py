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
    dN = v * N
    return dN


fig_st_e = go.Figure()
f = join("~", "toxicity_kuleuven", "data", "nut_gradient", "measurement_data.csv")
raw = pd.read_csv(f)


def filter_time(xs, ys, x0=2, x1=20):
    xs_f = []
    ys_f = []
    for x, y in zip(xs, ys):
        mask = (x >= x0) & (x <= x1)
        xs_f.append(x[mask])
        ys_f.append(y[mask])
    return xs_f, ys_f


def get_lgs(strain, conc):
    f = join(
        "~", "toxicity_kuleuven", "data", "nut_gradient", "pooled_df_joint_metadata.csv"
    )
    pooled_meta = pd.read_csv(f)
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
    return meta["linegroup"].to_list()


def st_1():
    conc = 1
    lgs = get_lgs("Salmonella Typhimurium mcherry:pGDPI CTX-M-15", conc)
    xs = [raw[raw["linegroup"] == lg]["time"].to_numpy() for lg in lgs]
    ys = [raw[raw["linegroup"] == lg]["measurement"].to_numpy() for lg in lgs]
    y_diff = [np.array([min(y) + 1e-6]) for y in ys]
    ys = [y - yd for y, yd in zip(ys, y_diff)]
    xs, ys = filter_time(
        xs,
        ys,
    )
    r = 1.7
    sim = odeint(logistic, [np.average([y[0] for y in ys])], xs[0], args=(r, 1.1))[:, 0]
    y = np.mean(ys, axis=0)
    x = np.mean(xs, axis=0)  # if all xs are identical, you could just use xs[0]
    y_lower = np.percentile(ys, 16, axis=0)
    y_upper = np.percentile(ys, 84, axis=0)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=np.concatenate([x, x[::-1]]),
            y=np.concatenate([y_upper, y_lower[::-1]]),
            fill="toself",
            fillcolor=colors["st_error"],
            line=dict(color=colors["st_error"]),
            hoverinfo="skip",
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name="<i>Salmonella Typhimurium</i> <br>mcherry:pGDPI CTX-M-15",
            line=dict(color=colors["st"]),
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=[s for s in sim if s <= max(y)],
            name="Logistic fit",
            mode="lines",
            line=dict(color=colors["st"], dash="dot"),
            showlegend=True,
        )
    )
    fig.add_annotation(
        x=10,
        y=-1.5,
        yanchor="top",
        text=f"µ: {r:.2f} 1/h",
        showarrow=False,
    )
    fig.update_layout(
        xaxis=dict(title="Time [h]"),
        yaxis=dict(title="OD", type="log"),
        width=width,
        height=height,
        showlegend=False,
    )
    fig = style_plot(
        fig,
        line_thickness=1,
        buttom_margin=0,
        left_margin=40,
        top_margin=0,
        right_margin=0,
    )
    fig.write_image("../plots/experiments/st_fit_conc_1.svg")


def st_025():
    conc = 0.25
    lgs = get_lgs("Salmonella Typhimurium mcherry:pGDPI CTX-M-15", conc)
    xs = [raw[raw["linegroup"] == lg]["time"].to_numpy() for lg in lgs]
    ys = [raw[raw["linegroup"] == lg]["measurement"].to_numpy() for lg in lgs]
    y_diff = [np.array([min(y) + 1e-6]) for y in ys]
    ys = [y - yd for y, yd in zip(ys, y_diff)]
    xs, ys = filter_time(
        xs,
        ys,
    )
    r = 1.4
    sim = odeint(logistic, [np.average([y[0] for y in ys])], xs[0], args=(r, 1.1))[:, 0]
    y = np.mean(ys, axis=0)
    x = np.mean(xs, axis=0)  # if all xs are identical, you could just use xs[0]
    y_lower = np.percentile(ys, 16, axis=0)
    y_upper = np.percentile(ys, 84, axis=0)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=np.concatenate([x, x[::-1]]),
            y=np.concatenate([y_upper, y_lower[::-1]]),
            fill="toself",
            fillcolor=colors["st_error"],
            line=dict(color=colors["st_error"]),
            hoverinfo="skip",
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name="<i>Salmonella Typhimurium</i> <br>mcherry:pGDPI CTX-M-15",
            line=dict(color=colors["st"]),
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=[s for s in sim if s <= max(y)],
            name="Logistic fit",
            mode="lines",
            line=dict(color=colors["st"], dash="dot"),
            showlegend=True,
        )
    )
    fig.add_annotation(
        x=10,
        y=-1.5,
        yanchor="top",
        text=f"µ: {r:.2f} 1/h",
        showarrow=False,
    )
    fig.update_layout(
        xaxis=dict(title="Time [h]"),
        yaxis=dict(title="OD", type="log"),
        width=width,
        height=height,
        showlegend=False,
    )
    fig = style_plot(
        fig,
        line_thickness=1,
        buttom_margin=0,
        left_margin=40,
        top_margin=0,
        right_margin=0,
    )
    fig.write_image("../plots/experiments/st_fit_conc_025.svg")


def ecoli_1():
    conc = 1
    lgs = get_lgs("Escherichia coli  sfGFP:pGDPI CAT", conc)
    xs = [raw[raw["linegroup"] == lg]["time"].to_numpy() for lg in lgs]
    ys = [raw[raw["linegroup"] == lg]["measurement"].to_numpy() for lg in lgs]
    y_diff = [np.array([min(y) + 1e-6]) for y in ys]
    ys = [y - yd for y, yd in zip(ys, y_diff)]
    xs, ys = filter_time(
        xs,
        ys,
    )
    r = 1.8
    sim = odeint(logistic, [np.average([y[0] for y in ys])], xs[0], args=(r, 1.1))[:, 0]
    y = np.mean(ys, axis=0)
    x = np.mean(xs, axis=0)  # if all xs are identical, you could just use xs[0]
    y_lower = np.percentile(ys, 16, axis=0)
    y_upper = np.percentile(ys, 84, axis=0)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=np.concatenate([x, x[::-1]]),
            y=np.concatenate([y_upper, y_lower[::-1]]),
            fill="toself",
            fillcolor=colors["ecoli_error"],
            line=dict(color=colors["ecoli_error"]),
            hoverinfo="skip",
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name="<i>Escherichia coli</i> <br>sfGFP:pGDPI CAT",
            line=dict(color=colors["ecoli"]),
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=[s for s in sim if s <= max(y)],
            name="Logistic fit",
            mode="lines",
            line=dict(color=colors["ecoli"], dash="dot"),
            showlegend=True,
        )
    )
    fig.add_annotation(
        x=10,
        y=-1.5,
        yanchor="top",
        text=f"µ: {r:.2f} 1/h",
        showarrow=False,
    )
    fig.update_layout(
        xaxis=dict(title="Time [h]"),
        yaxis=dict(title="OD", type="log"),
        width=width,
        height=height,
        showlegend=False,
    )
    fig = style_plot(
        fig,
        line_thickness=1,
        buttom_margin=0,
        left_margin=40,
        top_margin=0,
        right_margin=0,
    )
    fig.write_image("../plots/experiments/ecoli_fit_conc_1.svg")


def ecoli_025():
    conc = 0.25
    lgs = get_lgs("Escherichia coli  sfGFP:pGDPI CAT", conc)
    xs = [raw[raw["linegroup"] == lg]["time"].to_numpy() for lg in lgs]
    ys = [raw[raw["linegroup"] == lg]["measurement"].to_numpy() for lg in lgs]
    y_diff = [np.array([min(y) + 1e-6]) for y in ys]
    ys = [y - yd for y, yd in zip(ys, y_diff)]
    xs, ys = filter_time(
        xs,
        ys,
    )
    r = 1.3
    sim = odeint(logistic, [np.average([y[0] for y in ys])], xs[0], args=(r, 1.1))[:, 0]
    y = np.mean(ys, axis=0)
    x = np.mean(xs, axis=0)  # if all xs are identical, you could just use xs[0]
    y_lower = np.percentile(ys, 16, axis=0)
    y_upper = np.percentile(ys, 84, axis=0)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=np.concatenate([x, x[::-1]]),
            y=np.concatenate([y_upper, y_lower[::-1]]),
            fill="toself",
            fillcolor=colors["ecoli_error"],
            line=dict(color=colors["ecoli_error"]),
            hoverinfo="skip",
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name="<i>Escherichia coli</i> <br>sfGFP:pGDPI CAT",
            line=dict(color=colors["ecoli"]),
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=[s for s in sim if s <= max(y)],
            name="Logistic fit",
            mode="lines",
            line=dict(color=colors["ecoli"], dash="dot"),
            showlegend=True,
        )
    )
    fig.add_annotation(
        x=10,
        y=-1.5,
        yanchor="top",
        text=f"µ: {r:.2f} 1/h",
        showarrow=False,
    )
    fig.update_layout(
        xaxis=dict(title="Time [h]"),
        yaxis=dict(title="OD", type="log"),
        width=width,
        height=height,
        showlegend=False,
    )
    fig = style_plot(
        fig,
        line_thickness=1,
        buttom_margin=0,
        left_margin=40,
        top_margin=0,
        right_margin=0,
    )
    fig.write_image("../plots/experiments/ecoli_fit_conc_025.svg")


st_1()
st_025()
ecoli_1()
ecoli_025()
