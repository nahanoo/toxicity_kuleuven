import pandas as pd
import numpy as np
import plotly.graph_objects as go
from style import *
from scipy.optimize import curve_fit

repeats = ["Repeat_1", "Repeat_2", "Repeat_3", "Repeat_4", "Repeat_5", "Repeat_6"]
minutes = np.array([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360])

strains = ["EcN sfGFP", "ST mCherry", "EcN sfGFP CAT"]

f = "../data/cefotaxime_death_rate/Cefotaxime_death_rate_2.xlsx"

data = {strain: [] for strain in strains}

for repeat in repeats:
    df = pd.read_excel(f, sheet_name=repeat, index_col=0, header=None)
    for strain in strains:
        data[strain].append(df.loc[strain].to_numpy())


def death_model(t, N0, k):
    return N0 * np.exp(-k * t)


def plot_death_experiment_2():
    strain = "EcN sfGFP"
    fig = go.Figure()
    showlegend = False
    for i, repeat in enumerate(repeats):
        fig.add_trace(
            go.Scatter(
                x=minutes / 60,
                y=data[strain][i],
                showlegend=showlegend,
                line=dict(color=colors["ecoli"]),
            )
        )
        showlegend = False

    fig.update_layout(
        width=width / 1.3, height=height / 1.3, title="EcN sfGFP Experiment 2"
    )
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="CFUs/mL", type="log")
    fig = style_plot(fig, line_thickness=2, marker_size=4, left_margin=60)
    fig.write_image("../figures/fig4_ecoli_sfGFP_experiment2.svg")

    strain = "EcN sfGFP CAT"
    fig = go.Figure()
    showlegend = False
    for i, repeat in enumerate(repeats):
        fig.add_trace(
            go.Scatter(
                x=minutes / 60,
                y=data[strain][i],
                showlegend=showlegend,
                line=dict(color=colors["ecoli"]),
            )
        )
        showlegend = False

    fig.update_layout(
        width=width / 1.3, height=height / 1.3, title="EcN sfGFP CAT Experiment 2"
    )
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="CFUs/mL", type="log")
    fig = style_plot(fig, line_thickness=2, marker_size=4, left_margin=60)
    fig.write_image("../figures/fig4_ecoli_sfGFP_CAT_experiment2.svg")
    strain = "ST mCherry"
    fig = go.Figure()
    showlegend = False
    for i, repeat in enumerate(repeats):
        fig.add_trace(
            go.Scatter(
                x=minutes / 60,
                y=data[strain][i],
                showlegend=showlegend,
                line=dict(color=colors["st"]),
            )
        )
        showlegend = False

    fig.update_layout(
        width=width / 1.3, height=height / 1.3, title="ST mCherry Experiment 2"
    )
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="CFUs/mL", type="log")
    fig = style_plot(fig, line_thickness=2, marker_size=4, left_margin=60)
    fig.write_image("../figures/fig4_st_mCherry_experiment2.svg")


plot_death_experiment_2()
