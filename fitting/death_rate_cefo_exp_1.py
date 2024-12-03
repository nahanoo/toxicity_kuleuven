import pandas as pd
import numpy as np
import plotly.graph_objects as go
from style import *
from scipy.optimize import curve_fit

repeats = ["Repeat_1", "Repeat_2", "Repeat_3"]
minutes = np.array([0, 15, 30, 45, 60, 90, 120, 150, 180, 210, 240])

strains = ["EcN sfGFP", "EcN sfGFP CAT", "ST sfGFP CAT"]

f = "../data/cefotaxime_death_rate/Cefotaxime_death_rate_1.xlsx"

data = {strain: [] for strain in strains}

for repeat in repeats:
    df = pd.read_excel(
        f,
        sheet_name=repeat,
        index_col=0,
    )
    for strain in strains:
        data[strain].append(df.loc[strain].to_numpy())


def plot_death_experiment_1():
    strain = "EcN sfGFP"
    fig = go.Figure()
    for i, repeat in enumerate(repeats):
        fig.add_trace(
            go.Scatter(
                x=minutes / 60,
                y=data[strain][i],
                showlegend=False,
                line=dict(color=colors["ecoli"]),
            )
        )
    fig.update_layout(
        width=width / 1.3, height=height / 1.3, title="EcN sfGFP Experiment 1"
    )
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="CFUs/mL")
    fig = style_plot(fig, line_thickness=2, marker_size=4, left_margin=60)
    fig.write_image("../figures/fig4_ecoli_sfGFP_experiment1.svg")

    strain = "EcN sfGFP CAT"
    fig = go.Figure()
    for i, repeat in enumerate(repeats):
        fig.add_trace(
            go.Scatter(
                x=minutes / 60,
                y=data[strain][i],
                showlegend=False,
                line=dict(color=colors["ecoli"]),
            )
        )
    fig.update_layout(
        width=width / 1.3, height=height / 1.3, title="EcN sfGFP CAT Experiment 1"
    )
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="CFUs/mL")
    fig = style_plot(fig, line_thickness=2, marker_size=4, left_margin=60)
    fig.write_image("../figures/fig4_ecoli_sfGFP_CAT_experiment1.svg")
