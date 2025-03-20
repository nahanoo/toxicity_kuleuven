import pandas as pd
import numpy as np
import plotly.graph_objects as go
from style import *
from scipy.optimize import curve_fit
from scipy.stats import linregress

colors["EcN sfGFP"] = colors["ecoli"]
colors["EcN sfGFP CAT"] = colors["st"]
colors["ST mCherry"] = "#007c00"


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
    fig.write_image("../plots/experiments/death_rate_ecoli_sfGFP_experiment2.svg")

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
    fig.write_image("../plots/experiments/death_rate_ecoli_sfGFP_CAT_experiment2.svg")
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


def plot_death_experiment_linregress():
    repeats = ["Repeat_1", "Repeat_2", "Repeat_3"]

    strains = ["EcN sfGFP", "ST mCherry", "EcN sfGFP CAT"]

    f = "../data/cefotaxime_death_rate/Cefotaxime_death_rate_2.xlsx"

    data = {strain: [] for strain in strains}

    for repeat in repeats:
        df = pd.read_excel(f, sheet_name=repeat, index_col=0, header=None)
        for strain in strains:
            data[strain].append(df.loc[strain].to_numpy())
    strain = "EcN sfGFP"
    fig = go.Figure()
    showlegend = True
    for i, repeat in enumerate(repeats):
        fig.add_trace(
            go.Scatter(
                x=minutes / 60,
                y=data[strain][i],
                showlegend=showlegend,
                line=dict(color=colors["ecoli"]),
                name="Data",
            )
        )
        showlegend = False
    m, b = linregress(minutes / 60, np.log(np.nanmean(data[strain], axis=0)))[:2]
    print(m)
    fit = np.exp(np.array([m * t + b for t in minutes / 60]))
    fig.add_trace(
        go.Scatter(
            x=minutes / 60, y=fit, line=dict(color="black", dash="dash"), name="Model"
        )
    )
    fig.update_layout(
        width=width / 1.3, height=height / 1.3, title="EcN sfGFP CAT Experiment 2"
    )
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="CFUs/mL", type="log")
    fig = style_plot(fig, line_thickness=2, marker_size=4, left_margin=60)
    fig.write_image("../plots/experiments/death_rate_e_coli_sfGFP_CAT_fit.svg")


def plot_all_strains():
    fig = go.Figure()
    for strain in strains:
        strain_legend_dashed = True  # Legend flag for dashed lines
        strain_legend_solid = True  # Legend flag for solid lines

        for repeat, repeat_name in zip(data[strain], repeats):
            # Determine line style based on repeat name
            if ("1" in repeat_name) | ("2" in repeat_name) | ("3" in repeat_name):
                line_dash = "dash"
                show_legend = strain_legend_dashed  # Show legend for dashed line
                strain_legend_dashed = False  # Ensure it's only shown once
                legend_group = f"{strain}_batch1"  # Group by strain and batch
                legend_name = f"{strain} (Batch 1)"  # Add batch name to legend
            else:
                line_dash = "solid"
                show_legend = strain_legend_solid  # Show legend for solid line
                strain_legend_solid = False  # Ensure it's only shown once
                legend_group = f"{strain}_batch2"  # Group by strain and batch
                legend_name = f"{strain} (Batch 2)"  # Add batch name to legend

            # Add trace
            fig.add_trace(
                go.Scatter(
                    x=minutes / 60,
                    y=repeat,
                    line=dict(color=colors[strain], dash=line_dash),
                    showlegend=show_legend,
                    name=legend_name,  # Use batch-specific name
                    legendgroup=legend_group,  # Group traces by batch and strain
                )
            )

    fig.update_yaxes(title="CFUs/mL", type="log")
    fig.update_xaxes(title="Time [h]")
    fig = style_plot(fig)
    fig.update_layout(width=width * 1.5, height=height)
    fig.write_image("../plots/experiments/death_rate_all_strains.svg")
