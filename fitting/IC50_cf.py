import pandas as pd
import numpy as np
import plotly.graph_objects as go
from style import *
from scipy.optimize import curve_fit

concentrations = np.array(
    [
        0,
        0.0015625,
        0.003125,
        0.0045875,
        0.00625,
        0.015625,
        0.03125,
        0.045875,
        0.0625,
        0.125,
        0.25,
    ]
)
strains = [
    "EcN sfGFP",
    "EcN tdTomato",
    "EcN tdTomato CTX",
    "EcN sfGFP CAT",
    "ST sfGFP",
    "ST mCherry",
    "ST mCherry CTX",
    "ST sfGFP CAT",
]

data = {key: [] for key in strains}
response = {key: [] for key in strains}

repeats = ["Repeat_1", "Repeat_2", "Repeat_3"]
f = "../data/IC50_cefotaxime/Processed_IC50_data_Cf.xlsx"
for repeat in repeats:
    df = pd.read_excel(f, sheet_name=repeat, skiprows=[1], index_col="Strain")
    for strain in strains:
        data[strain].append(df.loc[strain].to_numpy())
        cf_0 = df.loc[strain].iloc[0]
        inhibition = cf_0 - df.loc[strain].to_numpy()[1:]
        inhibition = inhibition / cf_0
        response[strain].append(inhibition)


def st_resistance():
    strain = "ST mCherry CTX"
    fig = go.Figure()
    for i, repeat in enumerate(repeats):
        fig.add_trace(
            go.Scatter(
                x=concentrations,
                y=data[strain][i],
                mode="lines+markers",
                marker=dict(color=colors["st"]),
                showlegend=False,
            )
        )
    fig.update_xaxes(title="Cefotaxime [µL/mL]"), fig.update_yaxes(
        title="Abundance [OD]"
    )
    fig.update_layout(
        width=width / 1.3, height=height / 1.3, title="ST mCherry CTX resistance"
    )
    fig = style_plot(fig, marker_size=5, line_thickness=2, left_margin=60)
    fig.write_image("../figures/fig2_st_cefotaxime_resistance.svg")


def e_coli_sus():
    strain = "EcN sfGFP CAT"
    fig = go.Figure()
    for i, repeat in enumerate(repeats):
        fig.add_trace(
            go.Scatter(
                x=concentrations,
                y=data[strain][i],
                mode="lines+markers",
                marker=dict(color=colors["ecoli"]),
                showlegend=False,
            )
        )
    fig.add_vline(x=0.0625, annotation=dict(text="MIC"), line_color="black")
    fig.update_xaxes(title="Cefotaxime [µL/mL]"), fig.update_yaxes(
        title="Abundance [OD]"
    )
    fig.update_layout(
        width=width / 1.3, height=height / 1.3, title="EcN sfGFP CAT Susceptibility"
    )
    fig = style_plot(fig, marker_size=5, line_thickness=2, left_margin=60)
    fig.write_image("../figures/fig2_e_coli_susceptibility.svg")


def e_coli_response_curve():
    strain = "EcN sfGFP CAT"
    fig = go.Figure()
    for i, repeat in enumerate(repeats):
        fig.add_trace(
            go.Scatter(
                x=concentrations[1:],
                y=response[strain][i],
                mode="markers",
                marker=dict(color=colors["ecoli"]),
                showlegend=False,
            )
        )
    fig.update_xaxes(title="Cefotaxime [µL/mL]", type="log"), fig.update_yaxes(
        title="Response in %"
    )
    fig.update_layout(
        width=width / 1.3, height=height / 1.3, title="EcN sfGFP CAT Response curve"
    )
    fig = style_plot(fig, marker_size=5, line_thickness=2, left_margin=60)
    fig.write_image("../figures/fig2_e_coli_response_curve.svg")


def logistic_curve(log_dose, top, bottom, log_ic50, hill_slope):
    return bottom + (top - bottom) / (1 + 10 ** ((log_ic50 - log_dose) * hill_slope))


def fit_IC50():
    strain = "EcN sfGFP CAT"
    initial_guess = [1, 0, -4, 1]
    params, covariance = curve_fit(
        logistic_curve,
        np.log(concentrations[1:]),
        np.average(response[strain], axis=0),
        p0=initial_guess,
    )

    top, bottom, log_ic50, hill_slope = params
    ic50 = np.exp(log_ic50)  # Convert log_IC50 back to linear scale
    print(ic50)
    log_dose_fit = np.linspace(
        np.log(concentrations[1:]).min(), np.log(concentrations[1:]).max(), 100
    )
    response_fit = logistic_curve(log_dose_fit, *params)
    fig = go.Figure()
    show_legend = True
    for i, repeat in enumerate(repeats):
        fig.add_trace(
            go.Scatter(
                x=concentrations[1:],
                y=response[strain][i],
                mode="markers",
                marker=dict(color=colors["ecoli"]),
                name="Data",
                showlegend=show_legend,
            )
        )
        show_legend = False
    fig.add_trace(
        go.Scatter(
            x=np.exp(log_dose_fit),
            y=response_fit,
            line=dict(color=colors["ecoli"], dash="dash"),
            name="Model",
        )
    )
    fig.add_vline(
        x=ic50, line_color="black", line_dash="dash", annotation=dict(text="IC50")
    )
    fig.update_xaxes(title="Cefotaxime [µL/mL]", type="log"), fig.update_yaxes(
        title="Response in %"
    )
    fig.update_layout(
        width=width / 1.3, height=height / 1.3, title="EcN sfGFP CAT Response curve"
    )
    fig = style_plot(fig, marker_size=5, line_thickness=2, left_margin=60)
    fig.write_image("../figures/fig2_e_coli_response_curve_fit.svg")


fit_IC50()
