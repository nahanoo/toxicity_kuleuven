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


def four_param_logistic(dose, top, bottom, ic50, hill_slope):
    dose = np.asarray(dose, dtype=float)
    return bottom + (top - bottom) / (1.0 + (ic50 / dose) ** hill_slope)


strain = "EcN sfGFP CAT"
concs = concentrations[1:]  # exclude zero concentration
mean_response = np.mean(response[strain], axis=0)
# Initial guess
top0 = mean_response.max()
bottom0 = mean_response.min()
half_max = (top0 + bottom0) / 2

# Concentration closest to half-max as IC50 guess
ic50_0 = concs[np.argmin(np.abs(mean_response - half_max))]
p0 = [top0, bottom0, ic50_0, 1.0]

bounds = (
    [0.5, -0.6, concs.min() / 10, 0.05],  # top_lo, bottom_lo, ic50_lo, hill_lo
    [1.5, 0.5, concs.max() * 10, 5.0],
)

params, covariance = curve_fit(
    four_param_logistic,
    concs,
    mean_response,
    p0=p0,
    bounds=bounds,
    maxfev=100000,
)

top, bottom, ic50, hill_slope = params
print(
    f"top = {top:.2f}, bottom = {bottom:.2f}, IC50 = {ic50:.3f} µg/mL, hill = {hill_slope:.2f}"
)

x_fit = np.logspace(np.log10(concs.min()), np.log10(concs.max()), 400)
y_fit = four_param_logistic(x_fit, *params)

fig = go.Figure()
for resp in response[strain]:
    fig.add_trace(
        go.Scatter(
            x=concs,
            y=resp,
            mode="markers",
            marker=dict(color="#888578"),
            showlegend=False,
        )
    )

fig.add_trace(
    go.Scatter(
        x=x_fit,
        y=y_fit,
        mode="lines",
        name=f"4PL fit (IC50 ≈ {ic50:.2f} µg/mL)",
        line=dict(color="#888578"),
        showlegend=False,
    )
)

fig.update_layout(
    xaxis=dict(title="Cefotaxime [µg/mL]", type="log", tickformat="E", tickangle=45),
    yaxis=dict(title="Inhibition [%]"),
    width=250,
    height=200,
    title="IC50 for cefotaxime",
)
fig.add_vline(
    x=ic50,
    line_color="#888578",
    line_dash="dot",
    line_width=2,
    annotation=dict(
        text=f"IC50 = {ic50:.2f} µg/mL",
        showarrow=False,
        yref="y domain",
        xref="x domain",
        y=0.9,
        x=0.6,
        xanchor="right",
        yanchor="top",
    ),
)
fig = style_plot(fig, marker_size=6, line_thickness=2, font_size=12)
fig.write_image("../plots/experiments/cefotaxime_dose_response.svg")
