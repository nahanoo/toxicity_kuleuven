import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import plotly.graph_objects as go
from style import *


# --- 4-parameter logistic in concentration space ---
# dose: chloramphenicol concentration [µg/mL]
# top: maximum inhibition (high dose)
# bottom: minimum inhibition (low dose)
# ic50: concentration where response = (top + bottom)/2
# hill_slope: steepness of the transition
def four_param_logistic(dose, top, bottom, ic50, hill_slope):
    dose = np.asarray(dose, dtype=float)
    return bottom + (top - bottom) / (1.0 + (ic50 / dose) ** hill_slope)


# --- Load data ---
f = "../data/IC50_chloram/20251205_Cm2.xlsx"
df = pd.read_excel(f, sheet_name="parsing")

# concentrations are in the header from col 2 onwards
concs = df.columns[2:].astype(float).to_numpy()

# mean percentage inhibition across replicates
y_mat = df.iloc[:, 2:].to_numpy(dtype=float)
mean_response = y_mat.mean(axis=0)

# --- Initial guesses ---
top0 = mean_response.max()
bottom0 = mean_response.min()
half_max = (top0 + bottom0) / 2

# concentration closest to half-max as IC50 guess
ic50_0 = concs[np.argmin(np.abs(mean_response - half_max))]
p0 = [top0, bottom0, ic50_0, 1.0]

# optional parameter bounds (keeps fit sane)
bounds = (
    [0.0, 0.0, concs.min() / 10, 0.01],  # lower
    [120.0, 120.0, concs.max() * 10, 5.0],  # upper
)

# --- Fit ---
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

# --- Smooth curve for plotting ---
x_fit = np.logspace(np.log10(concs.min()), np.log10(concs.max()), 400)
y_fit = four_param_logistic(x_fit, *params)

# --- Plot ---
fig = go.Figure()

# raw replicates
for _, row in df.iterrows():
    fig.add_trace(
        go.Scatter(
            x=concs,
            y=row.iloc[2:].to_numpy(dtype=float),
            mode="markers",
            marker=dict(color=colors["light_black"]),
            showlegend=False,
        )
    )

# fitted curve
fig.add_trace(
    go.Scatter(
        x=x_fit,
        y=y_fit,
        mode="lines",
        name=f"4PL fit (IC50 ≈ {ic50:.2f} µg/mL)",
        line=dict(color=colors["light_black"]),
        showlegend=False,
    )
)

fig.update_layout(
    xaxis=dict(title="Chloramphenicol [µg/mL]", type="log"),
    yaxis=dict(title="Inhibition [%]"),
)
fig.add_vline(
    x=ic50,
    line_color=colors["light_black"],
    line_dash="dash",
    line_width=1,
    annotation=dict(text=f"IC50 = {ic50:.2f} µg/mL"),
    # annotation_position="center left",
)
fig = style_plot(fig, marker_size=8, line_thickness=2)
fig.write_image("../plots/experiments/chloramphenicol_dose_response_redo.svg")
