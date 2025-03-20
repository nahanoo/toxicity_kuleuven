import plotly.graph_objects as go
import numpy as np
from style import *

# Define the coolwarm colorscale explicitly
coolwarm_colorscale = [
    [0.0, "rgb(58, 76, 192)"],  # Dark blue
    [0.5, "rgb(221, 221, 221)"],  # White
    [1.0, "rgb(180, 4, 38)"],  # Dark red
]

J = np.linspace(0, 0.9, 100)
Cf = np.linspace(0, 0.25, 100)
J_grid, Cf_grid = np.meshgrid(np.array(J), np.array(Cf))
JD = J_grid - 1 * Cf_grid / (Cf_grid + 0.029)
fig = go.Figure()
fig.add_trace(
    go.Contour(
        x=J,
        y=Cf,
        z=JD,
        colorscale=coolwarm_colorscale,
        colorbar=dict(
            title="Netto growth rate [1/h]",
            title_side="right",
            dtick=0.2,
            tickfont=dict(size=font_size),
        ),
        contours=dict(showlines=False, start=-0.8, end=0.8, size=0.2),
    )
)
fig.update_layout(
    xaxis=dict(
        range=[0, max(J)],
        dtick=0.1,
        title="Per capita growth rate µ [1/h]",
    ),
    yaxis=dict(
        range=[0, max(Cf)],
        dtick=0.05,
        title="Cefotaxime concentration μg/mL",
    ),
    title="J = -1 1/h",
    width=width,
    height=height,
)
fig = style_plot(fig, line_thickness=1, left_margin=50)
fig.write_image("plots/isoclines/cefotaxime.svg")

fig = go.Figure()
J = np.linspace(0, 0.7, 100)
chloram = np.linspace(0, 16, 100)
J_grid, chloram_grid = np.meshgrid(np.array(J), np.array(chloram))
JD = J_grid * 1 / (1 + chloram_grid / 0.32)
fig = go.Figure()
fig.add_trace(
    go.Contour(
        x=J,
        y=chloram,
        z=JD,
        colorscale=coolwarm_colorscale,
        colorbar=dict(
            title="Netto growth rate [1/h]",
            title_side="right",
            dtick=0.2,
            tickfont=dict(size=font_size),
        ),
        contours=dict(showlines=False, start=0, end=0.8, size=0.2),
    )
)
fig.update_layout(
    xaxis=dict(
        range=[0, max(J)],
        dtick=0.1,
        title="Per capita growth rate µ [1/h]",
    ),
    yaxis=dict(
        range=[0, max(chloram)],
        dtick=2,
        title="Chloramphenicol concentration μg/mL",
    ),
    title="J = -1 1/h",
    width=width,
    height=height,
)
fig = style_plot(fig, line_thickness=1, left_margin=50)
fig.write_image("plots/isoclines/chloramphenicol.svg")
