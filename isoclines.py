import plotly.graph_objects as go
import numpy as np
from style import *
from model import Experiment

coolwarm_colorscale = [
    [0.0, "#3B4CC0"],
    [0.25, "#8CABD8"],
    [0.5, "#FFFFFF"],
    [0.75, "#D9907E"],
    [1.0, "#B40426"],
]


def cefotaxime():
    e = Experiment()
    J = np.linspace(0, e.E.r, 100)
    Cf = np.linspace(0, e.M_cefo, 100)
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
                tickfont=dict(size=11),
            ),
            showscale=False,
            contours=dict(showlines=False, start=-0.8, end=0.8, size=0.2),
        )
    )
    fig.update_layout(
        xaxis=dict(
            # range=[0, max(J)],
            # dtick=0.1,
            title="µ [1/h]",
            ticks="inside",
        ),
        yaxis=dict(
            ##range=[0, max(Cf)],
            # dtick=0.05,
            title="Cefotaxime [μg/mL]",
            ticks="inside",
        ),
        # title="J = -1 1/h",
        width=width,
        height=height,
        showlegend=False,
    )
    fig = style_plot(fig, line_thickness=1, left_margin=30, font_size=11)
    fig.write_image("plots/isoclines/cefotaxime.svg")


def chloramphenicol():
    e = Experiment()
    fig = go.Figure()
    J = np.linspace(0, e.S.r, 100)
    chloram = np.linspace(0, e.M_chloram, 100)
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
                tickfont=dict(size=11),
            ),
            contours=dict(showlines=False),
            showscale=False,
        )
    )
    fig.update_layout(
        xaxis=dict(
            title="µ [1/h]",
            ticks="inside",
        ),
        yaxis=dict(
            title="Chloram [μg/mL]",
            ticks="inside",
        ),
        # title="J = -1 1/h",
        width=width,
        height=height,
    )
    fig = style_plot(fig, line_thickness=1, left_margin=30, font_size=11)
    fig.write_image("plots/isoclines/chloramphenicol.svg")


chloramphenicol()
