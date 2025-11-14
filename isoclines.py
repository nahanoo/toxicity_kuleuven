import plotly.graph_objects as go
import numpy as np
from style import *
from model import Experiment
from joblib import Parallel, delayed
import sys

sys.path.append("/users/eulrich/toxicity_kuleuven")

coolwarm_colorscale = [
    [0, "#316142"],
    [1, "#98caaa"],
]


def net_growth_cefo():
    e = Experiment()
    J = np.linspace(0, e.E.r, 100)
    Cf = np.linspace(0, e.M_cefo, 100)
    J_grid, Cf_grid = np.meshgrid(np.array(J), np.array(Cf))
    JD = J_grid - e.E.u * Cf_grid / (Cf_grid + e.E.KT)
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
        width=180,
        height=180,
        showlegend=False,
    )
    fig = style_plot(fig, line_thickness=1, left_margin=30, font_size=11)
    fig.write_image("plots/isoclines/cefotaxime.svg")


def net_growth_chloramphenicol():
    e = Experiment()
    fig = go.Figure()
    J = np.linspace(0, e.S.r, 100)
    chloram = np.linspace(0, e.M_chloram, 100)
    J_grid, chloram_grid = np.meshgrid(np.array(J), np.array(chloram))
    JD = J_grid * 1 / (1 + chloram_grid / e.S.KT)
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
            showscale=True,
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
        width=180,
        height=180,
    )
    fig = style_plot(fig, line_thickness=1, left_margin=30, font_size=11)
    fig.write_image("plots/isoclines/chloramphenicol.svg")


def sweep_growth_rates():
    EXTINCT = 1e5

    def run_one(r_e, r_s):
        e = Experiment()  # construct in worker
        e.E.r, e.S.r = r_e, r_s
        e.E.u = 2.1
        e.M_cefo = 0
        # e.M_chloram = 0
        e.total_transfers = 4
        e.transfer()
        E = e.E.N[-1]
        S = e.S.N[-1]
        # extinction threshold
        E = None if E < EXTINCT else E
        S = None if S < EXTINCT else S
        if (E is None) and (S is not None):
            ratio = 0
        elif (E is not None) and (S is None):
            ratio = 1
        elif (E is None) and (S is None):
            ratio = None
        else:
            ratio = E / (E + S)
        return E, S, ratio

    rs_E = np.linspace(1, 1.8, 20)
    rs_S = np.linspace(0, 1.8, 20)
    n = len(rs_E)

    results = Parallel(
        n_jobs=-1, backend="loky", prefer="processes", batch_size="auto"
    )(delayed(run_one)(r_e, r_s) for r_e in rs_E for r_s in rs_S)

    E_flat, S_flat, R_flat = map(np.array, zip(*results))
    Z_E = E_flat.reshape(n, n).astype(np.float32)
    Z_S = S_flat.reshape(n, n).astype(np.float32)
    Z_ratio = R_flat.reshape(n, n).astype(np.float32)
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            z=Z_ratio,
            x=rs_S,
            y=rs_E,
            contours=dict(showlines=False),
            colorscale=coolwarm_colorscale,
            zmin=0,
            zmax=1,
            zauto=False,
            ncontours=50,
            colorbar=dict(len=0.5, y=0.2, thickness=12),
            showscale=False,
        )
    )
    d_x = rs_S[rs_S >= 1]
    d_y = rs_S[rs_S >= 1]
    fig.add_trace(
        go.Scatter(
            x=d_x,
            y=d_y,
            mode="lines",
            line=dict(color="black", width=1.5, dash="dash"),
            showlegend=False,
        )
    )

    fig.update_layout(
        xaxis=dict(
            title="Maximum growth rate ST (1/h)", ticks="inside", showgrid=False
        ),
        yaxis=dict(
            title="Maximum growth rate EcN (1/h)", ticks="inside", showgrid=False
        ),
        width=180,
        height=180,
    )
    fig = style_plot(fig, line_thickness=1.5, font_size=11)
    fig.write_image("plots/contours/s_protects_e/sweep_growth_rates.svg")


def sweep_death_vs_degradation():

    EXTINCT = 1e5

    def run_one(beta, delta):
        e = Experiment()  # construct in worker
        # e.M_chloram = 0
        e.E.r, e.S.r = 1.8, 1.5
        e.E.u = 1
        e.E.a = beta
        e.S.a = beta
        e.transfer()
        E = e.E.N[-1]
        S = e.S.N[-1]
        E = None if E < EXTINCT else E
        S = None if S < EXTINCT else S
        if (E is None) and (S is not None):
            ratio = 0
        elif (E is not None) and (S is None):
            ratio = 1
        elif (E is None) and (S is None):
            ratio = None
        else:
            ratio = E / (E + S)
        return E, S, ratio

    # Degradation rates
    n = 100
    betas = np.linspace(0, 1, n)

    # Death rates
    deltas = np.linspace(0, 2.5, n)

    results = Parallel(
        n_jobs=-1, backend="loky", prefer="processes", batch_size="auto"
    )(delayed(run_one)(beta, delta) for beta in betas for delta in deltas)

    E_flat, S_flat, R_flat = map(np.array, zip(*results))
    Z_E = E_flat.reshape(n, n).astype(np.float32)
    Z_S = S_flat.reshape(n, n).astype(np.float32)
    Z_ratio = R_flat.reshape(n, n).astype(np.float32)
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            z=Z_ratio,
            x=betas,
            y=betas,
            contours=dict(showlines=False),
            colorscale=coolwarm_colorscale,
            zmin=0,
            zmax=1,
            zauto=False,
            ncontours=50,
            colorbar=dict(len=0.5, y=0.2, thickness=12),
            showscale=False,
        )
    )
    fig.update_layout(
        xaxis=dict(title="Death rate EcN (1/h)", ticks="inside"),
        yaxis=dict(title="Degradation rate ST (1/h)", ticks="inside", showgrid=False),
        width=180,
        height=180,
    )
    fig = style_plot(fig, line_thickness=1.5, font_size=11)
    fig.write_image("plots/contours/two_sided/sweep_death_vs_degradation.svg")


sweep_death_vs_degradation()
