import plotly.graph_objects as go
import numpy as np
from style import *
from model import Experiment
from joblib import Parallel, delayed
from parameters import *


def run_grid_xy(fun, xs, ys):
    results = Parallel(
        n_jobs=-1, backend="loky", prefer="processes", batch_size="auto"
    )(delayed(fun)(x, y) for y in ys for x in xs)
    E, S = map(np.array, zip(*results))
    E = E.astype(float).reshape(len(ys), len(xs))
    S = S.astype(float).reshape(len(ys), len(xs))
    return E, S


def outcome(E, S, thr=1e5):
    alive_E = E > thr
    alive_S = S > thr
    both_alive = alive_E & alive_S
    E_excludes = alive_E & (~alive_S)
    S_excludes = (~alive_E) & alive_S
    z = np.full_like(E, np.nan, dtype=float)
    z[both_alive] = np.log10(E[both_alive] / S[both_alive])
    finite = np.isfinite(z)
    zmin = float(np.nanmin(z)) if finite.any() else -1
    zmax = float(np.nanmax(z)) if finite.any() else 1
    return z, zmin, zmax, E_excludes, S_excludes


def add_outcome_heatmaps(fig, x, y, z, zmin, zmax, E_excludes, S_excludes):
    c_E_excludes = "#E76F51"
    c_S_excludes = "#2A9D8F"
    fig.add_trace(
        go.Heatmap(
            x=x,
            y=y,
            z=z,
            colorscale="Cividis",
            zmin=zmin,
            zmax=zmax,
            hoverongaps=False,
            colorbar=dict(
                title="log10(E/S)",
                title_side="right",
                tickfont=dict(size=11),
                y=0.5,
                len=1,
                thickness=12,
            ),
        )
    )
    fig.add_trace(
        go.Heatmap(
            x=x,
            y=y,
            z=np.where(E_excludes, 1.0, np.nan),
            colorscale=[[0, c_E_excludes], [1, c_E_excludes]],
            showscale=False,
            hoverongaps=False,
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Heatmap(
            x=x,
            y=y,
            z=np.where(S_excludes, 1.0, np.nan),
            colorscale=[[0, c_S_excludes], [1, c_S_excludes]],
            showscale=False,
            hoverongaps=False,
            hoverinfo="skip",
        )
    )
    return fig


def sweep(xs, ys, transfer, x_title, y_title, title, out, pts=None):
    E, S = run_grid_xy(transfer, xs, ys)
    z, zmin, zmax, E_excludes, S_excludes = outcome(E, S)
    fig = go.Figure()
    fig = add_outcome_heatmaps(fig, xs, ys, z, zmin, zmax, E_excludes, S_excludes)
    if pts is not None:
        fig.add_trace(
            go.Scatter(
                x=pts[0],
                y=pts[1],
                mode="markers+text",
                marker=dict(color="red", size=10),
                text=pts[2],
                textfont=dict(color="white", size=11),
                textposition="top center",
                showlegend=False,
            )
        )
    fig.update_layout(
        xaxis=dict(title=x_title, ticks="inside", showgrid=False),
        yaxis=dict(title=y_title, ticks="inside", showgrid=False),
        title=title,
        width=220,
        height=180,
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    fig = style_plot(fig, line_thickness=1.5, font_size=11, marker_size=6)
    fig.write_image(out)


def growth_rate_sweep(
    rs_E, rs_S, title, out, M_cefo=None, M_chloram=None, annotate=False
):
    def transfer(r_s, r_e):
        e = Experiment()
        e.S.r = r_s
        e.E.r = r_e
        if M_cefo is not None:
            e.M_cefo = M_cefo
        if M_chloram is not None:
            e.M_chloram = M_chloram
        e.total_transfers = 4
        e.transfer()
        return e.E.N[-1], e.S.N[-1]

    pts = None
    if annotate:
        pts = (
            [ST_100_CFA, ST_20_CFA],
            [EC_100_CFA, EC_20_CFA],
            ["100% CFA", "20% CFA"],
        )

    sweep(
        rs_S,
        rs_E,
        transfer,
        x_title="Max. growth rate ST (1/h)",
        y_title="Max. growth rate EcN (1/h)",
        title=title,
        out=out,
        pts=pts,
    )


def death_vs_st_degradation(
    r_e, r_s, e_u, s_a, title, out, M_cefo=None, M_chloram=None
):
    def transfer(u_e, a_s):
        e = Experiment()
        e.E.r = r_e
        e.S.r = r_s
        e.E.u = u_e
        e.S.a = a_s
        if M_cefo is not None:
            e.M_cefo = M_cefo
        if M_chloram is not None:
            e.M_chloram = M_chloram
        e.total_transfers = 4
        e.transfer()
        return e.E.N[-1], e.S.N[-1]

    sweep(
        e_u,
        s_a,
        transfer,
        x_title="Death rate EcN (1/h)",
        y_title="Max. degradation rate ST (1/h)",
        title=title,
        out=out,
        pts=None,
    )


def death_vs_ec_degradation(
    r_e, r_s, e_u, e_a, title, out, M_cefo=None, M_chloram=None
):
    def transfer(u_e, a_e):
        e = Experiment()
        e.E.r = r_e
        e.S.r = r_s
        e.E.u = u_e
        e.E.a = a_e
        if M_cefo is not None:
            e.M_cefo = M_cefo
        if M_chloram is not None:
            e.M_chloram = M_chloram
        e.total_transfers = 4
        e.transfer()
        return e.E.N[-1], e.S.N[-1]

    sweep(
        e_u,
        e_a,
        transfer,
        x_title="Death rate EcN (1/h)",
        y_title="Max. degradation rate EcN (1/h)",
        title=title,
        out=out,
        pts=None,
    )


def plot_all():
    rs_01 = np.linspace(0, 1.0, 20)
    rs_12 = np.linspace(0, 1.2, 20)
    e_u = np.linspace(0, 2.0, 20)
    s_a = np.linspace(0, 1.0, 20)
    e_a = np.linspace(0, 1.0, 20)

    growth_rate_sweep(
        rs_01,
        rs_01,
        title="No antibiotics: final ratio EcN/ST",
        out="plots/contours/no_antibiotics/sweep_growth_rates.svg",
        M_cefo=0,
        M_chloram=0,
        annotate=False,
    )

    growth_rate_sweep(
        rs_12,
        rs_12,
        title="Cf",
        out="plots/contours/s_protects_e/sweep_growth_rates.svg",
        M_chloram=0,
        annotate=True,
    )

    death_vs_st_degradation(
        EC_100_CFA,
        ST_100_CFA,
        e_u,
        s_a,
        title="Cf",
        out="plots/contours/s_protects_e/death_vs_st_degradation_100CFA.svg",
        M_chloram=0,
    )
    death_vs_st_degradation(
        EC_20_CFA,
        ST_20_CFA,
        e_u,
        s_a,
        title="Cf",
        out="plots/contours/s_protects_e/death_vs_st_degradation_20CFA.svg",
        M_chloram=0,
    )

    growth_rate_sweep(
        rs_12,
        rs_12,
        title="Cm",
        out="plots/contours/e_protects_s/sweep_growth_rates.svg",
        M_cefo=0,
        annotate=True,
    )

    growth_rate_sweep(
        rs_12,
        rs_12,
        title="Cf+Cm",
        out="plots/contours/two_sided/sweep_growth_rates.svg",
        annotate=True,
    )

    death_vs_st_degradation(
        EC_100_CFA,
        ST_100_CFA,
        e_u,
        s_a,
        title="Cf+Cm",
        out="plots/contours/two_sided/death_vs_st_degradation_100CFA.svg",
    )
    death_vs_st_degradation(
        EC_20_CFA,
        ST_20_CFA,
        e_u,
        s_a,
        title="Cf+Cm",
        out="plots/contours/two_sided/death_vs_st_degradation_20CFA.svg",
    )
    death_vs_ec_degradation(
        EC_100_CFA,
        ST_100_CFA,
        e_u,
        e_a,
        title="Cf+Cm",
        out="plots/contours/two_sided/death_vs_ec_degradation_100CFA.svg",
    )
    death_vs_ec_degradation(
        EC_20_CFA,
        ST_20_CFA,
        e_u,
        e_a,
        title="Cf+Cm",
        out="plots/contours/two_sided/death_vs_ec_degradation_20CFA.svg",
    )


def net_growth_cefo():
    e = Experiment()
    J = np.linspace(0, EC_100_CFA, 100)
    Cf = np.linspace(0, e.M_cefo, 100)
    J_grid, Cf_grid = np.meshgrid(np.array(J), np.array(Cf))
    JD = J_grid - e.E.u * Cf_grid / (Cf_grid + e.E.KT)
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            x=J,
            y=Cf,
            z=JD,
            colorscale="Cividis",
            zmin=np.min(JD),
            zmax=np.max(JD),
            hoverongaps=False,
            colorbar=dict(
                title="Net growth rate [1/h]",
                title_side="right",
                tickfont=dict(size=11),
                y=0.5,
                len=1,
                thickness=12,
            ),
            showscale=False,
            contours=dict(showlines=False, start=-0.8, end=0.8, size=0.2),
        )
    )
    fig.update_layout(
        xaxis=dict(
            title="µ [1/h]",
            ticks="inside",
        ),
        yaxis=dict(
            title="Cefotaxime [μg/mL]",
            ticks="inside",
        ),
        width=200,
        height=200,
        showlegend=False,
    )
    fig = style_plot(fig, line_thickness=1, left_margin=30, font_size=11)
    fig.write_image("plots/isoclines/cefotaxime.svg")


def net_growth_chloramphenicol():
    e = Experiment()
    fig = go.Figure()
    J = np.linspace(0, ST_20_CFA, 100)
    chloram = np.linspace(0, e.M_chloram, 100)
    J_grid, chloram_grid = np.meshgrid(np.array(J), np.array(chloram))
    JD = J_grid * 1 / (1 + chloram_grid / e.S.KT)
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            x=J,
            y=chloram,
            z=JD,
            colorscale="Cividis",
            colorbar=dict(
                title="Net growth rate [1/h]",
                title_side="right",
                tickfont=dict(size=11),
                y=0.5,
                len=1,
                thickness=12,
            ),
            contours=dict(showlines=False),
            showscale=False,
        )
    )
    fig.update_layout(
        xaxis=dict(
            title="ST growth rate [1/h]",
            ticks="inside",
        ),
        yaxis=dict(
            title="Chloramphenicol [μg/mL]",
            ticks="inside",
        ),
        width=200,
        height=200,
    )
    fig = style_plot(fig, line_thickness=1, left_margin=30, font_size=11)
    fig.write_image("plots/isoclines/chloramphenicol.svg")
