import plotly.graph_objects as go
import numpy as np
from style import *
from model import Experiment
from joblib import Parallel, delayed
from parameters import *


def no_antibiotics():
    def transfer(r_e, r_s):
        e = Experiment()  # construct in worker
        e.E.r, e.S.r = r_e, r_s
        e.M_cefo = 0
        e.M_chloram = 0
        e.total_transfers = 4
        e.transfer()
        E = e.E.N[-1]
        S = e.S.N[-1]
        return E, S

    rs_E = np.linspace(0, 1, 20)
    rs_S = np.linspace(0, 1, 20)

    results = Parallel(
        n_jobs=-1, backend="loky", prefer="processes", batch_size="auto"
    )(delayed(transfer)(r_e, r_s) for r_e in rs_E for r_s in rs_S)
    E, S = map(np.array, zip(*results))

    # reshape to grids
    n = len(rs_E)
    E = E.astype(float).reshape(n, n).T
    S = S.astype(float).reshape(n, n).T
    thr = 1e5
    alive_E = E > thr
    alive_S = S > thr

    both_alive = alive_E & alive_S
    E_excludes = alive_E & (~alive_S)
    S_excludes = (~alive_E) & alive_S

    z = np.full_like(E, np.nan, dtype=float)
    z[both_alive] = np.log10(E[both_alive] / S[both_alive])

    # pick two colors (you can change these)
    c_E_excludes = "#E76F51"  # dark (cividis-ish) for "E wins"
    c_S_excludes = "#2A9D8F"  # bright (cividis-ish) for "S wins"

    fig = go.Figure()

    # 1) Continuous (Cividis) where both alive
    finite = np.isfinite(z)
    zmin = float(np.nanmin(z)) if finite.any() else -1
    zmax = float(np.nanmax(z)) if finite.any() else 1

    fig.add_trace(
        go.Heatmap(
            x=rs_E,
            y=rs_S,
            z=z,
            colorscale="Cividis",
            zmin=zmin,
            zmax=zmax,
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
            x=rs_E,
            y=rs_S,
            z=np.where(E_excludes, 1.0, np.nan),
            colorscale=[[0, c_E_excludes], [1, c_E_excludes]],
            showscale=False,
        )
    )

    fig.add_trace(
        go.Heatmap(
            x=rs_E,
            y=rs_S,
            z=np.where(S_excludes, 1.0, np.nan),
            colorscale=[[0, c_S_excludes], [1, c_S_excludes]],
            showscale=False,
        )
    )
    fig.update_layout(
        xaxis=dict(title="Max. growth rate EcN (1/h)", ticks="inside", showgrid=False),
        yaxis=dict(title="Max. growth rate ST (1/h)", ticks="inside", showgrid=False),
        title="No antibiotics: final ratio EcN/ST",
        width=220,
        height=180,
    )
    fig = style_plot(fig, line_thickness=1.5, font_size=11)
    fig.write_image("plots/contours/no_antibiotics/sweep_growth_rates.svg")


def ceftoaxime_growth_rate_sweep():
    def transfer(r_e, r_s):
        e = Experiment()  # construct in worker
        e.E.r, e.S.r = r_e, r_s
        e.M_chloram = 0
        e.total_transfers = 4
        e.transfer()
        E = e.E.N[-1]
        S = e.S.N[-1]
        return E, S

    rs_E = np.linspace(0, 1.2, 20)  # e_r (should go on y)
    rs_S = np.linspace(0, 1.2, 20)  # s_r (should go on x)

    ne = len(rs_E)
    ns = len(rs_S)

    results = Parallel(
        n_jobs=-1, backend="loky", prefer="processes", batch_size="auto"
    )(delayed(transfer)(r_e, r_s) for r_e in rs_E for r_s in rs_S)
    E, S = map(np.array, zip(*results))

    # rows: r_e (y), cols: r_s (x)
    E = E.astype(float).reshape(ne, ns)
    S = S.astype(float).reshape(ne, ns)
    thr = 1e5
    alive_E = E > thr
    alive_S = S > thr

    both_alive = alive_E & alive_S
    E_excludes = alive_E & (~alive_S)
    S_excludes = (~alive_E) & alive_S

    z = np.full_like(E, np.nan, dtype=float)
    z[both_alive] = np.log10(E[both_alive] / S[both_alive])

    c_E_excludes = "#E76F51"  # dark (cividis-ish) for "E wins"
    c_S_excludes = "#2A9D8F"  # bright (cividis-ish) for "S wins"

    fig = go.Figure()

    finite = np.isfinite(z)
    zmin = float(np.nanmin(z)) if finite.any() else -1
    zmax = float(np.nanmax(z)) if finite.any() else 1

    fig.add_trace(
        go.Heatmap(
            x=rs_S,
            y=rs_E,
            z=z,
            colorscale="Cividis",
            zmin=zmin,
            zmax=zmax,
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
            x=rs_S,
            y=rs_E,
            z=np.where(E_excludes, 1.0, np.nan),
            colorscale=[[0, c_E_excludes], [1, c_E_excludes]],
            showscale=False,
            hoverinfo="skip",
        )
    )

    fig.add_trace(
        go.Heatmap(
            x=rs_S,
            y=rs_E,
            z=np.where(S_excludes, 1.0, np.nan),
            colorscale=[[0, c_S_excludes], [1, c_S_excludes]],
            showscale=False,
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[ST_100_CFA, ST_20_CFA],
            y=[EC_100_CFA, EC_20_CFA],
            mode="markers+text",
            marker=dict(color="red", size=10),
            text=["100% CFA", "20% CFA"],
            textfont=dict(color="white", size=11),
            textposition="top center",
            showlegend=False,
        )
    )
    fig.update_layout(
        xaxis=dict(title="Max. growth rate ST (1/h)", ticks="inside", showgrid=False),
        yaxis=dict(title="Max. growth rate EcN (1/h)", ticks="inside", showgrid=False),
        title="Cf",
        width=220,
        height=180,
    )
    fig = style_plot(fig, line_thickness=1.5, font_size=11, marker_size=6)
    fig.write_image("plots/contours/s_protects_e/sweep_growth_rates.svg")


def cefotaxime_degradation_rate_sweep(r_e, r_s, fname):
    def transfer(r_e, r_s, s_a, e_u):
        e = Experiment()
        e.S.a = s_a
        e.E.u = e_u
        e.E.r = r_e
        e.S.r = r_s
        e.M_chloram = 0
        e.total_transfers = 4
        e.transfer()
        return e.E.N[-1], e.S.N[-1]

    s_a = np.linspace(0, 1, 20)
    e_u = np.linspace(0, 2, 20)

    results = Parallel(
        n_jobs=-1, backend="loky", prefer="processes", batch_size="auto"
    )(
        delayed(transfer)(r_e, r_s, s_a_val, e_u_val)
        for s_a_val in s_a
        for e_u_val in e_u
    )
    E, S = map(np.array, zip(*results))

    ns = len(s_a)
    nu = len(e_u)

    E = E.astype(float).reshape(ns, nu)  # rows: s_a (y), cols: e_u (x)
    S = S.astype(float).reshape(ns, nu)

    thr = 1e5
    alive_E = E > thr
    alive_S = S > thr

    both_alive = alive_E & alive_S
    E_excludes = alive_E & (
        ~alive_S
    )  # EcN survives, ST extinct? (depends on who is E/S)
    S_excludes = (~alive_E) & alive_S

    z = np.full_like(E, np.nan, dtype=float)
    z[both_alive] = np.log10(E[both_alive] / S[both_alive])

    c_E_excludes = "#E76F51"  # coral
    c_S_excludes = "#2A9D8F"  # teal

    finite = np.isfinite(z)
    zmin = float(np.nanmin(z)) if finite.any() else -1
    zmax = float(np.nanmax(z)) if finite.any() else 1

    fig = go.Figure()

    fig.add_trace(
        go.Heatmap(
            x=e_u,
            y=s_a,
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
            x=e_u,
            y=s_a,
            z=np.where(E_excludes, 1.0, np.nan),
            colorscale=[[0, c_E_excludes], [1, c_E_excludes]],
            showscale=False,
            hoverongaps=False,
        )
    )

    fig.add_trace(
        go.Heatmap(
            x=e_u,
            y=s_a,
            z=np.where(S_excludes, 1.0, np.nan),
            colorscale=[[0, c_S_excludes], [1, c_S_excludes]],
            showscale=False,
            hoverongaps=False,
        )
    )

    fig.update_layout(
        xaxis=dict(title="Death rate EcN (1/h)", ticks="inside", showgrid=False),
        yaxis=dict(
            title="Max. degradation rate ST (1/h)", ticks="inside", showgrid=False
        ),
        title="Cf",
        width=220,
        height=180,
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    fig = style_plot(fig, line_thickness=1.5, font_size=11)
    fig.write_image(f"plots/contours/s_protects_e/{fname}")


def plot_ceftoaxime_growth_rate_sweep():
    cefotaxime_degradation_rate_sweep(
        EC_100_CFA, ST_100_CFA, fname="sweep_degradation_rate_100CFA.svg"
    )
    cefotaxime_degradation_rate_sweep(
        EC_20_CFA, ST_20_CFA, fname="sweep_degradation_rate_20CFA.svg"
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


def chloramphenicol_growth_rate_sweep():
    def transfer(r_e, r_s):
        e = Experiment()  # construct in worker
        e.E.r, e.S.r = r_e, r_s
        e.M_cefo = 0
        e.total_transfers = 4
        e.transfer()
        E = e.E.N[-1]
        S = e.S.N[-1]
        return E, S

    rs_E = np.linspace(0, 1.2, 20)  # e_r (should go on y)
    rs_S = np.linspace(0, 1.2, 20)  # s_r (should go on x)

    ne = len(rs_E)
    ns = len(rs_S)

    results = Parallel(
        n_jobs=-1, backend="loky", prefer="processes", batch_size="auto"
    )(delayed(transfer)(r_e, r_s) for r_e in rs_E for r_s in rs_S)
    E, S = map(np.array, zip(*results))

    # rows: r_e (y), cols: r_s (x)
    E = E.astype(float).reshape(ne, ns)
    S = S.astype(float).reshape(ne, ns)
    thr = 1e5
    alive_E = E > thr
    alive_S = S > thr

    both_alive = alive_E & alive_S
    E_excludes = alive_E & (~alive_S)
    S_excludes = (~alive_E) & alive_S

    z = np.full_like(E, np.nan, dtype=float)
    z[both_alive] = np.log10(E[both_alive] / S[both_alive])

    c_E_excludes = "#E76F51"  # dark (cividis-ish) for "E wins"
    c_S_excludes = "#2A9D8F"  # bright (cividis-ish) for "S wins"

    fig = go.Figure()

    finite = np.isfinite(z)
    zmin = float(np.nanmin(z)) if finite.any() else -1
    zmax = float(np.nanmax(z)) if finite.any() else 1

    fig.add_trace(
        go.Heatmap(
            x=rs_S,
            y=rs_E,
            z=z,
            colorscale="Cividis",
            zmin=zmin,
            zmax=zmax,
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
            x=rs_S,
            y=rs_E,
            z=np.where(E_excludes, 1.0, np.nan),
            colorscale=[[0, c_E_excludes], [1, c_E_excludes]],
            showscale=False,
            hoverinfo="skip",
        )
    )

    fig.add_trace(
        go.Heatmap(
            x=rs_S,
            y=rs_E,
            z=np.where(S_excludes, 1.0, np.nan),
            colorscale=[[0, c_S_excludes], [1, c_S_excludes]],
            showscale=False,
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[ST_100_CFA, ST_20_CFA],
            y=[EC_100_CFA, EC_20_CFA],
            mode="markers+text",
            marker=dict(color="red", size=10),
            text=["100% CFA", "20% CFA"],
            textfont=dict(color="white", size=11),
            textposition="top center",
            showlegend=False,
        )
    )
    fig.update_layout(
        xaxis=dict(title="Max. growth rate ST (1/h)", ticks="inside", showgrid=False),
        yaxis=dict(title="Max. growth rate EcN (1/h)", ticks="inside", showgrid=False),
        title="Cm",
        width=220,
        height=180,
    )
    fig = style_plot(fig, line_thickness=1.5, font_size=11, marker_size=6)
    fig.write_image("plots/contours/e_protects_s/sweep_growth_rates.svg")


def cefotaxime_chloramphenicol_growth_rate_sweep():
    def transfer(r_e, r_s):
        e = Experiment()  # construct in worker
        e.E.r, e.S.r = r_e, r_s
        e.total_transfers = 4
        e.transfer()
        E = e.E.N[-1]
        S = e.S.N[-1]
        return E, S

    rs_E = np.linspace(0, 1.2, 20)  # e_r (should go on y)
    rs_S = np.linspace(0, 1.2, 20)  # s_r (should go on x)

    ne = len(rs_E)
    ns = len(rs_S)

    results = Parallel(
        n_jobs=-1, backend="loky", prefer="processes", batch_size="auto"
    )(delayed(transfer)(r_e, r_s) for r_e in rs_E for r_s in rs_S)
    E, S = map(np.array, zip(*results))

    # rows: r_e (y), cols: r_s (x)
    E = E.astype(float).reshape(ne, ns)
    S = S.astype(float).reshape(ne, ns)
    thr = 1e5
    alive_E = E > thr
    alive_S = S > thr

    both_alive = alive_E & alive_S
    E_excludes = alive_E & (~alive_S)
    S_excludes = (~alive_E) & alive_S

    z = np.full_like(E, np.nan, dtype=float)
    z[both_alive] = np.log10(E[both_alive] / S[both_alive])

    c_E_excludes = "#E76F51"  # dark (cividis-ish) for "E wins"
    c_S_excludes = "#2A9D8F"  # bright (cividis-ish) for "S wins"

    fig = go.Figure()

    finite = np.isfinite(z)
    zmin = float(np.nanmin(z)) if finite.any() else -1
    zmax = float(np.nanmax(z)) if finite.any() else 1

    fig.add_trace(
        go.Heatmap(
            x=rs_S,
            y=rs_E,
            z=z,
            colorscale="Cividis",
            zmin=zmin,
            zmax=zmax,
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
            x=rs_S,
            y=rs_E,
            z=np.where(E_excludes, 1.0, np.nan),
            colorscale=[[0, c_E_excludes], [1, c_E_excludes]],
            showscale=False,
            hoverinfo="skip",
        )
    )

    fig.add_trace(
        go.Heatmap(
            x=rs_S,
            y=rs_E,
            z=np.where(S_excludes, 1.0, np.nan),
            colorscale=[[0, c_S_excludes], [1, c_S_excludes]],
            showscale=False,
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[ST_100_CFA, ST_20_CFA],
            y=[EC_100_CFA, EC_20_CFA],
            mode="markers+text",
            marker=dict(color="red", size=10),
            text=["100% CFA", "20% CFA"],
            textfont=dict(color="white", size=11),
            textposition="top center",
            showlegend=False,
        )
    )
    fig.update_layout(
        xaxis=dict(title="Max. growth rate ST (1/h)", ticks="inside", showgrid=False),
        yaxis=dict(title="Max. growth rate EcN (1/h)", ticks="inside", showgrid=False),
        title="Cm",
        width=220,
        height=180,
    )
    fig = style_plot(fig, line_thickness=1.5, font_size=11, marker_size=6)
    fig.write_image("plots/contours/two_sided/sweep_growth_rates.svg")


cefotaxime_chloramphenicol_growth_rate_sweep()
