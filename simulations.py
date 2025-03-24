from model import Experiment
import numpy as np
from scipy.integrate import odeint
import plotly.graph_objects as go
from style import *
from joblib import Parallel, delayed
from timeit import default_timer as timer


def no_antibiotics():
    # Setting antibitoics inflow to 0
    e = Experiment()
    e.M_cefo = 0
    e.M_chloram = 0
    e.transfer()
    e.plot_N("no_antibiotics/abundance")


def s_protects_e():
    # Only adding cefotaxime thus setting chloramphenicol inflow to 0
    e = Experiment()
    e.M_chloram = 0
    # Case with no degradation, setting degradation rate to 0
    e.S.a = 0
    e.transfer()
    e.plot_N("s_protects_e/abundance_no_degradation")
    e.plot_cefo("s_protects_e/cefotaxime_no_degradation")

    # Case considering degradation
    e = Experiment()
    e.M_chloram = 0
    e.transfer()
    e.plot_N("s_protects_e/abundance")
    e.plot_cefo("s_protects_e/cefotaxime")


def e_protects_s():
    # Only adding chloramphenicol thus setting cefotaxime inflow to 0
    e = Experiment()
    e.M_cefo = 0
    # Case with no degradation, setting degradation rate to 0
    e.E.a = 0
    e.transfer()
    e.plot_N("e_protects_s/abundance_no_degradation")
    e.plot_chloram("e_protects_s/chloram_no_degradation")

    e = Experiment()
    e.M_cefo = 0
    # Case with degradation
    e.transfer()
    e.plot_N("e_protects_s/abundance")
    e.plot_chloram("e_protects_s/chloram")


def two_sided_protection():
    e = Experiment()
    # Case with no degradation, setting degradation rate to 0
    e.E.a = 0
    e.S.a = 0
    e.transfer()
    e.plot_N("two_sided_protection/abundance_no_degradation")

    e = Experiment()
    e.transfer()
    e.plot_N("two_sided_protection/abundance")
    e.plot_cefo("two_sided_protection/cefotaxime")
    e.plot_chloram("two_sided_protection/chloram")


def run_experiment(rE, rS):
    e = Experiment()
    e.E.r = rE
    e.S.r = rS
    e.total_transfers = 50
    e.transfer()
    return e.E.N[-1] / (e.E.N[-1] + e.S.N[-1])


def co_existence_s_protects_e():
    start = timer()
    rs = np.linspace(0, 1.5, 100)
    zs = np.zeros((len(rs), len(rs)))

    # Parallel computation
    results = Parallel(n_jobs=-1)(
        delayed(run_experiment)(rE, rS) for rE in rs for rS in rs
    )

    # Reshape back to matrix form
    zs = np.array(results).reshape(len(rs), len(rs))
    stop = timer()
    print(stop - start)

    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            x=rs,
            y=rs,
            z=zs,
            colorscale=coolwarm_colorscale,
            colorbar=dict(
                title="Relative abundance <i>E. coli</i>",
                title_side="right",
                dtick=0.2,
                tickfont=dict(size=font_size),
            ),
            contours=dict(showlines=False, start=0, end=1, size=0.2),
        )
    )
    fig.update_layout(
        xaxis=dict(
            range=[0, max(rs)],
            dtick=0.2,
            title="Maximum growth rate <i>E. coli</i> [1/h]",
        ),
        yaxis=dict(
            range=[0, max(rs)],
            dtick=0.2,
            title="Maximum growth rate <i>Salmonella</i> [1/h]",
        ),
        title="<i>Salmonella</i> protects <i>E. coli</i>",
        width=width,
        height=height,
    )
    fig = style_plot(fig, line_thickness=1, left_margin=50)
    fig.write_image("plots/coexistence/s_protects_e/growth_rates.svg")
