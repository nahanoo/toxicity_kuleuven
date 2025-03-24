#!/usr/bin/env python
#
# Reserve 1 CPUs for this job
#
#SBATCH --cpus-per-task=16
#SBATCH --mem=1G
#
# Request it to run this for DD:HH:MM with ?G per core
#
#SBATCH --time=00:10:00
#
import sys

sys.path.append("/users/eulrich/toxicity_kuleuven")
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


def plot_competition(rE, rS, zs,title):
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            x=rE,
            y=rS,
            z=zs,
            colorscale=coolwarm_colorscale,
            colorbar=dict(
                title="Relative abundance <i>E. coli</i>",
                title_side="right",
                dtick=0.1,
                tickfont=dict(size=font_size),
            ),
            contours=dict(showlines=False, start=0, end=1, size=0.1),
        )
    )
    fig.update_layout(
        xaxis=dict(
            range=[0, max(rE)],
            dtick=0.2,
            title="Maximum growth rate <i>E. coli</i> [1/h]",
        ),
        yaxis=dict(
            range=[0, max(rS)],
            dtick=0.2,
            title="Maximum growth rate <i>Salmonella</i> [1/h]",
        ),
        title="<i>Salmonella</i> protects <i>E. coli</i>",
        width=width,
        height=height,
    )
    fig = style_plot(fig, line_thickness=1, left_margin=50)
    return fig


def competition_simulations():
    def no_antibiotics(rE, rS):
        e = Experiment()
        e.M_cefo = 0
        e.M_chloram = 0
        e.E.r = rE
        e.S.r = rS
        e.total_transfers = 50
        e.transfer()
        return e.E.N[-1] / (e.E.N[-1] + e.S.N[-1])

    def s_protects_e(rE, rS):
        e = Experiment()
        e.M_chloram = 0
        e.E.r = rE
        e.S.r = rS
        e.total_transfers = 50
        e.transfer()
        return e.E.N[-1] / (e.E.N[-1] + e.S.N[-1])

    def e_protects_s(rE, rS):
        e = Experiment()
        e.M_cefo = 0
        e.E.r = rE
        e.S.r = rS
        e.total_transfers = 50
        e.transfer()
        return e.E.N[-1] / (e.E.N[-1] + e.S.N[-1])

    def two_sided_protection(rE, rS):
        e = Experiment()
        e.E.r = rE
        e.S.r = rS
        e.total_transfers = 50
        e.transfer()
        return e.E.N[-1] / (e.E.N[-1] + e.S.N[-1])

    rs = np.linspace(0, 1.5, 100)
    zs = np.zeros((len(rs), len(rs)))
    results = Parallel(n_jobs=-1)(
        delayed(no_antibiotics)(rE, rS) for rE in rs for rS in rs
    )
    zs = np.array(results).reshape(len(rs), len(rs))
    fig = plot_competition(rs, rs, zs,"No antibiotics")
    fig.write_image("plots/coexistence/no_antibiotics/growth_rates.svg")

    zs = np.zeros((len(rs), len(rs)))
    results = Parallel(n_jobs=-1)(
        delayed(s_protects_e)(rE, rS) for rE in rs for rS in rs
    )
    zs = np.array(results).reshape(len(rs), len(rs))
    fig = plot_competition(rs, rs, zs,"<i>Salmonella</i> protects <i>E. coli</i>")
    fig.write_image("plots/coexistence/s_protects_e/growth_rates.svg")

    zs = np.array(results).reshape(len(rs), len(rs))
    results = Parallel(n_jobs=-1)(
        delayed(e_protects_s)(rE, rS) for rE in rs for rS in rs
    )
    zs = np.array(results).reshape(len(rs), len(rs))
    fig = plot_competition(rs, rs, zs,"<i>E. coli</i> protects <i>Salmonella</i>")
    fig.write_image("plots/coexistence/e_protects_s/growth_rates.svg")

    zs = np.array(results).reshape(len(rs), len(rs))
    results = Parallel(n_jobs=-1)(
        delayed(two_sided_protection)(rE, rS) for rE in rs for rS in rs
    )
    zs = np.array(results).reshape(len(rs), len(rs))
    fig = plot_competition(rs, rs, zs,"Two sided protection")
    fig.write_image("plots/coexistence/two_sided/growth_rates.svg")


competition_simulations()