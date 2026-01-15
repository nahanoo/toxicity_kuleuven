#!/usr/bin/env python
#
# Reserve 1 CPUs for this job
#
# SBATCH --cpus-per-task=16
# SBATCH --mem=4G
#
# Request it to run this for DD:HH:MM with ?G per core
#
# SBATCH --time=00:30:00
#
import sys

sys.path.append("/users/eulrich/toxicity_kuleuven")
from model import Experiment
import numpy as np
import plotly.graph_objects as go
from style import *

coolwarm_colorscale = [
    [0.0, "#3B4CC0"],
    [0.25, "#8CABD8"],
    [0.5, "#FFFFFF"],
    [0.75, "#D9907E"],
    [1.0, "#B40426"],
]


def no_antibiotics():
    # Setting antibiotics inflow to 0
    def transfer(r_e, r_s):
        e = Experiment()
        e.E.r = r_e
        e.S.r = r_s
        e.total_transfers = 4
        e.M_cefo = 0
        e.M_chloram = 0
        e.transfer()
        return e.plot_N()

    # 100% CFA
    r_e = 1.048
    r_s = 0.87
    fig = transfer(r_e, r_s)
    fig.update_layout(
        title="No antibiotics: 100% CFA<br>r_EC: {:.3f}, r_ST: {:.3f}".format(r_e, r_s),
        showlegend=True,
        legend=dict(
            x=0.4,
            y=0.1,
            bgcolor="rgba(0, 0, 0, 0)",
        ),
    )
    fig = style_plot(fig, top_margin=40)
    fig.write_image("plots/transfers/no_antibiotics/ecn_st_100_cfa.svg")

    # 20% CFA
    r_e = 0.831
    r_s = 0.915
    fig = transfer(r_e, r_s)
    fig.update_layout(
        title="No antibiotics: 20% CFA<br>r_EC: {:.3f}, r_ST: {:.3f}".format(r_e, r_s),
        showlegend=True,
        legend=dict(
            x=0.4,
            y=0.1,
            bgcolor="rgba(0, 0, 0, 0)",
        ),
    )
    fig = style_plot(fig, top_margin=40)
    fig.write_image("plots/transfers/no_antibiotics/ecn_st_20_cfa.svg")


def s_protects_e():
    # No degradation
    e = Experiment()
    e.total_transfers = 1
    e.M_chloram = 0
    e.S.a = 0
    e.S.r = 0
    e.S.N0 = 0
    e.transfer()
    e.plot_N("s_protects_e/abundance_no_degradation")
    e.plot_cefo("s_protects_e/cefotaxime_no_degradation")

    # S mono-culture degradation
    e = Experiment()
    e.total_transfers = 1
    e.M_chloram = 0
    e.E.N0 = 0
    e.transfer()
    e.plot_cefo("s_protects_e/cefotaxime_degradation")

    e = Experiment()
    e.E.r = 0.25
    e.S.r = 1.8
    e.E.u = 1
    # e.S.a = 0.5
    e.total_transfers = 4
    e.M_chloram = 0
    e.transfer()
    e.plot_N("s_protects_e/abundance")


def e_protects_s():
    # No degradation
    e = Experiment()
    e.total_transfers = 1
    e.M_cefo = 0
    e.E.a = 0
    e.E.N0 = 0
    e.transfer()
    e.plot_N("e_protects_s/abundance_no_degradation")
    e.plot_chloram("e_protects_s/chloram_no_degradation")

    # mono-culture degradation
    e = Experiment()
    e.total_transfers = 1
    e.M_cefo = 0
    e.S.N0 = 0
    e.transfer()
    e.plot_chloram("e_protects_s/chloram_degradation")

    e = Experiment()
    e.E.r = 0.25
    e.S.r = 1.9
    e.total_transfers = 8
    e.M_cefo = 0
    e.transfer()
    e.plot_N("e_protects_s/abundance")


def two_sided_protection():
    e = Experiment()
    e.total_transfers = 4
    # e.E.r = 1.8
    # e.S.r = 1.5
    # e.E.u = 1.8
    # e.S.a = 0.1
    e.E.r = 1.8
    e.S.r = 1.5
    e.E.u = 1.9
    e.S.a = 0.1
    e.transfer()
    e.plot_N("two_sided_protection/abundance")
    e.plot_cefo("two_sided_protection/cefotaxime")
    e.plot_chloram("two_sided_protection/chloram")
