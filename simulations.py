from model import Experiment
import numpy as np
import plotly.graph_objects as go
from style import *
from parameters import *


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
    r_e = EC_100_CFA
    r_s = ST_100_CFA
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
    r_e = EC_20_CFA
    r_s = ST_20_CFA
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


def cefotaxime_degraddation():
    def transfer(r_e, r_s, s_q):
        e = Experiment()
        e.E.r = r_e
        e.S.r = r_s
        e.S.q = s_q
        e.E.N0 = 0
        e.total_transfers = 1
        e.M_chloram = 0
        e.transfer()
        return e.plot_cefo(), e.plot_N()

    r_e = EC_100_CFA
    r_s = ST_100_CFA
    s_q = ST_Q
    fig, _ = transfer(r_e, r_s, s_q)
    fig.update_layout(
        title="100% CFA, r_EC: {:.3f}, r_ST: {:.3f}<br>q_ST: {:.1e}, a: 1".format(
            r_e,
            r_s,
            s_q,
        ),
    )
    fig = style_plot(fig, top_margin=40)
    fig.write_image("plots/transfers/s_protects_e/cefotaxime_degradation_100CFA.svg")

    r_e = EC_20_CFA
    r_s = ST_20_CFA
    s_q = ST_Q
    fig, _ = transfer(r_e, r_s, s_q)
    fig.update_layout(
        title="20% CFA, r_EC: {:.3f}, r_ST: {:.3f}<br>q_ST: {:.1e}, a: 1".format(
            r_e,
            r_s,
            s_q,
        ),
    )
    fig = style_plot(fig, top_margin=40)
    fig.write_image("plots/transfers/s_protects_e/cefotaxime_degradation_20CFA.svg")


def chloramphenicol_degradation():
    def transfer(r_e, r_s, e_q):
        e = Experiment()
        e.E.r = r_e
        e.S.r = r_s
        e.E.q = e_q
        e.S.N0 = 0
        e.total_transfers = 1
        e.M_cefo = 0
        e.transfer()
        return e.plot_chloram(), e.plot_N()

    r_e = EC_100_CFA
    r_s = ST_100_CFA
    e_q = EC_Q
    fig, _ = transfer(r_e, r_s, e_q)
    fig.update_layout(
        title="100% CFA, r_EC: {:.3f}, r_ST: {:.3f}<br>q_EC: {:.1e}, a: 1".format(
            r_e,
            r_s,
            e_q,
        ),
    )
    fig = style_plot(fig, top_margin=40)
    fig.write_image(
        "plots/transfers/e_protects_s/chloramphenicol_degradation_100CFA.svg"
    )

    r_e = EC_20_CFA
    r_s = ST_20_CFA
    e_q = EC_Q
    fig, _ = transfer(r_e, r_s, e_q)
    fig.update_layout(
        title="20% CFA, r_EC: {:.3f}, r_ST: {:.3f}<br>q_EC: {:.1e}, a: 1".format(
            r_e,
            r_s,
            e_q,
        ),
    )
    fig = style_plot(fig, top_margin=40)
    fig.write_image(
        "plots/transfers/e_protects_s/chloramphenicol_degradation_20CFA.svg"
    )


def cefotaxime_susceptibility():
    def transfer(r_e, r_s, e_u):
        e = Experiment()
        e.E.r = r_e
        e.S.r = r_s
        e.E.u = e_u
        e.S.N0 = 0
        e.total_transfers = 1
        e.M_chloram = 0
        e.transfer()
        return e.plot_cefo(), e.plot_N()

    r_e = EC_100_CFA
    r_s = ST_100_CFA
    e_u = EC_U
    IC50 = Experiment().E.KT

    _, fig = transfer(r_e, r_s, e_u)
    fig.update_layout(
        title="100% CFA, r_EC: {:.3f}, r_ST: {:.3f}<br>u: {:.1f},IC50: {:.3f}".format(
            r_e,
            r_s,
            e_u,
            IC50,
        ),
    )
    fig = style_plot(fig, top_margin=40)
    fig.write_image("plots/transfers/s_protects_e/EcN_death_100CFA.svg")

    r_e = EC_20_CFA
    r_s = ST_20_CFA
    e_u = EC_U
    _, fig = transfer(r_e, r_s, e_u)
    fig.update_layout(
        title="20% CFA, r_EC: {:.3f}, r_ST: {:.3f}<br>u: {:.1f},IC50: {:.3f}".format(
            r_e,
            r_s,
            e_u,
            IC50,
        ),
    )
    fig = style_plot(fig, top_margin=40)
    fig.write_image("plots/transfers/s_protects_e/EcN_death_20CFA.svg")


def chloramphenicol_susceptibility():
    def transfer(
        r_e,
        r_s,
    ):
        e = Experiment()
        e.E.r = r_e
        e.S.r = r_s
        e.E.N0 = 0
        e.total_transfers = 1
        e.M_cefo = 0
        e.transfer()
        return e.plot_chloram(), e.plot_N()

    r_e = EC_100_CFA
    r_s = ST_100_CFA
    IC50 = Experiment().S.KT
    _, fig = transfer(r_e, r_s)
    fig.update_layout(
        title="100% CFA, r_EC: {:.3f}, r_ST: {:.3f}<br>IC50: {:.3f}".format(
            r_e,
            r_s,
            IC50,
        ),
    )
    fig = style_plot(fig, top_margin=40)
    fig.write_image("plots/transfers/e_protects_s/ST_death_100CFA.svg")

    r_e = EC_20_CFA
    r_s = ST_20_CFA
    e_u = EC_U
    _, fig = transfer(r_e, r_s)
    fig.update_layout(
        title="20% CFA, r_EC: {:.3f}, r_ST: {:.3f}<br>IC50: {:.3f}".format(
            r_e,
            r_s,
            IC50,
        ),
    )
    fig = style_plot(fig, top_margin=40)
    fig.write_image("plots/transfers/e_protects_s/ST_death_20CFA.svg")


def s_protects_e():
    e = Experiment()
    e.E.r = EC_100_CFA
    e.S.r = ST_100_CFA
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
