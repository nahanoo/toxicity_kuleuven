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
    # Setting antibitoics inflow to 0
    e = Experiment()
    e.E.r = 1.8
    e.S.r = 1.5
    e.total_transfers = 4
    e.M_cefo = 0
    e.M_chloram = 0
    e.transfer()
    e.plot_N("no_antibiotics/e_wins")

    e = Experiment()
    e.E.r = 1.5
    e.S.r = 1.8
    e.total_transfers = 4
    e.M_cefo = 0
    e.M_chloram = 0
    e.transfer()
    e.plot_N("no_antibiotics/s_wins")


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
    pass


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


two_sided_protection()
