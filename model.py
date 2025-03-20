import numpy as np
from scipy.integrate import odeint
import plotly.graph_objects as go
from style import *


class Salmonella:
    def __init__(self):
        # max. growth rate
        self.r = 0.7
        # monod constant
        self.K = 1.1
        # toxin monod constance
        self.KT = 0.32
        # toxin degradation rate
        self.a = 50
        # population density
        self.N = np.array([0.1])
        # passive uptake rate


class E_coli:
    def __init__(self):
        # max. growth rate
        self.r = 0.9
        # monod constant
        self.K = 1.1
        # toxin monod constance
        self.KT = 0.03
        # toxin degradation rate
        self.a = 1
        # cefo death rate
        self.u = 1
        # population density
        self.N = np.array([0.1])


# eq.e_min_death_rate()
class Experiment:

    def __init__(self):
        # transfers
        self.total_transfers = 4
        # dilution factor
        self.dilution_factor = 100
        # transfer period
        self.transfer_period = 24
        # cefotaxime concentration
        self.M_cefo = 0.25
        # chloramphenicol concentration
        self.M_chloram = 16
        # array for cefotaxime concentration
        self.cefo = np.array([self.M_cefo])
        # array for chlorampheicol concentration
        self.chloram = np.array([self.M_chloram])
        # array for time
        self.time = np.array([0])
        # create bugs
        self.S = Salmonella()
        self.E = E_coli()

    def transfer(self, model):
        for i in range(self.total_transfers):
            steps = 1000
            t = np.linspace(0, self.transfer_period, steps)
            Y = odeint(model, [self.E.N[-1], self.S.N[-1]], t)
            self.E.N = np.concatenate((self.E.N, Y[:, 0]))
            self.S.N = np.concatenate((self.S.N, Y[:, 1]))
            self.time = np.concatenate((self.time, self.time[-1] + t), axis=0)

            self.E.N[-1] = self.E.N[-1] / self.dilution_factor
            self.S.N[-1] = self.S.N[-1] / self.dilution_factor

    def no_antibiotics(self, y, t):
        E, S = y
        dE = self.E.r * (1 - (E + S) / self.E.K) * E
        dS = self.S.r * (1 - (S + E) / self.S.K) * S
        return dE, dS

    def transfer_no_antibiotics(self):
        self.transfer(self.no_antibiotics)

    def s_e_cefotaxime(self, y, t):
        E, S = y
        dE = (
            self.E.r * (1 - (E + S) / self.E.K) * E
            - self.E.u * self.M_cefo / (self.M_cefo + self.E.KT) * E
        )
        dS = self.S.r * (1 - (E + S) / self.S.K) * S
        return dE, dS

    def transfer_s_e_cefotaxime(self):
        self.transfer(self.s_e_cefotaxime)

    def s_e_chloram(self, y, t):
        E, S = y
        dE = self.E.r * (1 - (E + S) / self.E.K) * E
        dS = (
            self.S.r
            * (1 - (E + S) / self.S.K)
            * 1
            / (1 + (self.M_chloram / self.S.KT))
            * S
        )
        return dE, dS

    def transfer_s_e_chloram(self):
        self.transfer(self.s_e_chloram)

    def plot_N(self, fname):
        fig = go.Figure(
            go.Scatter(
                x=self.time,
                y=self.E.N,
                mode="lines",
                name="E_coli",
                line=dict(color=colors["ecoli"]),
            )
        )
        fig.add_trace(
            go.Scatter(
                x=self.time,
                y=self.S.N,
                mode="lines",
                name="Salmonella",
                line=dict(color=colors["st"]),
            )
        )
        fig.update_layout(
            xaxis=dict(
                range=[0, max(self.time)],
                dtick=12,
                title="Time [1/h]",
            ),
            yaxis=dict(
                range=[0, 1.2],
                dtick=0.2,
                title="Abundance",
            ),
            width=width,
            height=height,
        )
        fig = style_plot(fig, line_thickness=1.7)
        fig.write_image("plots/transfers/" + fname + ".svg")


e = Experiment()
e.transfer_no_antibiotics()
e.plot_N("no_antibiotics")

e = Experiment()
e.transfer_s_e_cefotaxime()
e.plot_N("s_e_cefotaxime")

e = Experiment()
e.transfer_s_e_chloram()
e.plot_N("s_e_chloram")
