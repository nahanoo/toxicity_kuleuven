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
        self.a = 1
        # detoxification factor
        self.q = 80
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
        # detoxification factor
        self.q = 1
        # cefo death rate
        self.u = 1
        # population density
        self.N = np.array([0.1])


# eq.e_min_death_rate()
class Experiment:

    def __init__(self):
        # transfers
        self.total_transfers = 10
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

    def transfer(self, model, y0):
        for i in range(self.total_transfers):
            steps = 100
            t = np.linspace(0, self.transfer_period, steps)
            Y = odeint(model, y0, t)
            self.E.N = np.concatenate((self.E.N, Y[:, 0]))
            self.S.N = np.concatenate((self.S.N, Y[:, 1]))
            if len(Y[-1]) == 3:
                self.cefo = np.concatenate((self.cefo, Y[:, 2]))
            if len(Y[-1]) == 4:
                self.chloram = np.concatenate((self.chloram, Y[:, 3]))
            self.time = np.concatenate((self.time, self.time[-1] + t), axis=0)

            self.E.N[-1] = self.E.N[-1] / self.dilution_factor
            self.S.N[-1] = self.S.N[-1] / self.dilution_factor
            y0[0], y0[1] = self.E.N[-1], self.S.N[-1]

    def no_antibiotics(self, y, t):
        E, S = y
        dE = self.E.r * (1 - (E + S) / self.E.K) * E
        dS = self.S.r * (1 - (S + E) / self.S.K) * S
        return dE, dS

    def s_e_chloram(self, y, t):
        E, S = y
        dE = self.E.r * (1 - (E + S) / self.E.K) * E
        dS = (
            self.S.r * (1 - (E + S) / self.S.K) / (1 + (self.M_chloram / self.S.KT)) * S
        )
        return dE, dS

    def s_e_cefotaxime(self, y, t):
        E, S = y
        dE = (
            self.E.r * (1 - (E + S) / self.E.K) * E
            - self.E.u * self.M_cefo / (self.M_cefo + self.E.KT) * E
        )
        dS = self.S.r * (1 - (E + S) / self.S.K) * S
        return dE, dS

    def s_e_chloram_cefo(self, y, t):
        E, S = y
        dE = (
            self.E.r * (1 - (E + S) / self.E.K) * E
            - self.E.u * self.M_cefo / (self.M_cefo + self.E.KT) * E
        )
        dS = (
            self.S.r * (1 - (E + S) / self.S.K) / (1 + (self.M_chloram / self.S.KT)) * S
        )
        return dE, dS

    def s_e_ceftaxime_degradation(self, y, t):
        E, S, cefo = y
        dE = (
            self.E.r * (1 - (E + S) / self.E.K) * E
            - self.E.u * cefo / (cefo + self.E.KT) * E
        )
        dS = self.S.r * (1 - (E + S) / self.S.K) * S
        dCefo = -self.S.a * S / self.S.q
        return dE, dS, dCefo

    def transfer_no_antibiotics(self):
        self.transfer(self.no_antibiotics, [self.E.N[0], self.S.N[0]])
        self.plot_N("no_antibiotics")

    def transfer_s_e_cefotaxime(self):
        self.transfer(self.s_e_cefotaxime, [self.E.N[0], self.S.N[0]])
        self.plot_N("s_e_cefotaxime")

    def transfer_s_e_chloram(self):
        self.transfer(self.s_e_chloram, [self.E.N[0], self.S.N[0]])
        self.plot_N("s_e_chloram")

    def transfer_s_e_chloram_cefo(self):
        self.transfer(self.s_e_chloram_cefo, [self.E.N[0], self.S.N[0]])
        self.plot_N("s_e_chloram_cefo")

    def transfer_s_e_cefotaxime_degradation(self):
        self.transfer(
            self.s_e_ceftaxime_degradation, [self.E.N[0], self.S.N[0], self.M_cefo]
        )
        self.plot_N("s_e_cefotaxime_degradation")
        self.plot_cefo("s_cefotaxime_degradation")

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

    def plot_cefo(self, fname):
        fig = go.Figure(
            go.Scatter(
                x=self.time,
                y=self.cefo,
                mode="lines",
                name="Cefotaxime",
                line=dict(color="black"),
            )
        )

        fig.update_layout(
            xaxis=dict(
                range=[0, max(self.time)],
                dtick=12,
                title="Time [1/h]",
            ),
            yaxis=dict(
                range=[0, 0.3],
                dtick=0.05,
                title="Cefotaxime concentration [µg/mL]",
            ),
            width=width,
            height=height,
        )
        fig = style_plot(fig, line_thickness=1.7)
        fig.write_image("plots/transfers/" + fname + ".svg")


e = Experiment()
e.transfer_no_antibiotics()

e = Experiment()
e.transfer_s_e_cefotaxime()

e = Experiment()
e.transfer_s_e_chloram()

e = Experiment()
e.transfer_s_e_chloram_cefo()

e = Experiment()
e.transfer_s_e_chloram_cefo()

e = Experiment()
e.transfer_s_e_cefotaxime_degradation()
