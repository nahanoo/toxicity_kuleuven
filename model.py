import numpy as np
from scipy.integrate import odeint
import plotly.graph_objects as go
from style import *


class Salmonella:
    def __init__(self):
        # max. growth rate
        self.r = 0.45
        # Carryin capacity
        self.K = 1
        # toxin monod constance
        self.KT = 0.32
        # toxin degradation rate
        self.a = 1
        # detoxification factor
        self.q = 0.05
        # Initial population density
        self.N0 = 0.1
        # Population density
        self.N = np.array([])
        # passive uptake rate


class E_coli:
    def __init__(self):
        # max. growth rate
        self.r = 0.9
        # Carrying capacity
        self.K = 1
        # toxin monod constance
        self.KT = 0.03
        # toxin degradation rate
        self.a = 1
        # detoxification factor
        self.q = 0.05
        # cefo death rate
        self.u = 1
        # Initial population density
        self.N0 = 0.1
        # population density
        self.N = np.array([])


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
        # Time after cefotaixme is added
        self.cefo_start = 2
        # chloramphenicol concentration
        self.M_chloram = 16
        # array for cefotaxime concentration
        self.cefo = np.array([0])
        # array for chlorampheicol concentration
        self.chloram = np.array([self.M_chloram])
        # array for time
        self.time = np.array([])
        # create bugs
        self.S = Salmonella()
        self.E = E_coli()

    def growth_e(self, E, S):
        # Growth function for e coli
        return self.E.r * (1 - (E + S) / self.E.K)

    def growth_s(self, E, S):
        # Growth function for salmonella
        return self.S.r * (1 - (E + S) / self.S.K)

    def death_e(self, cefo):
        # Death function for e coli
        return self.E.u * cefo / (cefo + self.E.KT)

    def death_s(self, chloram):
        # Death function for salmonella
        return 1 / (1 + (chloram / self.S.KT))

    def degradation_e(self, E):
        # Degradation of chloramphenicol by e coli
        return self.E.a * E / self.E.q

    def degradation_s(self, S):
        # Degradation of cefotaxime by salmonella
        return self.S.a * S / self.S.q

    def model(self, y, t):
        E, S, cefo, chloram = y
        if t >= self.cefo_start:
            # Considering that cefotaxime is added after two hours
            dE = E * (self.growth_e(E, S) - self.death_e(cefo))
            dCefo = -self.growth_s(E, S) * self.degradation_s(S) * cefo
        else:
            dE = self.growth_e(E, S) * E
            dCefo = 0
        dS = self.growth_s(E, S) * self.death_s(chloram) * S
        dChloram = -self.growth_e(E, S) * self.degradation_e(E) * chloram
        return dE, dS, dCefo, dChloram

    def transfer(self):
        for i in range(self.total_transfers):
            steps = 1000
            # For continuity of time we need to add previous time
            t0 = self.time[-1] if len(self.time) > 0 else 0
            t1 = np.linspace(t0, t0 + self.cefo_start, steps)
            y0 = [self.E.N0, self.S.N0, 0, self.M_chloram]
            # Integration over period before cefotaxime is added
            Y1 = odeint(self.model, y0, t1)
            self.time = np.concatenate((self.time, t1))
            self.E.N = np.concatenate((self.E.N, Y1[:, 0]))
            self.S.N = np.concatenate((self.S.N, Y1[:, 1]))
            self.cefo = np.concatenate((self.cefo, Y1[:, 2]))
            self.chloram = np.concatenate((self.chloram, Y1[:, 3]))
            # Integration after cefotaxime is added
            t2 = np.linspace(
                t1[-1], t1[-1] + (self.transfer_period - self.cefo_start), steps
            )
            y0 = [
                self.E.N[-1],
                self.S.N[-1],
                self.M_cefo,
                self.chloram[-1],
            ]  # Add antibiotic
            Y2 = odeint(self.model, y0, t2)

            self.time = np.concatenate((self.time, t2))
            self.E.N = np.concatenate((self.E.N, Y2[:, 0]))
            self.S.N = np.concatenate((self.S.N, Y2[:, 1]))
            self.cefo = np.concatenate((self.cefo, Y2[:, 2]))
            self.chloram = np.concatenate((self.chloram, Y2[:, 3]))

            # Dilution for next transfer cycle
            self.E.N0 = self.E.N[-1] / self.dilution_factor
            self.S.N0 = self.S.N[-1] / self.dilution_factor

    def plot_N(self, fname):
        # Plotting the abundance of e coli and salmonella
        fig = go.Figure(
            go.Scatter(
                x=self.time,
                y=self.E.N,
                mode="lines",
                name="<i>E. coli</i>",
                line=dict(color=colors["ecoli"]),
            )
        )
        fig.add_trace(
            go.Scatter(
                x=self.time,
                y=self.S.N,
                mode="lines",
                name="<i>ST</i>",
                line=dict(color=colors["st"]),
            )
        )
        fig.update_layout(
            xaxis=dict(
                # range=[0, max(self.time)],
                # dtick=12,
                title="Time [h]",
            ),
            yaxis=dict(
                # range=[0, 1.2],
                # dtick=0.2,
                title="Abundance [OD]",
            ),
            width=width,
            height=height,
        )
        fig = style_plot(
            fig,
            line_thickness=1.7,
            left_margin=30,
            top_margin=0,
            buttom_margin=30,
            right_margin=0,
        )
        fig.write_image("plots/transfers/" + fname + ".svg")

    def plot_cefo(self, fname):
        # Plotting the cefotaxime concentration
        fig = go.Figure(
            go.Scatter(
                x=self.time,
                y=self.cefo,
                mode="lines",
                name="Cefotaxime",
                line=dict(color=colors["cf"]),
            )
        )

        fig.update_layout(
            xaxis=dict(
                # range=[0, max(self.time)],
                # dtick=6,
                title="Time [h]",
            ),
            yaxis=dict(
                # range=[-10, 0.3],
                # dtick=0.05,
                title="Cf [µg/mL]",
            ),
            width=width,
            height=height,
        )
        fig = style_plot(
            fig,
            line_thickness=1.7,
            left_margin=30,
            top_margin=0,
            buttom_margin=25,
            right_margin=0,
        )
        fig.update_layout(width=width * 0.45)
        fig.write_image("plots/transfers/" + fname + ".svg")

    def plot_chloram(self, fname):
        # Plotting the chloramphenicol concentration
        fig = go.Figure(
            go.Scatter(
                x=self.time,
                y=self.chloram,
                mode="lines",
                name="Cefotaxime",
                line=dict(color=colors["ch"]),
            )
        )

        fig.update_layout(
            xaxis=dict(
                # range=[0, max(self.time)],
                # dtick=6,
                title="Time [h]",
            ),
            yaxis=dict(
                # range=[-10, 18],
                # dtick=4,
                title="Ch [µg/mL]",
            ),
            width=width,
            height=height,
        )
        fig = style_plot(
            fig,
            line_thickness=1.7,
            left_margin=25,
            top_margin=0,
            buttom_margin=25,
            right_margin=0,
        )
        fig.update_layout(width=width * 0.45)
        fig.write_image("plots/transfers/" + fname + ".svg")
