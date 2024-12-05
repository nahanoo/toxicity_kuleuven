import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from sympy import symbols, solve, Eq, init_printing
 

class Salmonella:

    def __init__(self):
        # max. growth rate
        self.r = 0.8
        # monod constant
        self.Km = 10
        # yield
        self.Y = 0.1
        # toxin monod constance
        self.KT = 0.07256629755179302
        # toxin degradation rate
        self.a = 50
        # population density
        self.N = np.array([0.1])
        # passive uptake rate

    # growth
    def J(self, R):
        return self.r * R / (R + self.Km)

    # death
    def D(self, chloram):
        # return self.KT / (self.KT + chloram)
        return 1 / (1 + chloram / self.KT)

    # detoxification
    def detox(self, cefo, S, R):
        return -cefo * self.a * S


class E_coli:

    def __init__(self):
        # max. growth rate
        self.r = 0.9
        # monod constant
        self.Km = 10
        # yield
        self.Y = 0.1
        # toxin monod constance
        self.KT = 0.00030537131107306284
        # toxin degradation rate
        self.a = 1
        # cefo death rate, mic = 0.7
        self.u = 0.5

        # population density
        self.N = np.array([0.1])

    # growth
    def J(self, R):
        return self.r * R / (R + self.Km)

    # death
    def D(self, cefo):
        return self.u * cefo / (cefo + self.KT)

    def detox(self, chloram, E, R):
        return -chloram * self.a * E


# eq.e_min_death_rate()
class Experiment:

    def __init__(self):
        # transfers
        self.total_transfers = 4
        # dilution factor
        self.dilution_factor = 100
        # transfer period
        self.transfer_period = 24
        # media concentration
        self.M = 10
        # cefotaxime concentration
        self.M_cefo = 0.25
        # chloramphenicol concentration
        self.M_chloram = 16
        # array for resource concentration
        self.R = np.array([self.M])
        # array for cefotaxime concentration
        self.cefo = np.array([self.M_cefo])
        # array for chlorampheicol concentration
        self.chloram = np.array([self.M_chloram])
        # array for time
        self.time = np.array([0])
        # create bugs
        self.S = Salmonella()
        self.E = E_coli()

    def cr(self, y, t):
        # simple consumer resource model with no toxins yet
        E, S, R = y
        dE = self.E.J(R) * E
        dS = self.S.J(R) * S
        dR = -self.E.J(R) / self.E.Y * E - self.S.J(R) / self.S.Y * S
        return dE, dS, dR

    def e_toxin(self, y, t):
        E, R = y
        dR = -self.E.J(R) / self.E.Y * E
        if t > 2:
            dE = (self.E.J(R) - self.E.D(self.M_cefo)) * E
        else:
            dE = self.E.J(R) * E
        return dE, dR

    def s_toxin(self, y, t):
        S, R = y
        dR = self.S.J(R) * S / self.S.Y
        dS = self.S.J(R) * self.S.D(self.M_chloram) * S
        return dS, dR

    def e_protects_s(self, y, t):
        E, S, R, chloram = y
        dE = self.E.J(R) * E
        dS = self.S.J(R) * self.S.D(chloram) * S
        dR = -self.E.J(R) * E / self.E.Y - self.S.J(R) * S / self.S.Y
        dChloram = self.E.detox(chloram, E, R)
        return dE, dS, dR, dChloram

    def s_protects_e(self, y, t):
        E, S, R, cefo = y
        if t < 2:
            dE = self.E.J(R) * E
        else:
            dE = (self.E.J(R) - self.E.D(cefo)) * E
        dS = self.S.J(R) * S
        dR = -self.E.J(R) * E / self.E.Y - self.S.J(R) * S / self.S.Y
        if t > 2:
            dCefo = self.S.detox(cefo, S, R)
        else:
            dCefo = 0
        return dE, dS, dR, dCefo

    def two_sided_protection(self, y, t):
        E, S, R, cefo, chloram = y
        if t < 2:
            dE = self.E.J(R) * E
        else:
            dE = (self.E.J(R) - self.E.D(cefo)) * E
        dS = self.S.J(R) * self.S.D(chloram) * S
        dR = -self.E.J(R) * E / self.E.Y - self.S.J(R) * S / self.S.Y
        if t > 2:
            dCefo = self.S.detox(cefo, S, R)
        else:
            dCefo = 0
        dChloram = self.E.detox(chloram, E, R)
        return dE, dS, dR, dCefo, dChloram

    def simulate_experiment(self, model):
        # function to simulate transfers
        for i in range(self.total_transfers):
            steps = 1000
            t = np.linspace(0, self.transfer_period, steps)
            if model == "consumer_resource":
                Y = odeint(self.cr, [self.E.N[-1], self.S.N[-1], self.M], t)
                self.R = np.concatenate((self.R, Y[:, 2]))
                self.E.N = np.concatenate((self.E.N, Y[:, 0]))
                self.S.N = np.concatenate((self.S.N, Y[:, 1]))
            # growth
            if model == "logistic":
                Y = odeint(self.logistic, [self.E.N[-1], self.S.N[-1]], t)
            # add population densities
            if model == "e_protects_s":
                Y = odeint(
                    self.e_protects_s,
                    [self.E.N[-1], self.S.N[-1], self.M, self.M_chloram],
                    t,
                )
                self.R = np.concatenate((self.R, Y[:, 2]))
                self.chloram = np.concatenate((self.chloram, Y[:, 3]))
                self.E.N = np.concatenate((self.E.N, Y[:, 0]))
                self.S.N = np.concatenate((self.S.N, Y[:, 1]))
            if model == "s_protects_e":
                Y = odeint(
                    self.s_protects_e,
                    [self.E.N[-1], self.S.N[-1], self.M, self.M_cefo],
                    t,
                )
                self.R = np.concatenate((self.R, Y[:, 2]))
                self.cefo = np.concatenate((self.cefo, Y[:, 3]))
                self.E.N = np.concatenate((self.E.N, Y[:, 0]))
                self.S.N = np.concatenate((self.S.N, Y[:, 1]))
            if model == "two_sided":
                Y = odeint(
                    self.two_sided_protection,
                    [self.E.N[-1], self.S.N[-1], self.M, self.M_cefo, self.M_chloram],
                    t,
                )
                self.R = np.concatenate((self.R, Y[:, 2]))
                self.cefo = np.concatenate((self.cefo, Y[:, 3]))
                self.chloram = np.concatenate((self.chloram, Y[:, 4]))
                self.E.N = np.concatenate((self.E.N, Y[:, 0]))
                self.S.N = np.concatenate((self.S.N, Y[:, 1]))
            if model == "e_toxin":
                Y = odeint(
                    self.e_toxin,
                    [self.E.N[-1], self.M],
                    t,
                )
                self.E.N = np.concatenate((self.E.N, Y[:, 0]))
                self.R = np.concatenate((self.R, Y[:, 1]))
            if model == "s_toxin":
                Y = odeint(
                    self.s_toxin,
                    [self.S.N[-1], self.M],
                    t,
                )
                self.S.N = np.concatenate((self.S.N, Y[:, 0]))
                self.R = np.concatenate((self.R, Y[:, 1]))
            # dilution
            self.E.N[-1] = self.E.N[-1] / self.dilution_factor
            self.S.N[-1] = self.S.N[-1] / self.dilution_factor

            # add time
            self.time = np.concatenate((self.time, self.time[-1] + t), axis=0)

    def plot_N(self):
        plt.plot(self.time, self.E.N, label="E_coli")
        plt.plot(self.time, self.S.N, label="Salmonella")
        plt.xlabel("Time [h]"), plt.ylabel("OD")
        plt.legend()
        plt.show()

    def plot_e_coli(self):
        plt.plot(self.time, self.E.N, label="E_coli")
        plt.xlabel("Time [h]"), plt.ylabel("OD")
        plt.legend()
        plt.show()

    def plot_salmonella(self):
        plt.plot(self.time, self.S.N, label="Salmonella [mM]")
        plt.xlabel("Time [h]"), plt.ylabel("OD")
        plt.legend()
        plt.show()

    def plot_cefo(self):
        plt.plot(self.time, self.cefo, label="Cefotaxime [mM]")
        plt.xlabel("Time [h]"), plt.ylabel("Cefotaxime [mM]")
        plt.legend()
        plt.show()

    def plot_chloram(self):
        plt.plot(self.time, self.chloram, label="Chloramphenicol [mM]")
        plt.xlabel("Time [h]"), plt.ylabel("Chloramphenicol [mM]")
        plt.legend()
        plt.show()

    def plot_cefo_chloram(self):
        plt.plot(self.time, self.cefo, label="Cefotaxime [mM]")
        plt.plot(self.time, self.chloram, label="Chloramphenicol [mM]")
        plt.xlabel("Time [h]"), plt.ylabel("Antibiotic concentration [mM]")
        plt.legend()
        plt.show()

    def plot_R(self):
        plt.plot(self.time, self.R, label="Resource [mM]")
        plt.xlabel("Time [h]"), plt.ylabel("Limitting resource concentration [mM]")
        plt.legend()
        plt.show()

    def isoclines(self):
        Rs = np.linspace(0, 30)
        E_isoclines = [self.E.J(R) for R in Rs]
        S_isoclines = [self.S.J(R) for R in Rs]
        plt.plot(Rs, E_isoclines, label="E_coli")
        plt.plot(Rs, S_isoclines, label="Salmonella")
        plt.xlabel("Resource [mM]")
        plt.ylabel("Per capita growth rate [1/h]")
        plt.legend()
        plt.show()

    def isoclines_inhibition_e_coli(self):
        R_grid, cefo_grid = np.meshgrid(
            np.linspace(0, self.M, 100), np.linspace(0, self.M_cefo, 100)
        )
        u = self.E.J(R_grid) - self.E.D(cefo_grid)
        plt.contourf(R_grid, cefo_grid, u, cmap="plasma", levels=100)
        plt.colorbar(label="μ [1/h]")
        plt.xlabel("Resource concentration [mM]")
        plt.ylabel("Cefotaxime concentration [mM]")
        plt.show()
        plt.close()

    def isoclines_inhibition_e_coli_salmonella(self):
        R_grid, cefo_grid = np.meshgrid(
            np.linspace(0, self.M, 100), np.linspace(0, self.M_cefo, 100)
        )
        uE = self.E.J(R_grid) - self.E.D(cefo_grid)
        uS = self.S.J(R_grid)
        plt.contourf(R_grid, cefo_grid, uE, cmap="plasma", levels=100)
        plt.colorbar(label="μ [1/h]")
        plt.contour(R_grid, cefo_grid, uS, levels=[0.2], colors="black")
        plt.contour(R_grid, cefo_grid, uE, levels=[0.2], colors="black")
        plt.xlabel("Resource concentration [mM]")
        plt.ylabel("Cefotaxime concentration [mM]")
        plt.show()
        plt.close()

    def isoclines_inhibition_salmonella(self):
        R_grid, chloram_grid = np.meshgrid(
            np.linspace(0, self.M, 100), np.linspace(0, self.M_chloram, 100)
        )
        u = self.S.J(R_grid) * self.S.D(chloram_grid)
        plt.contourf(R_grid, chloram_grid, u, cmap="plasma", levels=100)
        plt.colorbar(label="μ [1/h]")
        plt.xlabel("Resource concentration [mM]")
        plt.ylabel("Chloramphenicol concentration [mM]")
        plt.show()
        plt.close()

    def isoclines_inhibition_salmonella_e_coli(self):
        R_grid, chloram_grid = np.meshgrid(
            np.linspace(0, self.M, 100), np.linspace(0, self.M_chloram, 100)
        )
        uS = self.S.J(R_grid) * self.S.D(chloram_grid)
        uE = self.E.J(R_grid)
        plt.contourf(R_grid, chloram_grid, uS, cmap="plasma", levels=100)
        plt.colorbar(label="μ [1/h]")
        plt.contour(R_grid, chloram_grid, uS, levels=[0.2], colors="black")
        plt.contour(R_grid, chloram_grid, uE, levels=[0.2], colors="black")
        plt.xlabel("Resource concentration [mM]")
        plt.ylabel("Chloramphenicol concentration [mM]")
        plt.show()
        plt.close()

    def isoclines_inhibition_salmonella_inhibition_e_coli(self):
        R_grid, cefo_grid = np.meshgrid(
            np.linspace(0, self.M, 100), np.linspace(0, self.M_cefo, 100)
        )
        R_grid, chloram_grid = np.meshgrid(
            np.linspace(0, self.M, 100), np.linspace(0, self.M_chloram, 100)
        )
        uE = self.E.J(R_grid) - self.E.D(cefo_grid)
        plt.contour(R_grid, cefo_grid, uE, cmap="plasma", levels=[0.2])
        uS = self.S.J(R_grid) * self.S.D(chloram_grid)
        plt.contour(R_grid, chloram_grid, uS, cmap="plasma", levels=[0.2])
        plt.show()


def e_protects_s():
    exp = Experiment()
    exp.simulate_experiment("e_protects_s")
    exp.plot_N()
    # exp.plot_chloram()


def s_protects_e():
    exp = Experiment()
    exp.simulate_experiment("s_protects_e")
    exp.plot_N()
    # exp.plot_cefo()


def two_sided():
    exp = Experiment()
    exp.simulate_experiment("two_sided")
    exp.plot_N()
    # exp.plot_cefo_chloram()


def consumer_resource():
    exp = Experiment()
    exp.simulate_experiment("consumer_resource")
    exp.plot_N()


def e_coli_susceptible():
    exp = Experiment()
    exp.simulate_experiment("e_toxin")
    exp.plot_e_coli()


def salmonella_susceptible():
    exp = Experiment()
    exp.simulate_experiment("s_toxin")
    exp.plot_salmonella()


"""

import plotly.graph_objects as go

e = Experiment()
Rs = np.linspace(0, e.M, 100)
chlorams = np.linspace(0, e.M_chloram, 100)
cefos = np.linspace(0, e.M_cefo, 100)

R_grid, chloram_grid, cefo_grid = np.meshgrid(Rs, chlorams, cefos)
uS = e.S.J(R_grid) * e.S.D(chloram_grid)
uE = e.E.J(R_grid) - e.E.D(cefo_grid)
fig = go.Figure()
fig.add_trace(
    go.Isosurface(
        x=R_grid.flatten(),
        y=chloram_grid.flatten(),
        z=cefo_grid.flatten(),
        value=uS.flatten(),
        isomin=0.2,  # Set the value of interest (C)
        isomax=0.2,  # Keep isomin and isomax the same for one surface
        surface_count=1,  # Number of isosurfaces
        colorscale="Viridis",
        #caps=dict(x_show=False, y_show=False, z_show=False)
    ))
fig.add_trace(
    go.Isosurface(
        x=R_grid.flatten(),
        y=chloram_grid.flatten(),
        z=cefo_grid.flatten(),
        value=uE.flatten(),
        isomin=0.2,  # Set the value of interest (C)
        isomax=0.2,  # Keep isomin and isomax the same for one surface
        surface_count=1,  # Number of isosurfaces
        colorscale="Viridis",
        #caps=dict(x_show=False, y_show=False, z_show=False)
    ))
fig.update_layout(scene=dict(xaxis_title='Resource concentration',
                             yaxis_title='Chloramphenicol concentration [mM]',
                             zaxis_title='Cefotaxime concentration [mM]'), )
fig.show()
"""
