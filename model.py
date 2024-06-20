import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from sympy import symbols, solve, Eq, init_printing

R, chloram, cefo = symbols("R,chloram,cefo")


class Salmonella:
    def __init__(self):
        # max. growth rate
        self.r = 0.5
        # monod constant
        self.Km = 10
        # yield
        self.Y = 0.2
        # toxin monod constance
        self.KT = 16e-3 / 8
        # toxin degradation rate
        self.a = 500
        # chloram death rate, MIC=0.6
        self.u = 0.33
        # population density
        self.N = np.array([0.1])

        # functions
        # resource growth
        self.J = self.r * R / (R + self.Km)
        # death
        self.D = self.KT / (self.KT + chloram)


class E_coli:
    def __init__(self):
        # max. growth rate
        self.r = 0.55
        # monod constant
        self.Km = 10
        # yield
        self.Y = 0.2
        # toxin monod constance
        self.KT = 0.25e-3 / 8
        # toxin degradation rate
        self.a = 100
        # passive toxin uptake rate
        self.j = 0.1
        # cefo death rate, mic = 0.7
        self.u = 0.48

        # population density
        self.N = np.array([0.1])

        # functions
        # resource growth
        self.J = self.r * R / (R + self.Km)
        # death
        self.D = self.u * cefo / (cefo + self.KT)


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
        self.M_cefo = 0.25e-3
        # chloramphenicol concentration
        self.M_chloram = 16e-3
        # array for resource concentration
        self.R = np.array([self.M])
        # array for cefotaxime concentration
        self.cefo = np.array([self.M_cefo])
        # array for chlorampheicol concentration
        self.chloram = np.array([self.M_chloram])
        # array for time
        self.time = np.array([0])
        # create bugs
        self.s = Salmonella()
        self.e = E_coli()

    def cr(self, y, t):
        # simple consumer resource model with no toxins yet
        e, s, R = y
        dE = self.e.r * R / (R + self.e.Km) * e
        dS = self.s.r * R / (R + self.s.Km) * s
        dR = (
            -self.e.r * R / (R + self.e.Km) * e / self.e.Y
            - self.s.r * R / (R + self.s.Km) * s / self.s.Y
        )
        return dE, dS, dR

    def e_toxin(self, y, t):
        e, R = y
        J = self.e.J.subs({"R": R})
        D = self.e.D.subs({"cefo": self.M_cefo})
        dR = -J * e / self.e.Y
        if t > 2:
            dE = (J - D) * e
        else:
            dE = J * e
        return dE, dR

    def s_toxin(self, y, t):
        s, R = y
        J = self.s.J.subs({"R": R})
        D = self.s.D.subs({"chloram": self.M_chloram})
        dR = -J * s / self.s.Y
        dS = (J * D) * s
        return dS, dR

    def e_protects_s(self, y, t):
        e, s, R, chloram = y
        JE = self.e.J.subs({"R": R})
        JS = self.s.J.subs({"R": R})
        DS = self.s.D.subs({"chloram": chloram})
        dE = JE * e
        dS = (JS * DS) * s
        dR = -JE * e / self.e.Y - JS * s / self.s.Y
        dChloram = -chloram * self.e.a * e
        return dE, dS, dR, dChloram

    def s_protects_e(self, y, t):
        e, s, R, cefo = y
        JE = self.e.J.subs({"R": R})
        JS = self.s.J.subs({"R": R})
        DE = self.e.D.subs({"cefo": cefo})
        if t < 2:
            dE = JE * e
        else:
            dE = (JE - DE) * e
        dS = JS * s
        dR = -JE * e / self.e.Y - JS * s / self.s.Y
        if t > 2:
            dCefo = -cefo * self.s.a * s
        else:
            dCefo = 0
        return dE, dS, dR, dCefo

    def two_sided_protection(self, y, t):
        e, s, R, cefo, chloram = y
        JE = self.e.J.subs({"R": R})
        JS = self.s.J.subs({"R": R})
        DE = self.e.D.subs({"cefo": cefo})
        DS = self.s.D.subs({"chloram": chloram})
        if t < 2:
            dE = JE * e
        else:
            dE = (JE - DE) * e
        dS = (JS * DS) * s
        dR = -JE * e / self.e.Y - JS * s / self.s.Y
        if t > 2:
            dCefo = -cefo * self.s.a * s
        else:
            dCefo = 0
        dChloram = -chloram * self.e.a * e
        return dE, dS, dR, dCefo, dChloram

    def simulate_experiment(self, model):
        # function to simulate transfers
        for i in range(self.total_transfers):
            steps = 1000
            t = np.linspace(0, self.transfer_period, steps)
            if model == "consumer_resource":
                Y = odeint(self.cr, [self.e.N[-1], self.s.N[-1], self.M], t)
                self.R = np.concatenate((self.R, Y[:, 2]))
                self.e.N = np.concatenate((self.e.N, Y[:, 0]))
                self.s.N = np.concatenate((self.s.N, Y[:, 1]))
            # growth
            if model == "logistic":
                Y = odeint(self.logistic, [self.e.N[-1], self.s.N[-1]], t)
            # add population densities
            if model == "e_protects_s":
                Y = odeint(
                    self.e_protects_s,
                    [self.e.N[-1], self.s.N[-1], self.M, self.M_chloram],
                    t,
                )
                self.R = np.concatenate((self.R, Y[:, 2]))
                self.chloram = np.concatenate((self.chloram, Y[:, 3]))
                self.e.N = np.concatenate((self.e.N, Y[:, 0]))
                self.s.N = np.concatenate((self.s.N, Y[:, 1]))
            if model == "s_protects_e":
                Y = odeint(
                    self.s_protects_e,
                    [self.e.N[-1], self.s.N[-1], self.M, self.M_cefo],
                    t,
                )
                self.R = np.concatenate((self.R, Y[:, 2]))
                self.cefo = np.concatenate((self.cefo, Y[:, 3]))
                self.e.N = np.concatenate((self.e.N, Y[:, 0]))
                self.s.N = np.concatenate((self.s.N, Y[:, 1]))
            if model == "two_sided":
                Y = odeint(
                    self.two_sided_protection,
                    [self.e.N[-1], self.s.N[-1], self.M, self.M_cefo, self.M_chloram],
                    t,
                )
                self.R = np.concatenate((self.R, Y[:, 2]))
                self.cefo = np.concatenate((self.cefo, Y[:, 3]))
                self.chloram = np.concatenate((self.chloram, Y[:, 4]))
                self.e.N = np.concatenate((self.e.N, Y[:, 0]))
                self.s.N = np.concatenate((self.s.N, Y[:, 1]))
            if model == "e_toxin":
                Y = odeint(
                    self.e_toxin,
                    [self.e.N[-1], self.M],
                    t,
                )
                self.e.N = np.concatenate((self.e.N, Y[:, 0]))
                self.R = np.concatenate((self.R, Y[:, 1]))
            if model == "s_toxin":
                Y = odeint(
                    self.s_toxin,
                    [self.s.N[-1], self.M],
                    t,
                )
                self.s.N = np.concatenate((self.s.N, Y[:, 0]))
                self.R = np.concatenate((self.R, Y[:, 1]))
            # dilution
            self.e.N[-1] = self.e.N[-1] / self.dilution_factor
            self.s.N[-1] = self.s.N[-1] / self.dilution_factor

            # add time
            self.time = np.concatenate((self.time, self.time[-1] + t), axis=0)

    def isoclines(self):
        Rs = np.linspace(0, 30)
        J_e_coli = [self.e.r * R / (R + self.e.Km) for R in Rs]
        J_salmonella = [self.s.r * R / (R + self.s.Km) for R in Rs]
        plt.plot(Rs, J_e_coli, label="E_coli")
        plt.plot(Rs, J_salmonella, label="Salmonella")
        plt.xlabel("Resource [mM]")
        plt.ylabel("Per capita growth rate [1/h]")
        plt.legend()
        plt.show()

    def plot_N(self):
        plt.plot(self.time, self.e.N, label="E_coli")
        plt.plot(self.time, self.s.N, label="Salmonella")
        plt.xlabel("Time [h]"), plt.ylabel("OD")
        plt.legend()
        plt.show()

    def plot_e_coli(self):
        plt.plot(self.time, self.e.N, label="E_coli")
        plt.xlabel("Time [h]"), plt.ylabel("OD")
        plt.legend()
        plt.show()

    def plot_salmonella(self):
        plt.plot(self.time, self.s.N, label="Salmonella [mM]")
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


from matplotlib import pyplot as plt


def s_protects_e():
    e = Experiment()
    C, T = np.meshgrid(np.linspace(0, 10, 50), np.linspace(0, 0.25e-3, 50))
    J_E = e.e.r * C / (C + e.e.Km) - e.e.u * T / (T + e.e.KT)
    J_S = e.s.r * C / (C + e.s.Km)
    plt.contourf(C, T, J_E, levels=50)
    plt.show()
    plt.contour(C, T, J_E, levels=np.linspace(0, 0.5, 10))
    plt.contour(C, T, J_S, levels=np.linspace(0, 0.5, 10))
    plt.show()


def e_protects_s():
    e = Experiment()
    C, T = np.meshgrid(np.linspace(0, 30, 50), np.linspace(0, 16e-3, 50))
    J_E = e.e.r * C / (C + e.e.Km)
    J_S = e.s.r * C / (C + e.s.Km) * e.s.KT / (e.s.KT + T)
    plt.contourf(C, T, J_S, levels=50)
    plt.show()
    plt.contour(C, T, J_E, levels=np.linspace(0, 0.5, 10))
    plt.contour(C, T, J_S, levels=np.linspace(0, 0.5, 10))

    plt.show()


e = Experiment()
e.total_transfers = 1
e.transfer_period = 24
e.dilution_factor = 1
e.simulate_experiment("two_sided")
e.plot_N()
e.plot_cefo()
