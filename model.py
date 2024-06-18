import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt


class Salmonella:
    def __init__(self):
        # max. growth rate
        self.r = 0.5
        # monod constant
        self.Km = 10
        # yield
        self.Y = 0.2
        # carrying capacity
        self.K = 0.8
        # energy investment
        self.f = 0.1
        # toxin monod constance
        self.KT = 16e-3 / 8
        # toxin degradation rate
        self.a = 500
        # passive toxin uptake rate
        self.j = 0.1
        # chloram death rate, MIC=0.6
        self.u = 0.33
        # population density
        self.N = np.array([0.1])


class E_coli:
    def __init__(self):
        # max. growth rate
        self.r = 0.55
        # monod constant
        self.Km = 10
        # yield
        self.Y = 0.2
        # carrying capacity
        self.K = 1
        # energy investment
        self.f = 0.1
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


# eq.e_min_death_rate()
class Experiment:
    def __init__(self):
        # transfers
        self.total_transfers = 8
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

    def logistic(self, y, t):
        # simple logistic model
        e, s = y
        dE = self.e.r * (self.e.K - e) / self.e.K
        dS = self.s.r * (self.s.K - s) / self.s.K
        return dE, dS

    def e_toxin(self, y, t):
        e, R, cefo = y
        je = self.e.r * R / (R + self.e.Km)
        ue = self.e.u * self.M_cefo / (self.M_cefo + self.e.KT)
        dR = -je * e / self.e.Y
        if t < 2:
            dE = je * e
            dCefo = 0
        else:
            dE = je * e - ue * e
            dCefo = -cefo * self.e.j * e
        return dE, dR, dCefo

    def s_toxin(self, y, t):
        s, R, chloram = y
        js = self.s.r * R / (R + self.s.Km) * self.s.KT / (self.s.KT + chloram)
        dR = -js * s / self.s.Y
        dS = js * s
        dChloram = -chloram * self.s.j * s
        return dS, dR, dChloram

    def e_protects_s(self, y, t):
        e, s, R, chloram = y
        je = self.e.r * R / (R + self.e.Km)
        js = self.s.r * R / (R + self.s.Km) * self.s.KT / (self.s.KT + chloram)
        dE = ((1 - self.e.f) * je) * e
        dS = js * s
        dR = -je * e / self.e.Y - js * s / self.s.Y
        dChloram = -chloram * (self.e.f * self.e.a * je + self.e.j + self.s.j) * e
        return dE, dS, dR, dChloram

    def s_protects_e(self, y, t):
        e, s, R, cefo = y
        je = self.e.r * R / (R + self.e.Km)
        js = self.s.r * R / (R + self.s.Km)
        ue = self.e.u * cefo / (cefo + self.e.KT)
        if t < 2:
            dE = je * e
        else:
            dE = je * e - ue * e
        dS = ((1 - self.s.f) * js) * s
        dR = -je * e / self.e.Y - js * s / self.s.Y
        if t > 2:
            dCefo = -cefo * (self.s.f * self.s.a * js + self.s.j + self.e.j) * s
        else:
            dCefo = 0
        return dE, dS, dR, dCefo

    def two_sided_protection(self, y, t):
        e, s, R, cefo, chloram = y
        je = self.e.r * R / (R + self.e.Km)
        js = self.s.r * R / (R + self.s.Km) * self.s.KT / (self.s.KT + chloram)
        ue = self.e.u * cefo / (cefo + self.e.KT)
        if t < 2:
            dE = je * e
        else:
            dE = je * e - ue * e
        dS = js * s
        dR = -je * e / self.e.Y - js * s / self.s.Y
        if t > 2:
            dCefo = -cefo * (self.s.f * self.s.a * js + self.s.j + self.e.j) * s
        else:
            dCefo = 0
        dChloram = -chloram * (self.e.f * self.e.a * je + self.e.j + self.s.j) * e
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
                    [self.e.N[-1], self.M, self.M_cefo],
                    t,
                )
                self.e.N = np.concatenate((self.e.N, Y[:, 0]))
                self.R = np.concatenate((self.R, Y[:, 1]))
                self.cefo = np.concatenate((self.cefo, Y[:, 2]))

            if model == "s_toxin":
                Y = odeint(
                    self.s_toxin,
                    [self.s.N[-1], self.M, self.M_chloram],
                    t,
                )
                self.s.N = np.concatenate((self.s.N, Y[:, 0]))
                self.R = np.concatenate((self.R, Y[:, 1]))
                self.chloram = np.concatenate((self.chloram, Y[:, 2]))

            # dilution
            self.e.N[-1] = self.e.N[-1] / self.dilution_factor
            self.s.N[-1] = self.s.N[-1] / self.dilution_factor

            # add time
            self.time = np.concatenate((self.time, self.time[-1] + t), axis=0)

    def isoclines(self):
        Rs = np.linspace(0,30)
        J_e_coli = [self.e.r * R / (R + self.e.Km) for R in Rs]
        J_salmonella = [self.s.r * R / (R + self.s.Km) for R in Rs]
        plt.plot(Rs,J_e_coli,label='E_coli')
        plt.plot(Rs,J_salmonella,label='Salmonella')
        plt.xlabel('Resource [mM]')
        plt.ylabel('Per capita growth rate [1/h]')
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
        plt.legend()
        plt.show()

    def plot_chloram(self):
        plt.plot(self.time, self.chloram, label="Chloramphenicol [mM]")
        plt.legend()
        plt.show()

    def plot_cefo_chloram(self):
        plt.plot(self.time, self.cefo, label="Cefotaxime [mM]")
        plt.plot(self.time, self.chloram, label="Chloramphenicol [mM]")
        plt.legend()
        plt.show()

    def plot_R(self):
        plt.plot(self.time, self.R, label="Resource [mM]")
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


# e_coli_susceptible()
# salmonella_susceptible()
# consumer_resource()
# s_protects_e()
# e_protects_s()
# two_sided()
