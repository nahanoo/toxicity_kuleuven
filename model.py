import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt


class Salmonella:
    def __init__(self):
        # max. growth rate
        self.r = 0.4
        # monod constant
        self.Km = 10
        # yield
        self.Y = 0.025
        # carrying capacity
        self.K = 0.8
        # population density
        self.N = np.array([0.05])


class E_coli:
    def __init__(self):
        # max. growth rate
        self.r = 0.45
        # monod constant
        self.Km = 10
        # yield
        self.Y = 0.025
        # carrying capacity
        self.K = 1
        # population density
        self.N = np.array([0.05])


class Experiment:
    def __init__(self):
        # transfers
        self.total_transfers = 8
        # dilution factor
        self.dilution_factor = 100
        # transfer period
        self.transfer_period = 24
        # media concentration
        self.M = 40
        # array for resource concentration
        self.R = np.array([self.M])
        # array for time
        self.time = np.array([0])
        # create bugs
        self.s = Salmonella()
        self.e = E_coli()

    def cr(self, y, t):
        # simple consumer resource model with no toxins yet
        e, s, R = y[0], y[1], y[2]
        dE = self.e.r * R / (R + self.e.Km) * e
        dS = self.s.r * R / (R + self.s.Km) * s
        dR = (
            -self.e.r * R / (R + self.e.Km) * e / self.e.Y
            - self.s.r * R / (R + self.s.Km) * s / self.s.Y
        )
        return dE, dS, dR

    def logistic(self, y, t):
        # simple logistic model
        e, s = y[0], y[1]
        dE = self.e.r * (self.e.K - e) / self.e.K
        dS = self.s.r * (self.s.K - s) / self.s.K
        return dE, dS

    def simulate_experiment(self, model):
        # function to simulate transfers
        for i in range(self.total_transfers):
            t = np.linspace(0, self.transfer_period, 100)
            if model == "consumer_resource":
                Y = odeint(self.cr, [self.e.N[-1], self.s.N[-1], self.M], t)
                self.R = np.concatenate((self.R, Y[:, 2]))
            # growth
            if model == "logistic":
                Y = odeint(self.logistic, [self.e.N[-1], self.s.N[-1]], t)
            # add population densities
            self.e.N = np.concatenate((self.e.N, Y[:, 0]))
            self.s.N = np.concatenate((self.s.N, Y[:, 1]))

            # dilution
            self.e.N[-1] = self.e.N[-1] / self.dilution_factor
            self.s.N[-1] = self.s.N[-1] / self.dilution_factor

            # add time
            self.time = np.concatenate((self.time, self.time[-1] + t), axis=0)

    def plotting(self):
        plt.plot(self.time, self.e.N, label="E_coli")
        plt.plot(self.time, self.s.N, label="Salmonella")
        plt.xlabel("Time [h]"), plt.ylabel("OD")
        plt.legend()
        plt.show()


def main(model):
    exp = Experiment()
    exp.simulate_experiment(model)
    exp.plotting()


main("logistic")
main("consumer_resource")
