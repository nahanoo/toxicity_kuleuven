from model import Experiment
import numpy as np
from scipy.integrate import odeint


class Simulations:
    def __init__(self):
        self.e = Experiment()

    def transfer(self):
        for i in range(self.e.total_transfers):
            steps = 1000
            t = np.linspace(0, self.e.cefo_start, steps)
            Y = odeint(
                self.e.model,
                [self.e.E.N[-1], self.e.S.N[-1], 0, self.e.M_chloram],
                t,
            )
            self.e.time = np.concatenate((self.e.time, self.e.time[-1] + t), axis=0)
            self.e.E.N = np.concatenate((self.e.E.N, Y[:, 0]))
            self.e.S.N = np.concatenate((self.e.S.N, Y[:, 1]))
            self.e.cefo = np.concatenate((self.e.cefo, Y[:, 2]))
            self.e.chloram = np.concatenate((self.e.chloram, Y[:, 3]))
            t = np.linspace(self.e.cefo_start, self.e.transfer_period, steps)
            Y = odeint(
                self.e.model,
                [self.e.E.N[-1], self.e.S.N[-1], self.e.M_cefo, self.e.chloram[-1]],
                t,
            )
            self.e.time = np.concatenate(
                (self.e.time, self.e.time[-1] + t - self.e.cefo_start), axis=0
            )
            self.e.E.N = np.concatenate((self.e.E.N, Y[:, 0]))
            self.e.S.N = np.concatenate((self.e.S.N, Y[:, 1]))
            self.e.cefo = np.concatenate((self.e.cefo, Y[:, 2]))
            self.e.chloram = np.concatenate((self.e.chloram, Y[:, 3]))
            self.e.E.N[-1] = self.e.E.N[-1] / self.e.dilution_factor
            self.e.S.N[-1] = self.e.S.N[-1] / self.e.dilution_factor


def no_antibiotics():
    s = Simulations()
    s.e.M_cefo = 0
    s.e.M_chloram = 0
    s.transfer()
    s.e.plot_N("no_antibiotics/abundance")


def s_protects_e():
    s = Simulations()
    s.e.M_chloram = 0
    s.e.E.a = 0
    s.e.S.a = 0

    s.transfer()
    s.e.plot_N("s_protects_e/abundance_no_degradation")
    s.e.plot_cefo("s_protects_e/cefotaxime_no_degradation")

    s = Simulations()
    s.e.M_chloram = 0
    s.transfer()
    s.e.plot_N("s_protects_e/abundance")
    s.e.plot_cefo("s_protects_e/cefotaxime")


def e_protects_s():
    s = Simulations()
    s.e.M_cefo = 0
    s.transfer()
    s.e.plot_N("e_protects_s/abundance_no_degradation")
    s.e.plot_chloram("e_protects_s/chloram_no_degradation")

    s = Simulations()
    s.e.M_cefo = 0
    s.transfer()
    s.e.plot_N("e_protects_s/abundance")
    s.e.plot_chloram("e_protects_s/chloram")


def two_sided_protection():
    s = Simulations()
    s.e.E.a = 0
    s.e.S.a = 0
    s.transfer()
    s.e.plot_N("two_sided_protection/abundance_no_degradation")

    s = Simulations()
    s.transfer()
    s.e.plot_N("two_sided_protection/abundance")
    s.e.plot_cefo("two_sided_protection/cefotaxime")
    s.e.plot_chloram("two_sided_protection/chloram")


s_protects_e()
