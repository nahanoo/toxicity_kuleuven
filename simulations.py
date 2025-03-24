from model import Experiment
import numpy as np
from scipy.integrate import odeint


class Simulations:
    def __init__(self):
        self.e = Experiment()

    def transfer(self):
        for i in range(self.e.total_transfers):
            steps = 1000
            # For continuity of time we need to add previous time
            t0 = self.e.time[-1] if len(self.e.time) > 0 else 0
            t1 = np.linspace(t0, t0 + self.e.cefo_start, steps)
            y0 = [self.e.E.N0, self.e.S.N0, 0, self.e.M_chloram]
            # Integration over period before cefotaxime is added
            Y1 = odeint(self.e.model, y0, t1)
            self.e.time = np.concatenate((self.e.time, t1))
            self.e.E.N = np.concatenate((self.e.E.N, Y1[:, 0]))
            self.e.S.N = np.concatenate((self.e.S.N, Y1[:, 1]))
            self.e.cefo = np.concatenate((self.e.cefo, Y1[:, 2]))
            self.e.chloram = np.concatenate((self.e.chloram, Y1[:, 3]))
            # Integration after cefotaxime is added
            t2 = np.linspace(
                t1[-1], t1[-1] + (self.e.transfer_period - self.e.cefo_start), steps
            )
            y0 = [
                self.e.E.N[-1],
                self.e.S.N[-1],
                self.e.M_cefo,
                self.e.chloram[-1],
            ]  # Add antibiotic
            Y2 = odeint(self.e.model, y0, t2)

            self.e.time = np.concatenate((self.e.time, t2))
            self.e.E.N = np.concatenate((self.e.E.N, Y2[:, 0]))
            self.e.S.N = np.concatenate((self.e.S.N, Y2[:, 1]))
            self.e.cefo = np.concatenate((self.e.cefo, Y2[:, 2]))
            self.e.chloram = np.concatenate((self.e.chloram, Y2[:, 3]))

            # Dilution for next transfer cycle
            self.e.E.N0 = self.e.E.N[-1] / self.e.dilution_factor
            self.e.S.N0 = self.e.S.N[-1] / self.e.dilution_factor


def no_antibiotics():
    # Setting antibitoics inflow to 0
    s = Simulations()
    s.e.M_cefo = 0
    s.e.M_chloram = 0
    s.transfer()
    s.e.plot_N("no_antibiotics/abundance")


def s_protects_e():
    # Only adding cefotaxime thus setting chloramphenicol inflow to 0
    s = Simulations()
    s.e.M_chloram = 0
    # Case with no degradation, setting degradation rate to 0
    s.e.S.a = 0
    s.transfer()
    s.e.plot_N("s_protects_e/abundance_no_degradation")
    s.e.plot_cefo("s_protects_e/cefotaxime_no_degradation")

    # Case considering degradation
    s = Simulations()
    s.e.M_chloram = 0
    s.transfer()
    s.e.plot_N("s_protects_e/abundance")
    s.e.plot_cefo("s_protects_e/cefotaxime")


def e_protects_s():
    # Only adding chloramphenicol thus setting cefotaxime inflow to 0
    s = Simulations()
    s.e.M_cefo = 0
    # Case with no degradation, setting degradation rate to 0
    s.e.E.a = 0
    s.transfer()
    s.e.plot_N("e_protects_s/abundance_no_degradation")
    s.e.plot_chloram("e_protects_s/chloram_no_degradation")

    s = Simulations()
    s.e.M_cefo = 0
    # Case with degradation
    s.transfer()
    s.e.plot_N("e_protects_s/abundance")
    s.e.plot_chloram("e_protects_s/chloram")


def two_sided_protection():
    s = Simulations()
    # Case with no degradation, setting degradation rate to 0
    s.e.E.a = 0
    s.e.S.a = 0
    s.transfer()
    s.e.plot_N("two_sided_protection/abundance_no_degradation")

    s = Simulations()
    s.transfer()
    s.e.plot_N("two_sided_protection/abundance")
    s.e.plot_cefo("two_sided_protection/cefotaxime")
    s.e.plot_chloram("two_sided_protection/chloram")
