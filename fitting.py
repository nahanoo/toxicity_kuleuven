import pandas as pd
from os.path import join
import plotly.express as px
import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt
from scipy.integrate import odeint
import scipy.optimize as optimize


f = join("~", "curves", "metadata", "pooled_df_joint_metadata.csv")
pooled_meta = pd.read_csv(f)
f = join("~", "curves", "export", "240903_fran", "measurement_data.csv")
raw = pd.read_csv(f)

filter = (pooled_meta["project"] == "240903_fran") & (
    pooled_meta["Experiment description"].isin(
        [
            "Repeat 1 of TSB gradient experiment with all strains.",
            "Repeat 2 of TSB gradient experiment with all strains.",
            "Repeat 3 of TSB gradient experiment with all strains.",
        ]
    )
)

meta = pooled_meta[filter]

params = pd.DataFrame(columns=["growth_rate", "yield"])


def max_growth_rate():
    for species in set(meta["species"]):
        tmp = meta[meta["species"] == species]
        tsb_1 = tmp[tmp["cs_conc"] == 1]
        rs = []
        for lg in set(tsb_1["linegroup"]):
            xs = raw[raw["linegroup"] == lg]["time"].to_numpy()
            ys = raw[raw["linegroup"] == lg]["measurement"].to_numpy()
            for j, (x, y) in enumerate(zip(xs, ys)):
                if y < 0:
                    i = j
                if x < 5:
                    k = j
            x_window, y_window = xs[i + 1 : k], ys[i + 1 : k]
            rs.append(linregress(x_window, np.log(y_window))[0])
        params.at[species, "growth_rate"] = np.average(rs)


def get_yield():
    for species in set(meta["species"]):
        tmp = meta[meta["species"] == species]
        tsb_025 = tmp[tmp["cs_conc"] == 0.25]
        qs = []
        for lg in set(tsb_025["linegroup"]):
            max_OD = max(raw[raw["linegroup"] == lg]["measurement"].to_numpy())
            qs.append(max_OD / 0.25)
        params.at[species, "yield"] = np.average(qs)


def monod(y, t, v, Km, q):
    n, c = y
    dndt = v * c / (Km + c) * n
    dcdt = -v * c / (Km + c) * n / q
    return np.array([dndt, dcdt])


def simulate_monod(Km, v, t, q, n, c0, n0):
    y = odeint(monod, [n0, c0], t, args=(v, Km[0], q))
    return np.sum((n - y[:, 0]) ** 2)


for species in set(meta["species"]):
    tmp = meta[meta["species"] == species]
    tsb_1 = tmp[tmp["cs_conc"] == 1]
    rs = []
    for lg in set(tsb_1["linegroup"]):
        xs = raw[raw["linegroup"] == lg]["time"].to_numpy()
        ys = raw[raw["linegroup"] == lg]["measurement"].to_numpy()
        for j, (x, y) in enumerate(zip(xs, ys)):
            if y < 0:
                i = j
        xs, ys = xs[j:], ys[j:]
