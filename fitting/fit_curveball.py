import curveball
import numpy as np
import pandas as pd
from os.path import join
import plotly.graph_objects as go

model_names = [
    "Model(Richards)",
    "Model(RichardsLag1)",
    "Model(BaranyiRoberts)",
    "Model(Logistic)",
    "Model(LogisticLag1)",
    "Model(LogisticLag2)",
]


def filter_time(xs, ys, x0=0, x1=5):
    xs_f = []
    ys_f = []
    for xi, yi in zip(xs, ys):
        mask = (xi >= x0) & (xi <= x1)
        xs = xi[mask]
        ys = yi[mask]
        ys = ys - np.min(ys) + 0.0025
        xs_f.append(xs)
        ys_f.append(ys)
    return xs_f, ys_f


def get_lgs(strain, conc):
    f = join(
        "~", "toxicity_kuleuven", "data", "nut_gradient", "pooled_df_joint_metadata.csv"
    )
    pooled_meta = pd.read_csv(f)
    filter = (pooled_meta["project"] == "240903_fran") & (
        pooled_meta["Experiment description"].isin(
            [
                "Repeat 1 of TSB gradient experiment with all strains.",
                "Repeat 2 of TSB gradient experiment with all strains.",
                "Repeat 3 of TSB gradient experiment with all strains.",
            ]
        )
        & (pooled_meta["cs_conc"] == conc)
        & (pooled_meta["species"] == strain)
    )
    meta = pooled_meta[filter]
    return meta["linegroup"].to_list()


def mask_curve(xs, ys, x0=0, x1=5):
    xs_f = []
    ys_f = []
    for xi, yi in zip(xs, ys):
        mask = (xi >= x0) & (xi <= x1)
        xs = xi[mask]
        ys = yi[mask]
        ys = ys - np.min(ys) + 0.0025
        xs_f.append(xs)
        ys_f.append(ys)
    return xs_f, ys_f


def fit_curveball(strain, conc, fix_K=True, K=1):
    f = join("~", "toxicity_kuleuven", "data", "nut_gradient", "measurement_data.csv")
    raw = pd.read_csv(f)
    lgs = get_lgs(strain, conc)
    xs = [raw[raw["linegroup"] == lg]["time"].to_numpy() for lg in lgs]
    ys = [raw[raw["linegroup"] == lg]["measurement"].to_numpy() for lg in lgs]
    xs, ys = mask_curve(xs, ys)
    dfs = []
    for i in range(3):
        df = pd.DataFrame(
            {"Time": xs[i], "OD": ys[i], "Well": str(i), "Strain": strain}
        )
        dfs.append(df)
    if fix_K:
        models = curveball.models.fit_model(
            pd.concat(dfs),
            PLOT=False,
            PRINT=False,
            param_guess={"K": K, "y0": 0.0025},
            param_fix=["K", "y0"],
        )
    else:
        models = curveball.models.fit_model(
            pd.concat(dfs),
            PLOT=True,
            PRINT=False,
            param_guess={"y0": 0.0025},
            param_fix=["y0"],
        )
    models = sorted(models, key=lambda r: r.model.name)
    return models


m = fit_curveball(
    "Salmonella Typhimurium mcherry:pGDPI CTX-M-15", conc=1, fix_K=True, K=1
)
