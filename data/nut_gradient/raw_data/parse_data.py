import pandas as pd
from os.path import join
from datetime import timedelta

raw_rep_1 = pd.read_excel("growthrate_data.xlsx", sheet_name="Repeat 1")
meta_rep_1_p1 = pd.read_excel("replicate_1_plate_1.xlsx", index_col=0)
meta_rep_1_p2 = pd.read_excel("replicate_1_plate_2.xlsx", index_col=0)
data_rep_1_p1 = "data_rep_1_plate_1.xlsx"
data_rep_1_p2 = "data_rep_1_plate_2.xlsx"


raw_rep_2 = pd.read_excel("growthrate_data.xlsx", sheet_name="Repeat 2")
meta_rep_2_p1 = pd.read_excel("replicate_2_plate_1.xlsx", index_col=0)
meta_rep_2_p2 = pd.read_excel("replicate_2_plate_2.xlsx", index_col=0)
data_rep_2_p1 = "data_rep_2_plate_1.xlsx"
data_rep_2_p2 = "data_rep_2_plate_2.xlsx"


raw_rep_3 = pd.read_excel("growthrate_data.xlsx", sheet_name="Repeat 3")
meta_rep_3_p1 = pd.read_excel("replicate_3_plate_1.xlsx", index_col=0)
meta_rep_3_p2 = pd.read_excel("replicate_3_plate_2.xlsx", index_col=0)
data_rep_3_p1 = "data_rep_3_plate_1.xlsx"
data_rep_3_p2 = "data_rep_3_plate_2.xlsx"

rep_1 = [raw_rep_1, [meta_rep_1_p1, meta_rep_1_p2], [data_rep_1_p1, data_rep_1_p2]]
rep_2 = [raw_rep_2, [meta_rep_2_p1, meta_rep_2_p2], [data_rep_2_p1, data_rep_2_p2]]
rep_3 = [raw_rep_3, [meta_rep_3_p1, meta_rep_3_p2], [data_rep_3_p1, data_rep_3_p2]]

reps = [rep_1, rep_2, rep_3]

for rep in reps:
    raw, metas, datas = rep
    for meta, data in zip(metas, datas):
        meta = meta.dropna(axis=0, how="all")
        meta = meta.dropna(axis=1, how="all")
        wells = []
        for c in meta.columns:
            for i in meta.index:
                wells.append(i + str(c))
        df = pd.DataFrame(columns=["Time"] + wells)
        df["Time"] = raw["Time"]
        for c in meta.columns:
            for i in meta.index:
                cname = meta.at[i, c]
                well = i + str(c)
                df[well] = raw[str(cname)]
        df = df.reindex(["Time"] + sorted(df.columns[1:]), axis=1)
        df.to_excel(data, index=False)
