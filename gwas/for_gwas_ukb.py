import numpy as np
import pandas as pd


# GP data path
path = "/data/GP/RecordLevel/"
# Prescription data path
path2 = "/data/prescription/"


def extract_GP_Biomarker(
    read2s=[],
    read3s=[],
    read2_path=path + "read2.csv",
    read3_path=path + "read3.csv",
    gp=path + "clinical_processed_recordlevel2.csv",
):
    if len(read2s) == 0 and len(read3s) == 0:
        return
    else:
        read2 = pd.read_csv(read2_path)
        read3 = pd.read_csv(read3_path)
        read2 = read2.dropna().reset_index(drop=True)
        read3 = read3.dropna().reset_index(drop=True)
        clin = pd.read_csv(gp)

        if len(read2s) == 0:
            bio_all = clin[clin["read_3"].isin(read3s)]
        elif len(read3s) == 0:
            bio_all = clin[clin["read_2"].isin(read2s)]
        else:
            bio_all = pd.concat(
                [clin[clin["read_2"].isin(read2s)], clin[clin["read_3"].isin(read3s)]]
            )

        bio_all["value1"] = pd.to_numeric(bio_all["value1"], errors="coerce")
        bio_all["value2"] = pd.to_numeric(bio_all["value2"], errors="coerce")
        bio_all["value3"] = pd.to_numeric(bio_all["value3"], errors="coerce")

        bio_all["value1"] = bio_all["value1"].fillna(bio_all["value2"])
        bio_all["value1"] = bio_all["value1"].fillna(bio_all["value3"])
        bio_all = bio_all[["eid", "event_dt", "value1", "term_description"]]
        bio_all = bio_all.dropna().reset_index(drop=True)
        bio_all = bio_all.sort_values(by=["eid", "event_dt"]).reset_index(drop=True)
        bio_all["event_dt"] = pd.to_datetime(bio_all["event_dt"])
        bio_all = bio_all[bio_all["event_dt"] <= "2015-12-31"].reset_index(drop=True)
        bio_all.rename(columns={"event_dt": "blood_time"}, inplace=True)
        return bio_all


def extract_prescription_ATC(
    code, prescription_path=path2 + "Merged_append4_top100_v4.csv"
):
    pres = pd.read_csv(prescription_path, low_memory=False)
    pres["ATCl"] = pres["ATC"].str[: len(code)]
    pres["event_dt"] = pd.to_datetime(pres["event_dt"])
    pres = pres[pres["event_dt"] <= "2015-12-31"].reset_index(drop=True)
    pres = pres[pres["ATCl"] == code].reset_index(drop=True)
    return pres


def create_intolerance_file(
    pres, bio_all, pc_path="/data/UKBB/PC/PEDMASTER_UNRELATED_WhiteBritish.txt"
):
    df = pd.merge(pres, bio_all, on="eid")
    df["event_dt"] = pd.to_datetime(df["event_dt"])
    df["blood_time"] = pd.to_datetime(df["blood_time"])
    df["days"] = (df["blood_time"] - df["event_dt"]).dt.days
    dfs = df[["eid", "days", "value1", "blood_time", "event_dt"]]
    dfp = dfs[(dfs["days"] >= 0) & (dfs["days"] < 90)]
    df_max = dfp.loc[dfp.groupby("eid").value1.idxmax()].reset_index(drop=True)

    pc = pd.read_csv(pc_path, sep="\t")
    pc = pc.drop(["IID", "birthYear", "X153", "X193", "X365", "X411"], axis=1)
    df_output_count = pd.merge(pc, dfp, left_on="FID", right_on="eid")
    df_output_count = df_output_count.drop(["eid"], axis=1)
    df_output_count = df_output_count.sort_values(by=["FID"]).reset_index(drop=True)
    df_output_max = pd.merge(pc, df_max, left_on="FID", right_on="eid")
    df_output_max = df_output_max.drop(["eid"], axis=1)
    df_output_max = df_output_max.sort_values(by=["FID"]).reset_index(drop=True)

    return df_output_count, df_output_max


def check_intolerance(df_output):
    df = df_output
    filtered_data = df[
        ((df["Sex"] == 1) & (df["value1"] >= 180))
        | ((df["Sex"] == 2) & (df["value1"] >= 120))
    ]
    count = len(filtered_data)
    print(f"statin intolerance case: {count}")
    return filtered_data


def check_SRM(df_output):
    df = df_output
    filtered_data = df[
        ((df["Sex"] == 1) & (df["value1"] >= 4 * 180))
        | ((df["Sex"] == 2) & (df["value1"] >= 4 * 120))
    ]
    count = len(filtered_data)
    print(f"statin SRM case: {count}")
    return filtered_data


def check_SRSR(df_output):
    df = df_output
    filtered_data = df[
        ((df["Sex"] == 1) & (df["value1"] >= 10 * 180))
        | ((df["Sex"] == 2) & (df["value1"] >= 10 * 120))
    ]
    count = len(filtered_data)
    print(f"statin SRSR case: {count}")
    return filtered_data


read2 = ["44H4.", "44HE.", "44HG.", "44HL."]
read3 = ["44H4.", "X80Dl", "X80Dn", "X80N4", "XE28Y", "XaES5", "XaESB", "XaLWP", "XaX4A", "XaY88"]
ck = extract_GP_Biomarker(read2, read3)
null_idx = ck[ck["value1"].isnull()].index
ck.drop(null_idx, inplace=True)
# check the distribution of the biomarker values when required
# vals = ck['value1'].values
# plt.scatter(np.arange(vals.shape[0]), vals, s = 5)

statin_all = extract_prescription_ATC("C10AA")
# If phenotype-type is 'Treatment-modification phenotype':
# Get the statin users who have taken statin for at least 2 times
eid_counts = statin_all.groupby("eid").size()
statin_users_index = eid_counts[eid_counts >= 2].index.tolist()
statin_users = statin_all[statin_all["eid"].isin(statin_users_index)]

# Generic drug name and corresponding word list
generics = {
    "simvastatin": ["simvastatin", "zocor", "simvador"],
    "lovastatin": ["lovastatin"],
    "pravastatin": ["pravastatin", "lipostat"],
    "fluvastatin": ["fluvastatin", "lescol"],
    "atorvastatin": ["atorvastatin", "lipitor"],
    "cerivastatin": ["cerivastatin", "lipobay"],
    "rosuvastatin": ["rosuvastatin", "crestor"],
    "pitavastatin": ["pitavastatin"],
}

statin_users.loc[:, "generic"] = statin_users["mapping_term"]
for generic, keywords in generics.items():
    for keyword in keywords:
        statin_users.loc[
            statin_users["mapping_term"].str.contains(keyword), "generic"
        ] = generic

switching_3 = statin_users.sort_values(by=["eid", "event_dt", "generic"])
switching_3["generic_after"] = switching_3.groupby("eid")["generic"].shift(-1)
switching_3["count"] = 0
switching_year = 2012
switching_3.loc[
    (~switching_3["generic_after"].isna())
    & (switching_3["generic"] != switching_3["generic_after"]),
    "count",
] += 1
switching_3.loc[
    (switching_3["generic"] == "simvastatin")
    & (switching_3["event_dt"].dt.year >= switching_year)
    & (switching_3["generic_after"] == "atorvastatin"),
    "count",
] = 0
switching_3_eid = (
    switching_3.groupby("eid").filter(lambda x: x["count"].sum() >= 3)["eid"].unique()
)
stop_eids = (
    statin_users.groupby("eid")
    .filter(
        lambda group: (pd.Timestamp("2015-12-31") - group["event_dt"].max()).days >= 270
    )["eid"]
    .unique()
)
combined_eids = list(set(stop_eids.tolist() + switching_3_eid.tolist()))

treatment_count, treatment_max = create_intolerance_file(statin_users, ck)
# ck_only_count, ck_only_max = create_intolerance_file(statin_all, ck)

GSI = check_intolerance(treatment_max)
GSI_case = GSI[GSI["FID"].isin(combined_eids)]
population = treatment_max.copy()
population["GSI"] = 0
population.loc[GSI.index, "GSI"] = np.nan
pheno = GSI_case["FID"].tolist()
population.loc[population["FID"].isin(pheno), "GSI"] = 1
population = population[["FID", "Sex", "Age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "GSI"]]
population.rename(columns={"FID": "IID"}, inplace=True)
population.to_csv("gwas_input_GSI.csv", index=False)


# =====================For sensitivity analysis=========================
PATH_PRES = path2 + "Merged_append4_top100_v4.csv"
DATE_CUTOFF = "2015-12-31"
pres = pd.read_csv(PATH_PRES, low_memory=False)
pres["eid"] = pres["eid"].astype(str)
pres["ATC3"] = pres["ATC"].str[:3]
pres["ATC5"] = pres["ATC"].str[:5]
pres["event_dt"] = pd.to_datetime(pres["event_dt"], errors="coerce")
pres = pres[pres["event_dt"] <= DATE_CUTOFF].reset_index(drop=True)

# 1. HTN  (Rx ≥1)
#    C02   — Antihypertensives (centrally acting, alpha-blockers, vasodilators)
#    C03AA — Thiazide diuretics
#    C03BA — Thiazide-like diuretics
#    C08CA — Dihydropyridine CCBs
#    C09   — ACEi, ARBs, combinations, renin inhibitors
#
#    Excluded: C07 (beta-blockers), C03CA (loop diuretics),
#              C03DA (aldosterone antagonists), C08DA/DB (non-DHP CCBs)
# =============================================================
htn_rx = pres[
    (pres["ATC3"] == "C02")
    | (pres["ATC5"] == "C03AA")
    | (pres["ATC5"] == "C03BA")
    | (pres["ATC5"] == "C08CA")
    | (pres["ATC3"] == "C09")
]
htn_eids = set(htn_rx["eid"].unique())

# 2. DM  (Rx ≥1)
#    A10 — All drugs used in diabetes
# =============================================================
dm_rx = pres[pres["ATC3"] == "A10"]
dm_eids = set(dm_rx["eid"].unique())

# 3. Non-statin LLT  (Rx ≥1)
#    C10 — All lipid modifying agents
#    Excluding:
#      C10AA — Statins (plain)
#      C10BX — Statin + non-lipid drug combos (aspirin, amlodipine, etc.)
#    Keeping:
#      C10BA — Statin + lipid agent combos (contains non-statin LLT exposure)
# =============================================================
llt = pres[
    (pres["ATC3"] == "C10")
    & (pres["ATC5"] != "C10AA")
    & (pres["ATC5"] != "C10BX")  # statins plain  # statin + non-lipid drug combos
]
llt_eids = set(llt["eid"].unique())

# 4. PPI  (Rx ≥2)
#    A02BC — Proton pump inhibitors (plain)
# =============================================================
ppi = pres[pres["ATC5"] == "A02BC"]
ppi_counts = ppi.groupby("eid").size()
ppi_eids = set(ppi_counts[ppi_counts >= 2].index)

# Build single covariate file
all_eids = sorted(set(pres["eid"].unique()))
flags = pd.DataFrame({"eid": all_eids})
flags["htn"] = flags["eid"].isin(htn_eids).astype(int)
flags["dm"] = flags["eid"].isin(dm_eids).astype(int)
flags["non_statin_llt"] = flags["eid"].isin(llt_eids).astype(int)
flags["ppi"] = flags["eid"].isin(ppi_eids).astype(int)

flags.to_csv("covariate_flags.csv", index=False)
print(f"Date cutoff: {DATE_CUTOFF}")
print(f"HTN (Rx≥1):             {flags['htn'].sum():,}")
print(f"DM  (Rx≥1):             {flags['dm'].sum():,}")
print(f"Non-statin LLT (Rx≥1):  {flags['non_statin_llt'].sum():,}")
print(f"PPI (Rx≥2):             {flags['ppi'].sum():,}")
print(f"Total eids:             {len(flags):,}")


flags = flags.rename(columns={"eid": "IID"})
flags["IID"] = flags["IID"].astype(int)
sensitivity_input = population.merge(flags, on="IID", how="left")
sensitivity_input[["htn", "dm", "non_statin_llt", "ppi"]] = (
    sensitivity_input[["htn", "dm", "non_statin_llt", "ppi"]].fillna(0).astype(int)
)
sensitivity_input.to_csv("sensitivity_input_GSI.csv", index=False)
