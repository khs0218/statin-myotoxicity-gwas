"""
UK Biobank — Table 1 Anchored to Max Eligible CK Date
======================================================
DOB/Sex: from participant CSV (year_of_birth, month_of_birth, sex)
BMI:     assessment center (from participant CSV, with blood date)
         + GP clinical records (Read codes, with event_dt)
         → pick nearest to max eligible CK date

Inputs:
  - test_max: from create_intolerance_file(), one row per person
              columns: eid, measurement_date, value_final, SEX
  - covariate_flags: from sensitivity analysis input,
              columns: eid, htn, dm, non_statin_llt, ppi
"""

import os
import pandas as pd

PATH_PARTICIPANT = "/data/participants.csv"
PATH_GP_CLINICAL = "data/GP/RecordLevel/clinical_processed_recordlevel2.csv"
PATH_COVARIATE = "covariate_flags.csv"
PATH_POPULATION = "test_max.csv"


def remove_tz(s):
    s = pd.to_datetime(s)
    if s.dt.tz is not None:
        s = s.dt.tz_localize(None)
    return s


def load_participant(path):
    """Load participant file, extract DOB, sex, BMI, assessment date."""
    cols = ["eid", "sex", "month_of_birth", "year_of_birth", "BMI",
            "blood0", "blood1", "blood2", "blood3", "blood4", "blood5", "blood6"]
    df = pd.read_csv(path, usecols=cols, low_memory=False)
    df["eid"] = df["eid"].astype(str)
    df["date_of_birth"] = pd.to_datetime(
        df["year_of_birth"].astype(int).astype(str)
        + "-"
        + df["month_of_birth"].astype(int).astype(str).str.zfill(2)
        + "-15",
        errors="coerce",
    )

    blood_cols = ["blood0", "blood1", "blood2", "blood3", "blood4", "blood5", "blood6"]
    for col in blood_cols:
        df[col] = pd.to_datetime(df[col], errors="coerce")
        df[col] = df[col].where(df[col] != pd.Timestamp("1900-01-01"))

    df["assessment_date"] = df[blood_cols].min(axis=1)
    df["bmi_ac"] = pd.to_numeric(df["BMI"], errors="coerce")

    return df[["eid", "sex", "date_of_birth", "assessment_date", "bmi_ac"]]


BMI_READ2 = ["22K."]
BMI_READ3 = ["X76CO", "XaJJH"]


def extract_bmi_gp(path_gp_clinical):
    """Extract BMI from GP clinical records with dates."""
    gp = pd.read_csv(path_gp_clinical, low_memory=False)
    bmi_r2 = gp[gp["read_2"].isin(BMI_READ2)]
    bmi_r3 = gp[gp["read_3"].isin(BMI_READ3)]
    bmi = pd.concat([bmi_r2, bmi_r3])

    bmi["value1"] = pd.to_numeric(bmi["value1"], errors="coerce")
    bmi["value2"] = pd.to_numeric(bmi["value2"], errors="coerce")
    bmi["value3"] = pd.to_numeric(bmi["value3"], errors="coerce")
    bmi["value1"] = bmi["value1"].fillna(bmi["value2"]).fillna(bmi["value3"])

    bmi = bmi[["eid", "event_dt", "value1"]].dropna()
    bmi.columns = ["eid", "measurement_date", "bmi_value"]
    bmi["eid"] = bmi["eid"].astype(str)
    bmi["measurement_date"] = pd.to_datetime(bmi["measurement_date"], errors="coerce")
    bmi = bmi[(bmi["bmi_value"] >= 10) & (bmi["bmi_value"] <= 80)].dropna()

    return bmi


def build_table1_ukb(
    path_population, path_participant, path_gp_clinical, path_covariate
):
    population = pd.read_csv(path_population, low_memory=False)
    df = population[["FID", "blood_time", "value1", "Sex"]].copy()
    df.rename(
        columns={
            "FID": "eid",
            "blood_time": "anchor_date",
            "value1": "ck_max",
            "Sex": "sex_code",
        },
        inplace=True,
    )
    df["eid"] = df["eid"].astype(str)
    df["anchor_date"] = pd.to_datetime(df["anchor_date"])
    print("Loading participant data (DOB, sex, assessment BMI)...")
    participant = load_participant(path_participant)
    df = df.merge(participant, on="eid", how="left")

    df["anchor_date"] = remove_tz(df["anchor_date"])
    df["date_of_birth"] = remove_tz(df["date_of_birth"])
    df["age_at_ck"] = (df["anchor_date"] - df["date_of_birth"]).dt.days / 365.25

    print("Extracting GP BMI records...")
    bmi_gp = extract_bmi_gp(path_gp_clinical)

    bmi_ac = (
        df[["eid", "assessment_date", "bmi_ac"]]
        .dropna(subset=["bmi_ac", "assessment_date"])
        .copy()
    )
    bmi_ac.columns = ["eid", "measurement_date", "bmi_value"]
    bmi_all = pd.concat([bmi_gp, bmi_ac], ignore_index=True)

    print("Finding BMI nearest to CK date...")
    bmi_merged = bmi_all.merge(df[["eid", "anchor_date"]], on="eid", how="inner")
    bmi_merged["measurement_date"] = remove_tz(bmi_merged["measurement_date"])
    bmi_merged["anchor_date"] = remove_tz(bmi_merged["anchor_date"])
    bmi_merged["days_diff"] = (
        bmi_merged["measurement_date"] - bmi_merged["anchor_date"]
    ).dt.days.abs()

    bmi_nearest = bmi_merged.loc[bmi_merged.groupby("eid")["days_diff"].idxmin()]
    bmi_nearest = bmi_nearest[["eid", "bmi_value", "days_diff"]].rename(
        columns={"bmi_value": "bmi_nearest_ck", "days_diff": "bmi_days_from_ck"}
    )
    df = df.merge(bmi_nearest, on="eid", how="left")

    print("Merging covariate flags...")
    flags = pd.read_csv(path_covariate, dtype=str)
    flags["eid"] = flags["eid"].astype(str)
    for col in ["htn", "dm", "non_statin_llt", "ppi"]:
        flags[col] = pd.to_numeric(flags[col], errors="coerce").fillna(0).astype(int)
    df = df.merge(flags, on="eid", how="left")
    for col in ["htn", "dm", "non_statin_llt", "ppi"]:
        df[col] = df[col].fillna(0).astype(int)

    df["sex"] = df["sex"].astype(str)
    df["is_male"] = (df["sex_code"] == 1).astype(int)

    print(f"Table 1 dataframe built: {len(df):,} participants")
    return df


def summarize_table1(df, cohort_name="UK Biobank"):
    n = len(df)
    males = df[df["is_male"] == 1]
    females = df[df["is_male"] == 0]

    def med_iqr(s):
        s = s.dropna()
        if len(s) == 0:
            return "N/A"
        return f"{s.median():.1f} ({s.quantile(0.25):.1f}, {s.quantile(0.75):.1f})"

    def n_pct(count, total=n):
        return f"{count} ({100*count/total:.1f}%)"

    print(f"\n{'='*65}")
    print(f"Table 1: Study population — {cohort_name} (n={n:,})")
    print(f"{'='*65}")

    print(f"\nAge at CK ascertainment (years)")
    print(f"  Male:   {med_iqr(males['age_at_ck'])}")
    print(f"  Female: {med_iqr(females['age_at_ck'])}")
    print(f"  All:    {med_iqr(df['age_at_ck'])}")

    print(f"\nSex, n(%)")
    print(f"  Male:   {n_pct(df['is_male'].sum())}")
    print(f"  Female: {n_pct((1 - df['is_male']).sum())}")

    print(f"\nBMI nearest to CK ascertainment (kg/m²)")
    print(f"  Male:   {med_iqr(males['bmi_nearest_ck'])}")
    print(f"  Female: {med_iqr(females['bmi_nearest_ck'])}")
    print(f"  All:    {med_iqr(df['bmi_nearest_ck'])}")
    bmi_miss = df["bmi_nearest_ck"].isna().sum()
    if bmi_miss > 0:
        print(f"  Missing: {bmi_miss} ({100*bmi_miss/n:.1f}%)")
    bmi_med_days = df["bmi_days_from_ck"].dropna().median()
    print(f"  Median days from CK: {bmi_med_days:.0f}")

    print(f"\nMaximum eligible CK level (IU/L)")
    print(f"  Male:   {med_iqr(males['ck_max'])}")
    print(f"  Female: {med_iqr(females['ck_max'])}")
    print(f"  All:    {med_iqr(df['ck_max'])}")

    print(
        f"\nEver use of non-statin lipid-lowering medication: {n_pct(df['non_statin_llt'].sum())}"
    )
    print(f"Ever use of blood pressure-lowering medication:   {n_pct(df['htn'].sum())}")
    print(f"Ever use of diabetes medication:                  {n_pct(df['dm'].sum())}")
    print(f"Ever use of proton pump inhibitors (≥2 Rx):       {n_pct(df['ppi'].sum())}")


table1_df = build_table1_ukb(
    path_population=PATH_POPULATION,
    path_participant=PATH_PARTICIPANT,
    path_gp_clinical=PATH_GP_CLINICAL,
    path_covariate=PATH_COVARIATE,
)
summarize_table1(table1_df, "UK Biobank (UKB)")
table1_df.to_csv("ukb_table1_anchored.csv", index=False)


"""
All of Us — Table 1 Anchored to Max Eligible CK Date
=====================================================
Uses sex, age from AoU person table, BMI from measurement table.

Inputs:
  - test_max: from create_intolerance_file(), one row per person
              columns: person_id, measurement_date, value_final, SEX
  - covariate_flags: from sensitivity analysis input,
              columns: person_id, htn, dm, non_statin_llt, ppi

Run in All of Us Researcher Workbench notebook.
"""

DATASET = os.environ["WORKSPACE_CDR"]


def extract_demographics_aou():
    sql = f"""
    SELECT
        person.person_id,
        person.gender_concept_id,
        p_gender_concept.concept_name AS gender,
        person.birth_datetime AS date_of_birth,
        person.sex_at_birth_concept_id,
        p_sex_at_birth_concept.concept_name AS sex_at_birth
    FROM `{DATASET}.person` person
    LEFT JOIN `{DATASET}.concept` p_gender_concept
        ON person.gender_concept_id = p_gender_concept.concept_id
    LEFT JOIN `{DATASET}.concept` p_sex_at_birth_concept
        ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id
    """
    demo = pd.read_gbq(sql, dialect="standard")
    demo["date_of_birth"] = pd.to_datetime(demo["date_of_birth"]).dt.tz_localize(None)
    for col in ["gender", "sex_at_birth"]:
        demo[col] = demo[col].astype("category")
    return demo


def extract_bmi_aou():
    sql = f"""
    SELECT
      person_id,
      measurement_date,
      value_as_number AS bmi_value
    FROM `{DATASET}.measurement`
    WHERE measurement_concept_id = 3038553
      AND value_as_number IS NOT NULL
      AND value_as_number BETWEEN 10 AND 80
    """
    return pd.read_gbq(sql, dialect="standard")


def build_table1_aou(test_max, covariate_flags):
    df = test_max[["person_id", "measurement_date", "value_final", "SEX"]].copy()
    df.rename(
        columns={"measurement_date": "anchor_date", "value_final": "ck_max"},
        inplace=True,
    )
    df["anchor_date"] = pd.to_datetime(df["anchor_date"])
    df["SEX"] = df["SEX"].astype(str)
    print("Extracting demographics (DOB, sex_at_birth, gender)...")
    demo = extract_demographics_aou()
    df = df.merge(
        demo[["person_id", "date_of_birth", "sex_at_birth", "gender"]],
        on="person_id",
        how="left",
    )

    df["age_at_ck"] = (df["anchor_date"] - df["date_of_birth"]).dt.days / 365.25

    print("Extracting all BMI measurements...")
    bmi_all = extract_bmi_aou()
    bmi_all["measurement_date"] = pd.to_datetime(bmi_all["measurement_date"])
    bmi_merged = bmi_all.merge(
        df[["person_id", "anchor_date"]], on="person_id", how="inner"
    )
    bmi_merged["days_diff"] = (
        bmi_merged["measurement_date"] - bmi_merged["anchor_date"]
    ).dt.days.abs()
    bmi_nearest = bmi_merged.loc[bmi_merged.groupby("person_id")["days_diff"].idxmin()]
    bmi_nearest = bmi_nearest[["person_id", "bmi_value", "days_diff"]].rename(
        columns={"bmi_value": "bmi_nearest_ck", "days_diff": "bmi_days_from_ck"}
    )
    df = df.merge(bmi_nearest, on="person_id", how="left")

    print("Merging covariate flags...")
    flags = covariate_flags.copy()
    for col in ["htn", "dm", "non_statin_llt", "ppi"]:
        if col in flags.columns:
            flags[col] = (
                pd.to_numeric(flags[col], errors="coerce").fillna(0).astype(int)
            )
    df = df.merge(flags, on="person_id", how="left")
    for col in ["htn", "dm", "non_statin_llt", "ppi"]:
        df[col] = df[col].fillna(0).astype(int)

    print(f"Table 1 dataframe built: {len(df):,} participants")
    return df


def summarize_table1(df, cohort_name="All of Us"):
    n = len(df)

    def med_iqr(s):
        s = s.dropna()
        if len(s) == 0:
            return "N/A"
        return f"{s.median():.1f} ({s.quantile(0.25):.1f}, {s.quantile(0.75):.1f})"

    def n_pct(count, total=n):
        return f"{count} ({100*count/total:.1f}%)"

    print(f"\n{'='*65}")
    print(f"Table 1: Study population — {cohort_name} (n={n:,})")
    print(f"{'='*65}")

    print(f"\nSex at birth, n(%)")
    for val in df["sex_at_birth"].dropna().unique():
        c = (df["sex_at_birth"] == val).sum()
        print(f"  {val}: {n_pct(c)}")

    print(f"\nGender identity, n(%)")
    for val in df["gender"].dropna().unique():
        c = (df["gender"] == val).sum()
        print(f"  {val}: {n_pct(c)}")

    males = df[df["sex_at_birth"] == "Male"]
    females = df[df["sex_at_birth"] == "Female"]

    print(f"\nAge at CK ascertainment (years)")
    print(f"  Male:   {med_iqr(males['age_at_ck'])}")
    print(f"  Female: {med_iqr(females['age_at_ck'])}")
    print(f"  All:    {med_iqr(df['age_at_ck'])}")

    print(f"\nBMI nearest to CK ascertainment (kg/m²)")
    print(f"  Male:   {med_iqr(males['bmi_nearest_ck'])}")
    print(f"  Female: {med_iqr(females['bmi_nearest_ck'])}")
    print(f"  All:    {med_iqr(df['bmi_nearest_ck'])}")
    bmi_miss = df["bmi_nearest_ck"].isna().sum()
    if bmi_miss > 0:
        print(f"  Missing: {bmi_miss} ({100*bmi_miss/n:.1f}%)")
    bmi_med_days = df["bmi_days_from_ck"].dropna().median()
    print(f"  Median days from CK: {bmi_med_days:.0f}")

    print(f"\nMaximum eligible CK level (IU/L)")
    print(f"  Male:   {med_iqr(males['ck_max'])}")
    print(f"  Female: {med_iqr(females['ck_max'])}")
    print(f"  All:    {med_iqr(df['ck_max'])}")

    print(
        f"\nEver use of non-statin lipid-lowering medication: {n_pct(df['non_statin_llt'].sum())}"
    )
    print(f"Ever use of blood pressure-lowering medication:   {n_pct(df['htn'].sum())}")
    print(f"Ever use of diabetes medication:                  {n_pct(df['dm'].sum())}")
    print(f"Ever use of proton pump inhibitors (≥2 Rx):       {n_pct(df['ppi'].sum())}")


table1_df = build_table1_aou(
    test_max=test_max,  # from create_intolerance_file()
    covariate_flags=pd.read_csv("aou_covariate_flags.csv"),
)
summarize_table1(table1_df, "All of Us (AoU)")
table1_df.to_csv("aou_table1_anchored.csv", index=False)
