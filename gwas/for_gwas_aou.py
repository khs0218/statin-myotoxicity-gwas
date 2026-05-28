import os
import numpy as np
import pandas as pd

CDR = os.environ["WORKSPACE_CDR"]
COHORT_QUERY = f"SELECT person_id FROM `{CDR}.person`"
if not os.path.exists("gwas_common_files/demographics.parquet.gz"):
    demographics_sql = (
        """
        SELECT
            person.person_id,
            person.gender_concept_id,
            p_gender_concept.concept_name as gender,
            person.birth_datetime as date_of_birth,
            person.race_concept_id,
            p_race_concept.concept_name as race,
            person.ethnicity_concept_id,
            p_ethnicity_concept.concept_name as ethnicity,
            person.sex_at_birth_concept_id,
            p_sex_at_birth_concept.concept_name as sex_at_birth 
        FROM
            `"""
        + os.environ["WORKSPACE_CDR"]
        + """.person` person 
        LEFT JOIN
            `"""
        + os.environ["WORKSPACE_CDR"]
        + """.concept` p_gender_concept 
                ON person.gender_concept_id = p_gender_concept.concept_id 
        LEFT JOIN
            `"""
        + os.environ["WORKSPACE_CDR"]
        + """.concept` p_race_concept 
                ON person.race_concept_id = p_race_concept.concept_id 
        LEFT JOIN
            `"""
        + os.environ["WORKSPACE_CDR"]
        + """.concept` p_ethnicity_concept 
                ON person.ethnicity_concept_id = p_ethnicity_concept.concept_id 
        LEFT JOIN
            `"""
        + os.environ["WORKSPACE_CDR"]
        + """.concept` p_sex_at_birth_concept 
                ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id  
        WHERE
            person.PERSON_ID IN (
                SELECT
                    distinct person_id  
                FROM
                    `"""
        + os.environ["WORKSPACE_CDR"]
        + """.cb_search_person` cb_search_person  
                WHERE
                    cb_search_person.person_id IN (
                        SELECT
                            person_id 
                        FROM
                            `"""
        + os.environ["WORKSPACE_CDR"]
        + """.cb_search_person` p 
                        WHERE
                            has_whole_genome_variant = 1 
                    ) 
                )"""
    )

    demographics_df = pd.read_gbq(
        demographics_sql,
        dialect="standard",
        use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
        progress_bar_type="tqdm_notebook",
    )

    demographics_df = demographics_df[["person_id", "gender", "date_of_birth", "race", "ethnicity", "sex_at_birth"]]
    for col in ["gender", "race", "ethnicity", "sex_at_birth"]:
        demographics_df[col] = demographics_df[col].astype("category")

else:
    demographics_df = pd.read_parquet("gwas_common_files/demographics.parquet.gz")

ances = pd.read_csv("ancestry_preds.tsv", sep="\t")
pca = ances["pca_features"].tolist()
pca = [x.strip("][").split(", ") for x in pca]
pcs_df = pd.DataFrame(pca, columns=["PC" + str(x) for x in range(1, 17)])
cov = ances[["research_id", "ancestry_pred"]]
cov = pd.concat([cov, pcs_df], axis=1)
cov = cov.merge(demographics_df, left_on="research_id", right_on="person_id")
cov = cov.rename(columns={"research_id": "IID", "ancestry_pred": "ancestry"})
cov["SEX"] = [
    "0" if m == "Male" else ("1" if m == "Female" else "") for m in cov["sex_at_birth"]
]
cov_eur = cov[cov["ancestry"] == "eur"]
cov_eur = cov_eur[["IID", "SEX", "AGE", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"]]
sex_age = cov_eur[["IID", "SEX", "AGE"]]
sex_age = sex_age.rename(columns={"IID": "person_id"})


CK_measurement = pd.io.gbq.read_gbq(
    f"""
    SELECT person_id, 
    measurement_concept_id, 
    measurement_date, 
    measurement_type_concept_id,
    operator_concept_id,
    value_as_number,
    value_as_concept_id, 
    unit_concept_id 
    FROM `{CDR}.measurement` 
    WHERE measurement_concept_id IN (3015531, 3007220, 3029790, 3016913, 3008994, 3016070, 3016913) """,
    dialect="standard",
    progress_bar_type="tqdm_notebook",
)

CK_measurement_tmp = CK_measurement
CK_measurement_tmp[["value_as_concept_id"]] = CK_measurement_tmp[
    ["value_as_concept_id"]
].fillna(9999999)
nums = CK_measurement_tmp.loc[
    (CK_measurement_tmp["value_as_number"].notnull())
    & (CK_measurement_tmp["value_as_number"] > 0)
]["value_as_number"].tolist()
CK_measurement_tmp.loc[
    (CK_measurement_tmp["value_as_number"].notnull())
    & (CK_measurement_tmp["value_as_number"] > 0),
    "value_final",
] = nums
CK_measurement_filtered = CK_measurement_tmp.loc[
    (CK_measurement_tmp["value_final"].notnull())
    & (CK_measurement_tmp["value_final"] > 0)
].reset_index()

prescription = pd.io.gbq.read_gbq(
    f"""
    SELECT *
    FROM {CDR}.drug_exposure 
    INNER JOIN {CDR}.concept_ancestor ON drug_concept_id = descendant_concept_id
    INNER JOIN {CDR}.concept c ON ancestor_concept_id = c.concept_id
    WHERE c.concept_id = 21601855;"""
)

id_counts = prescription.groupby("person_id").size()
statin_users_index = id_counts[id_counts >= 2].index.tolist()
statin_users = prescription[prescription["person_id"].isin(statin_users_index)]
statin_users["drug_exposure_end_date"] = statin_users["drug_exposure_end_date"].fillna(
    statin_users["drug_exposure_start_date"]
)
statin_users["drug_exposure_start_date"] = pd.to_datetime(
    statin_users["drug_exposure_start_date"]
)
statin_users["drug_exposure_end_date"] = pd.to_datetime(
    statin_users["drug_exposure_end_date"]
)

for_RxNorm = np.unique(prescription["drug_concept_id"])
for_RxNorm_str = ",".join(map(str, for_RxNorm))
RxNorm = pd.io.gbq.read_gbq(
    f"""
    SELECT *
    FROM {CDR}.concept c 
    WHERE c.concept_id in ({for_RxNorm_str});"""
)

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

RxNorm.loc[:, "generic"] = RxNorm["concept_name"]

for generic, keywords in generics.items():
    for keyword in keywords:
        RxNorm.loc[RxNorm["concept_name"].str.contains(keyword), "generic"] = generic

merge = pd.merge(statin_users, RxNorm, left_on="drug_concept_id", right_on="concept_id")
selected_columns = [
    "drug_exposure_id",
    "person_id",
    "drug_concept_id",
    "drug_exposure_start_date",
    "drug_exposure_end_date",
    "drug_source_value",
    "drug_source_concept_id",
    "ancestor_concept_id",
    "concept_id",
    "concept_name",
    "generic",
]

prescription_selected = prescription[selected_columns]
merge = merge[selected_columns]
switching_3 = merge.sort_values(
    by=["person_id", "drug_exposure_start_date", "drug_exposure_end_date", "generic"]
)
switching_3["generic_after"] = switching_3.groupby("person_id")["generic"].shift(-1)
switching_year = 2012
switching_3["count"] = 0
switching_3.loc[
    (~switching_3["generic_after"].isna())
    & (switching_3["generic"] != switching_3["generic_after"]),
    "count",
] += 1
switching_3.loc[
    (switching_3["generic"] == "simvastatin")
    & (switching_3["drug_exposure_start_date"].dt.year >= switching_year)
    & (switching_3["generic_after"] == "atorvastatin"),
    "count",
] = 0

switching_3_id = (
    switching_3.groupby("person_id")
    .filter(lambda x: x["count"].sum() >= 3)["person_id"]
    .unique()
)
latest_end_date = statin_users["drug_exposure_end_date"].max()
statin_users = statin_users.sort_values(
    by=["drug_exposure_end_date", "drug_exposure_start_date"]
)
stop_eids = (
    statin_users.groupby("person_id")
    .filter(
        lambda group: (
            pd.Timestamp(latest_end_date) - group["drug_exposure_end_date"].max()
        ).days
        >= 270
    )["person_id"]
    .unique()
)
combined_ids = list(set(stop_eids.tolist() + switching_3_id.tolist()))


def create_intolerance_file(pres, bio_all, sex_age):
    df = pd.merge(pres, bio_all, on="person_id")
    df = pd.merge(df, sex_age, on="person_id")
    df["measurement_date"] = pd.to_datetime(df["measurement_date"])
    df["drug_exposure_start_date"] = pd.to_datetime(df["drug_exposure_start_date"])
    df["drug_exposure_end_date"] = pd.to_datetime(df["drug_exposure_end_date"])
    df["days_1"] = (df["measurement_date"] - df["drug_exposure_start_date"]).dt.days
    df["days_2"] = (df["measurement_date"] - df["drug_exposure_end_date"]).dt.days

    # Get data containing information within three months after prescription
    dfs = df[
        [
            "person_id",
            "SEX",
            "AGE",
            "drug_source_concept_id",
            "value_final",
            "days_1",
            "days_2",
            "drug_exposure_start_date",
            "drug_exposure_end_date",
            "measurement_date",
        ]
    ]
    dfp = dfs[(dfs["days_1"] >= 0) & (dfs["days_2"] <= 90)]
    df_max = dfp.loc[dfp.groupby("person_id").value_final.idxmax()].reset_index(
        drop=True
    )

    return dfp, df_max


def check_intolerance(df_output):
    df = df_output
    filtered_data = df[
        ((df["SEX"] == "0") & (df["value_final"] > 180))
        | ((df["SEX"] == "1") & (df["value_final"] > 120))
    ]
    count = len(filtered_data)
    print(f"statin intolerance case: {count}")

    return filtered_data


def check_SRM(df_output):
    df = df_output
    filtered_data = df[
        ((df["SEX"] == "0") & (df["value_final"] > 4 * 180))
        | ((df["SEX"] == "1") & (df["value_final"] > 4 * 120))
    ]
    count = len(filtered_data)
    print(f"statin SRM case: {count}")

    return filtered_data


def check_SRSR(df_output):
    df = df_output
    filtered_data = df[
        ((df["SEX"] == "0") & (df["value_final"] > 10 * 180))
        | ((df["SEX"] == "1") & (df["value_final"] > 10 * 120))
    ]
    count = len(filtered_data)
    print(f"statin SRSR case: {count}")

    return filtered_data


treatment_count, treatment_max = create_intolerance_file(
    statin_users, CK_measurement_filtered, sex_age
)
# ck_only_count, ck_only_max = create_intolerance_file(prescription_selected, CK_measurement_filtered, sex_age)

GSI = check_intolerance(treatment_max)
GSI_case = GSI[GSI["person_id"].isin(combined_ids)]
population = treatment_max.copy()
population["GSI"] = 0
pheno = GSI["person_id"].tolist()
population.loc[population["person_id"].isin(pheno), "GSI"] = np.nan
pheno = GSI_case["person_id"].tolist()
population.loc[population["person_id"].isin(pheno), "GSI"] = 1
population.drop(columns=["SEX", "AGE"], inplace=True)
merge_pc = pd.merge(cov_eur, population, left_on="IID", right_on="person_id")
merge_pc = merge_pc[["IID", "SEX", "AGE", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "GSI"]]
merge_pc.to_csv("gwas_input_GSI.csv", index=False)


# =====================For sensitivity analysis=========================
"""
All of Us — Covariate Flags for GWAS Sensitivity Analysis
==========================================================
All covariates defined by medication records (ATC-based):
  HTN:              Rx ≥1  (C02, C03AA, C03BA, C08CA, C09)
  DM:               Rx ≥1  (A10)
  Non-statin LLT:   Rx ≥1  (C10 excl. C10AA, C10BX)
  PPI:              Rx ≥2  (A02BC)
"""

sql = f"""
WITH
htn_drug_concepts AS (
  SELECT DISTINCT ca.descendant_concept_id AS concept_id
  FROM `{CDR}.concept_ancestor` ca
  JOIN `{CDR}.concept` anc ON ca.ancestor_concept_id = anc.concept_id
  WHERE anc.vocabulary_id = 'ATC'
    AND (
      -- C02 (all, ATC 2nd level)
      (anc.concept_code = 'C02' AND anc.concept_class_id = 'ATC 2nd')
      -- C03AA, C03BA (ATC 4th level)
      OR (anc.concept_code IN ('C03AA', 'C03BA') AND anc.concept_class_id = 'ATC 4th')
      -- C08CA (ATC 4th level)
      OR (anc.concept_code = 'C08CA' AND anc.concept_class_id = 'ATC 4th')
      -- C09 (all, ATC 2nd level)
      OR (anc.concept_code = 'C09' AND anc.concept_class_id = 'ATC 2nd')
    )
),

person_htn AS (
  SELECT DISTINCT person_id
  FROM `{CDR}.drug_exposure`
  WHERE drug_concept_id IN (SELECT concept_id FROM htn_drug_concepts)
),


dm_drug_concepts AS (
  SELECT DISTINCT ca.descendant_concept_id AS concept_id
  FROM `{CDR}.concept_ancestor` ca
  JOIN `{CDR}.concept` anc ON ca.ancestor_concept_id = anc.concept_id
  WHERE anc.vocabulary_id = 'ATC'
    AND anc.concept_code = 'A10'
    AND anc.concept_class_id = 'ATC 2nd'
),

person_dm AS (
  SELECT DISTINCT person_id
  FROM `{CDR}.drug_exposure`
  WHERE drug_concept_id IN (SELECT concept_id FROM dm_drug_concepts)
),


c10_all AS (
  SELECT DISTINCT ca.descendant_concept_id AS concept_id
  FROM `{CDR}.concept_ancestor` ca
  JOIN `{CDR}.concept` anc ON ca.ancestor_concept_id = anc.concept_id
  WHERE anc.vocabulary_id = 'ATC'
    AND anc.concept_code = 'C10'
    AND anc.concept_class_id = 'ATC 2nd'
),

statin_exclude AS (
  SELECT DISTINCT ca.descendant_concept_id AS concept_id
  FROM `{CDR}.concept_ancestor` ca
  JOIN `{CDR}.concept` anc ON ca.ancestor_concept_id = anc.concept_id
  WHERE anc.vocabulary_id = 'ATC'
    AND anc.concept_code IN ('C10AA', 'C10BX')
    AND anc.concept_class_id = 'ATC 4th'
),


llt_concepts AS (
  SELECT concept_id FROM c10_all
  EXCEPT DISTINCT
  SELECT concept_id FROM statin_exclude
),

person_llt AS (
  SELECT DISTINCT person_id
  FROM `{CDR}.drug_exposure`
  WHERE drug_concept_id IN (SELECT concept_id FROM llt_concepts)
),


ppi_concepts AS (
  SELECT DISTINCT ca.descendant_concept_id AS concept_id
  FROM `{CDR}.concept_ancestor` ca
  JOIN `{CDR}.concept` anc ON ca.ancestor_concept_id = anc.concept_id
  WHERE anc.vocabulary_id = 'ATC'
    AND anc.concept_code = 'A02BC'
    AND anc.concept_class_id = 'ATC 4th'
),

person_ppi AS (
  SELECT person_id, COUNT(*) AS n_rx
  FROM `{CDR}.drug_exposure`
  WHERE drug_concept_id IN (SELECT concept_id FROM ppi_concepts)
  GROUP BY person_id
  HAVING COUNT(*) >= 2
)


SELECT
  p.person_id,
  CASE WHEN h.person_id   IS NOT NULL THEN 1 ELSE 0 END AS htn,
  CASE WHEN d.person_id   IS NOT NULL THEN 1 ELSE 0 END AS dm,
  CASE WHEN l.person_id   IS NOT NULL THEN 1 ELSE 0 END AS non_statin_llt,
  CASE WHEN pp.person_id  IS NOT NULL THEN 1 ELSE 0 END AS ppi
FROM `{CDR}.person` p
LEFT JOIN person_htn h  ON p.person_id = h.person_id
LEFT JOIN person_dm  d  ON p.person_id = d.person_id
LEFT JOIN person_llt l  ON p.person_id = l.person_id
LEFT JOIN person_ppi pp ON p.person_id = pp.person_id
"""


flags = pd.read_gbq(sql, dialect="standard")
flags.to_csv("aou_covariate_flags_pres.csv", index=False)

print(f"HTN (Rx≥1):             {flags['htn'].sum():,}")
print(f"DM  (Rx≥1):             {flags['dm'].sum():,}")
print(f"Non-statin LLT (Rx≥1):  {flags['non_statin_llt'].sum():,}")
print(f"PPI (Rx≥2):             {flags['ppi'].sum():,}")
print(f"Total persons:          {len(flags):,}")


flags = flags.rename(columns={"person_id": "IID"})
sensitivity_input = merge_pc.merge(flags, on="IID", how="left")
sensitivity_input[["htn", "dm", "non_statin_llt", "ppi"]] = (
    sensitivity_input[["htn", "dm", "non_statin_llt", "ppi"]].fillna(0).astype(int)
)
sensitivity_input.to_csv("sensitivity_input_GSI.csv", index=False)
