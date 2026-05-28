#!/bin/bash

EQTL_DIR=~/MetaXcan/eqtl/mashr
OUTPUT=~/MetaXcan/summary-gwas-imputation/output
SPX_OUT_DIR=~/MetaXcan/output
GWAS_FILE=${OUTPUT}/imputation/imputed_CK_GSI.txt.gz

for MODEL_DB in ${EQTL_DIR}/*.db; do
  BASE=$(basename "${MODEL_DB}" .db)
  TISSUE=${BASE#mashr_}
  OUT_FILE=${SPX_OUT_DIR}/CK_GSI_${TISSUE}_spredixcan.txt
  COV_FILE=${EQTL_DIR}/${BASE}.txt.gz

  if [[ ! -f "${COV_FILE}" ]]; then
    echo "Warning: covariance file not found for model ${MODEL_DB} (expected ${COV_FILE}), skipping..."
    continue
  fi

  echo "Running SPrediXcan for model: ${BASE}"
  echo "  Model DB : ${MODEL_DB}"
  echo "  Covariance: ${COV_FILE}"
  echo "  Output    : ${OUT_FILE}"

  python3 ~/MetaXcan/software/SPrediXcan.py \
    --model_db_path "${MODEL_DB}" \
    --covariance   "${COV_FILE}" \
    --gwas_file    "${GWAS_FILE}" \
    --snp_column           panel_variant_id \
    --effect_allele_column effect_allele \
    --non_effect_allele_column non_effect_allele \
    --zscore_column        zscore \
    --pvalue_column        pvalue \
    --keep_non_rsid \
    --additional_output \
    --model_db_snp_key varID \
    --throw \
    --output_file "${OUT_FILE}"

  echo "Done: ${BASE}"
  echo "---------------------------------------------"
done