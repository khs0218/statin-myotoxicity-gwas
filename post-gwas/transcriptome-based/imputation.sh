#!/bin/bash
export GWAS_TOOLS=~/MetaXcan/summary-gwas-imputation/src
export DATA=~/MetaXcan/summary-gwas-imputation/data
export OUTPUT=~/MetaXcan/summary-gwas-imputation/output
export GWAS=~/data/


for chr in {1..22}; do
  for sub_batch in {0..9}; do
    python3 ${GWAS_TOOLS}/gwas_summary_imputation.py \
      -by_region_file ${DATA}/eur_ld.bed.gz \
      -gwas_file ${OUTPUT}/harmonized_CK_GSI.txt.gz \
      -parquet_genotype ${DATA}/reference_panel_1000G/chr${chr}.variants.parquet \
      -parquet_genotype_metadata ${DATA}/reference_panel_1000G/variant_metadata.parquet \
      -parsimony 7 \
      -regularization 0.1 \
      -frequency_filter 0.01 \
      -chromosome ${chr} \
      -sub_batch ${sub_batch} \
      -sub_batches 10 \
      --standardise_dosages \
      -output ${OUTPUT}/imputation/UKBB_chr${chr}_sb${sub_batch}.txt.gz
  done
done


python3 ${GWAS_TOOLS}/gwas_summary_imputation_postprocess.py \
  -gwas_file ${OUTPUT}/harmonized_CK_GSI.txt.gz \
  -folder ${OUTPUT}/imputation \
  -pattern "chr.*" \
  -parsimony 7 \
  -output ${OUTPUT}/imputation/imputed_CK_GSI.txt.gz
