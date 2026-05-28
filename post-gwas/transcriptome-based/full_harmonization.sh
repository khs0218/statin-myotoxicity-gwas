#!/bin/bash
  
export GWAS_TOOLS=~/MetaXcan/summary-gwas-imputation/src
export DATA=~/MetaXcan/summary-gwas-imputation/data
export OUTPUT=~/MetaXcan/summary-gwas-imputation/output
export GWAS=~/data/


# The hg19 -> hg38 liftover step applies to the UK Biobank summary statistics,
python3 ${GWAS_TOOLS}/gwas_parsing.py \
  -gwas_file ${GWAS}/CK_GSI.txt.gz \
  -liftover ${DATA}/liftover/hg19ToHg38.over.chain.gz \
  -snp_reference_metadata ${DATA}/reference_panel_1000G/variant_metadata.txt.gz METADATA \
  -output_column_map CHR chromosome \
  -output_column_map POS position \
  -output_column_map Allele1 non_effect_allele \
  -output_column_map Allele2 effect_allele \
  -output_column_map AF_Allele2 frequency \
  -output_column_map BETA effect_size \
  -output_column_map SE standard_error \
  -output_column_map p.value pvalue \
  --chromosome_format \
  --insert_value sample_size 9994 --insert_value n_cases 2927 \
  -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
  -output ${OUTPUT}/harmonized_CK_GSI.txt.gz
  
  