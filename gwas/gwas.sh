#!/bin/bash
# Single-variant GWAS with SAIGE (Scalable and Accurate Implementation of
# Generalized mixed model), run on the UK Biobank Research Analysis Platform
# via the DNAnexus Swiss Army Knife app.

# SAIGE step 1: fit the null model
dx run swiss-army-knife \
  -iin=UKB_Main:/Tutorial/example/output/gwas_input_GSI.csv \
  -iin=UKB_Main:/pruned/prune_all.bed \
  -iin=UKB_Main:/pruned/prune_all.bim \
  -iin=UKB_Main:/pruned/prune_all.fam \
  -icmd="step1_fitNULLGLMM.R \
    --plinkFile=prune_all\
    --phenoFile=gwas_input_GSI.csv \
    --phenoCol=GSI \
    --covarColList=Age,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --sampleIDColinphenoFile=IID \
    --traitType=binary \
    --outputPrefix=step1 \
    --IsOverwriteVarianceRatioFile=TRUE \
    --invNormalize=False" \
  -iimage_file="/docker_images/saige_1.1.8.tar.gz" \
  --destination UKB_Main:/Tutorial/example/output \
  --yes \
  --name=GSI_step1 \
  --instance-type=mem3_ssd1_v2_x8


# SAIGE step 2: single-variant association tests per chromosome
for i in {1..22}
do
dx run swiss-army-knife \
  -iin=UKB_Main:/Bulk/Imputation/ukb_c${i}_b0_v3.bgen \
  -iin=UKB_Main:/Bulk/Imputation/ukb_c${i}_b0_v3.bgen.bgi \
  -iin=UKB_Main:/Bulk/Imputation/ukb_c${i}_b0_v3.sample \
  -iin=UKB_Main:/Tutorial/example/output/step1.rda \
  -iin=UKB_Main:/Tutorial/example/output/step1.varianceRatio.txt \
  -icmd="step2_SPAtests.R \
    --bgenFile=ukb_c${i}_b0_v3.bgen \
    --bgenFileIndex=ukb_c${i}_b0_v3.bgen.bgi \
    --sampleFile=ukb_c${i}_b0_v3.sample \
    --chrom=${i}  \
    --AlleleOrder=ref-first \
    --minMAF=0 \
    --minMAC=20 \
    --GMMATmodelFile=step1.rda \
    --varianceRatioFile=step1.varianceRatio.txt \
    --SAIGEOutputFile=step2_chr${i} \
    --LOCO=FALSE \
    --is_fastTest=TRUE" \
  -iimage_file=/docker_images/saige_1.1.8.tar.gz \
  --destination UKB_Main:/Tutorial/example/output \
  --instance-type mem2_ssd2_v2_x4 \
  --name=GSI_step2_chr${i}  \
  --yes
done
