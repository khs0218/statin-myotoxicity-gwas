#!/bin/bash
# Exome-wide gene-based association tests with SAIGE-GENE+, run on the
# UK Biobank Research Analysis Platform via the DNAnexus Swiss Army Knife app.

# SAIGE-GENE+ step 1 reuses the null model fitted in the single-variant GWAS

# SAIGE-GENE+ step 2: gene-based tests per chromosome
for i in {1..22}
do
dx run swiss-army-knife \
  -iin=UKB_Main:/WES_470k/Plink_QC/plink_files/merged_by_chr${i}.norm.bed \
  -iin=UKB_Main:/WES_470k/Plink_QC/plink_files/merged_by_chr${i}.norm.bim \
  -iin=UKB_Main:/WES_470k/Plink_QC/plink_files/merged_by_chr${i}.norm.fam \
  -iin=UKB_Main:/WES_470k/group_files_LOFTEE/loftee_groupfile_chr${i}.withQC.txt \
  -iin=UKB_Main:/Tutorial/example/output/gene+/step1.rda \
  -iin=UKB_Main:/Tutorial/example/output/gene+/step1.varianceRatio.txt \
  -icmd="step2_SPAtests.R \
    --bedFile=merged_by_chr${i}.norm.bed \
    --bimFile=merged_by_chr${i}.norm.bim \
    --famFile=merged_by_chr${i}.norm.fam \
    --chrom=${i}  \
    --AlleleOrder=alt-first \
    --minMAF=0 \
    --minMAC=0.5 \
    --GMMATmodelFile=step1.rda \
    --varianceRatioFile=step1.varianceRatio.txt \
    --groupFile=loftee_groupfile_chr${i}.withQC.txt \
    --annotation_in_groupTest='lof,missense:lof,missense:lof:synonymous' \
    --SAIGEOutputFile=step2_gene+_chr${i} \
    --LOCO=FALSE \
    --is_fastTest=TRUE" \
  -iimage_file=/docker_images/saige_1.1.8.tar.gz \
  --destination UKB_Main:/Tutorial/example/output/gene+ \
  --instance-type mem2_ssd2_v2_x4 \
  --name=step2_gene+_chr${i}  \
  --yes
done
