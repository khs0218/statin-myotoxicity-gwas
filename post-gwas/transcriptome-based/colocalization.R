library(coloc)

df <- read.table("coloc_input_ANO5.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gwas_snp  = df$panel_variant_id
gwas_beta = df$effect_size
gwas_se   = df$standard_error
gwas_maf  = df$frequency
gwas_N    = df$sample_size[1]
gwas_case    = df$n_cases[1]
gwas_s    = gwas_case / gwas_N

eqtl_beta = df$slope
eqtl_se   = df$slope_se
eqtl_maf  = df$maf
eqtl_N <- 619  # GTEx specific tissue sample size

res <- coloc.abf(
    dataset1 = list(snp=gwas_snp, beta=eqtl_beta, varbeta=eqtl_se^2, N=eqtl_N, MAF=eqtl_maf, type="quant"),
    dataset2 = list(snp=gwas_snp, beta=gwas_beta, varbeta=gwas_se^2, N=gwas_N, MAF=gwas_maf, type="cc", s=gwas_s)
)

print(res$summary)
cat("\nPosterior for shared causal variant (H4):\n")
print(res$summary["PP.H4.abf"])