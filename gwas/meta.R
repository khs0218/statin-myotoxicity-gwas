library(data.table)
library(dplyr)
library(VGAM)

# TODO: set the per-cohort sample sizes
N1 <- NA  # UK Biobank sample size
N2 <- NA  # All of Us sample size

ukb<-fread(paste0('GSI_UKB.txt'))
aou<-fread(paste0('GSI_AoU.txt'))
colnames(ukb) <- c('CHR', 'BP', 'ALLELE0','ALLELE1', 'BETA','SE', 'pval')
colnames(aou) <- c('CHR', 'BP', 'ALLELE0','ALLELE1', 'BETA','SE', 'pval')
ukb$BP = as.numeric(ukb$BP)
aou$BP <- as.numeric(aou$BP)

df <-inner_join(ukb,aou,by=c('CHR','BP'))
colSums(is.na(df)) # no NA checked
df$beta_meta <- df$BETA.x * (((1 / df$SE.x)^2) / ((1 / df$SE.x)^2 + (1 / df$SE.y)^2)) +
                df$BETA.y * (((1 / df$SE.y)^2) / ((1 / df$SE.x)^2 + (1 / df$SE.y)^2))
df$var_meta <- 1 / ((1 / df$SE.x)^2 + (1 / df$SE.y)^2)
df$P_meta_reg <- exp(pchisq((df$beta_meta^2) / df$var_meta, df = 1, lower.tail = FALSE, log.p = TRUE))
df$P_meta_Z <- 2 * exp(pnorm(-abs((sqrt(N1) * probitlink(df$pval.x/2, bvalue = .Machine$double.eps) * sign(df$BETA.x) +
                                   sqrt(N2) * probitlink(df$pval.y/2, bvalue = .Machine$double.eps) * sign(df$BETA.y)) /
                                  sqrt(N1 + N2)), log.p = TRUE))
df$Q <- (((df$BETA.x - df$beta_meta)^2) / (df$SE.x)^2) + (((df$BETA.y - df$beta_meta)^2) / (df$SE.y^2))
df$I_squared <- pmax(0, (df$Q - 1) / df$Q) * 100
write.table(df, 'GSI_Meta.txt', row.names = FALSE, quote = FALSE)