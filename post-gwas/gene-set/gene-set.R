# Gene-set analysis with GAUSS (Gene-set Analysis Using Sparse Signals).

library(GAUSS)
library(data.table)
library(dplyr)
library(foreach)

gmt <- '/gauss/c2.cp.kegg.v2023.1.Hs.symbols.comma_parsed.txt'  #change to gene set file path
phenotypes <- c("GSI", "SRM", "SRSR")

# Process each phenotype
for (pheno in phenotypes) {
    tryCatch({	
        GAUSS_All(
            summary_file = paste('/summ/', pheno, '.txt', sep = ''),
            gene_name = 1,
            pv_name = 2,
            output_file = paste('/summ/c2/', pheno, '_gauss.txt', sep = ''),
            gmt = gmt,
            ags = "def",
            verbose = TRUE,
            parallel = FALSE,
            is.appx = TRUE
        )
        closeAllConnections()
    }, error = function(e) {
        cat("Error occurred for phenotype:", pheno, "\n", conditionMessage(e), "\n")
    })
}
