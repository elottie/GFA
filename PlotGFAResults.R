suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize)) # for color functions, if desired
suppressPackageStartupMessages(library(biomaRt))

setwd('./gfa_results/')

loading_files <- list.files(pattern = "gls_loadings.*\\d{1,2}\\.RDS$")

all_loadings <- lapply(loading_files, readRDS)
combined_loadings <- do.call(rbind, all_loadings)

# select rows of interest
combinedZsc <- combined_loadings[combined_loadings$factor1.p < 0.05 |
                                 combined_loadings$factor2.p < 0.05 |
                                 combined_loadings$factor3.p < 0.05,]
# select columns of interest
combinedZscCol <- as.matrix(combinedZsc[,c('factor1.z',
                                        'factor2.z',
                                        'factor3.z')])
rownames(combinedZscCol) <- combinedZsc$snp

# select part for plotting
quantile(combinedZscCol[,'factor1.z'])
sel_rows <- which(combinedZscCol[, "factor1.z"] > -1.15 & combinedZsc[, "factor1.z"] < 1.15)
plotCombZsc <- combinedZscCol[sel_rows[1:20], ]

setwd('/nfs/turbo/sph-jvmorr/GFA_metabolites_2025/')
png("GFA_FactorLoadings_Heatmap_20MiddleSNPs.png", width = 1200, height = 900, res = 150)
Heatmap(
  plotCombZsc,
  name = "Loadings",
  column_title = "Factor Loadings: 20 SNPs, in 2nd & 3rd Quants of Fac1 Zsc, Sig Assoc w/ 1+ Fac",
  row_title = "SNP",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_heatmap_legend = TRUE,
  col = colorRamp2(c(min(plotCombZsc, na.rm=TRUE), 0, max(plotCombZsc, na.rm=TRUE)),
                   c("blue", "white", "red"))
)
dev.off()

# okay, what do these snps do
mart <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
# rsids <- rownames(plotCombZsc)
rsids <- rownames(combinedZscCol)[sel_rows]
attrs <- c("refsnp_id", "chr_name", "chrom_start", "allele", "minor_allele", "minor_allele_freq",
           "associated_gene", "consequence_type_tv", "phenotype_name", "clinical_significance")
annot <- getBM(attributes=attrs, filters="snp_filter", values=rsids, mart=mart)
# print(annot)

# whichmin_rc <- which(combinedZsc == min(combinedZsc), arr.ind = TRUE)
# # Or to get the entire row:
# combinedZsc[whichmin_rc[1], ]
