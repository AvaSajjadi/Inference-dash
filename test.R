library(QuaternaryProd)

df <- data.frame(
  entrez = as.integer(c(1,2,3,4,5)),
  fc = as.numeric(c(1.5, -1.2, 0.3, -2.0, 1.8)),
  pvalue = as.numeric(c(0.01, 0.02, 0.5, 0.001, 0.04)),
  stringsAsFactors = FALSE
)

res <- RunCRE_HSAStringDB(
  gene_expression_data = df,
  method = "Quaternary",
  fc.thresh = log(1.3),
  pval.thresh = 0.05
)

print(head(res$regulators))
