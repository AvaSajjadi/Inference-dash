args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript run_cie.R expr.csv network.csv out.csv")
}

expr_file <- args[1]
net_file  <- args[2]
out_file  <- args[3]

suppressPackageStartupMessages({
  library(QuaternaryProd)
})

# TODO: replace this with the real pipeline you and your professor use.

expr <- read.csv(expr_file, stringsAsFactors = FALSE)
net  <- read.csv(net_file, stringsAsFactors = FALSE)

# Temporary fake result so the Dash app can run end-to-end.
# You will later replace this with the actual CIE/QuaternaryProd analysis.
res <- data.frame(
  TF      = head(expr[[1]], 10),
  score   = seq_len(min(10, nrow(expr))),
  p_value = seq(0.001, 0.01, length.out = min(10, nrow(expr)))
)

write.csv(res, out_file, row.names = FALSE)
