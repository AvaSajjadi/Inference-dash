#!/usr/bin/env Rscript

# ------------------------------------------------------------
# run_cie.R — Prof-match CIE runner (tcChIP default, Fisher default)
# + Reproducibility upgrades:
#   - writes <out>.meta.json with parameters, hashes, network/signature stats, package versions
#   - robust signature col normalization to EXACT CIE-required: entrez, fc, pvalue
#   - robust DGEs argument name detection (DGEs vs DEGs)
# ------------------------------------------------------------

# Ensure dplyr::n() exists even inside parallel workers
base_pkgs <- c("datasets","utils","grDevices","graphics","stats","methods","dplyr","magrittr")
Sys.setenv(R_DEFAULT_PACKAGES = paste(base_pkgs, collapse=","))

suppressPackageStartupMessages({
  library(optparse)
  library(tools)
  library(dplyr)
  library(magrittr)
  library(rjson)  # for metadata JSON
})

# -----------------------------
# CLI options
# -----------------------------
option_list <- list(
  make_option(c("-s","--sig"), type="character", help="Signature file (tsv/txt/csv)"),
  make_option(c("-o","--out"), type="character", default=NULL, help="Output path (tsv recommended)"),

  # Database selection (prof defaults)
  make_option(c("--db"), type="character", default="tcChIP",
              help="Database type: tcChIP (default), ChIP, TRRUST, STRINGdb, TRED"),
  make_option(c("--tissue"), type="character", default="all",
              help="tcChIP tissue (default: all). If not all, will look for <tissue>.rels in tcChIP dir."),

  # If using tcChIP locally, point to the folder containing all_tissues.rels and ChIPfilter.ents
  make_option(c("--tcchip_dir"), type="character", default=NULL,
              help="Path to tissueCorrectedChIP folder (contains all_tissues.rels and ChIPfilter.ents)."),

  # Enrichment parameters (prof defaults)
  make_option(c("-m","--method"), type="character", default="Fisher",
              help="Method: Fisher (default), Enrichment, Ternary, Quaternary"),
  make_option(c("-p","--pval"), type="double", default=0.05, help="p-value threshold (default 0.05)"),
  make_option(c("-f","--fc"), type="double", default=1.5, help="fold-change threshold (default 1.5)"),
  make_option(c("-u","--use_log_fc"), type="integer", default=1, help="Use log fold change? 1/0 (default 1)"),
  make_option(c("-b","--log_base"), type="double", default=NA, help="log base (optional; usually NA)"),
  make_option(c("-c","--cores"), type="integer", default=1, help="cores (default 1)"),

  make_option(c("--debug"), type="integer", default=0, help="debug 1/0")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$sig) || !nzchar(opt$sig)) {
  cat("ERROR: -s/--sig is required\n", file=stderr())
  quit(status=2)
}

# Default output name (prof-like) if not provided
if (is.null(opt$out) || !nzchar(opt$out)) {
  sig_base <- basename(opt$sig)
  opt$out <- file.path("results", paste0(opt$method, opt$db, sig_base, ".txt.tsv"))
}

# Normalize db name (case-insensitive)
db_in <- opt$db
db_norm <- tolower(db_in)
db_type <- switch(db_norm,
  "tcchip" = "tcChIP",
  "chip"   = "ChIP",
  "trrust" = "TRRUST",
  "stringdb" = "STRINGdb",
  "tred"   = "TRED",
  db_in
)

tissue <- opt$tissue
if (is.null(tissue) || !nzchar(tissue)) tissue <- "all"

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

say <- function(...) cat(..., "\n")
die <- function(...) { cat("ERROR: ", ..., "\n", file=stderr()); quit(status=1) }

md5_file <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NA)
  # tools::md5sum returns named vector
  unname(tools::md5sum(path))
}

safe_int <- function(x) {
  y <- suppressWarnings(as.integer(as.character(x)))
  ifelse(is.na(y), as.character(x), y)
}

infer_delim <- function(first_line) {
  if (grepl("\t", first_line)) "\t" else if (grepl(",", first_line)) "," else "\t"
}

read_signature <- function(path) {
  # Reads signature and normalizes output to EXACT CIE-required cols:
  #   entrez, fc, pvalue

  first_line <- readLines(path, n=1, warn=FALSE)
  sep_guess <- infer_delim(first_line)

  df <- tryCatch(
    read.table(path, header=TRUE, sep=sep_guess, quote="", comment.char="",
               stringsAsFactors=FALSE, check.names=FALSE),
    error=function(e) {
      read.table(path, header=TRUE, sep="\t", quote="", comment.char="",
                 stringsAsFactors=FALSE, check.names=FALSE)
    }
  )

  if (ncol(df) < 2) die("Signature file has too few columns: ", path)

  # --- Entrez / gene column ---
  gene_col <- NULL
  cand_gene <- c("entrez","ENTREZ","Entrez","entrezid","EntrezID",
                 "gene","Gene","GENE","symbol","SYMBOL","id","ID","Name","name")
  for (c in cand_gene) if (c %in% names(df)) { gene_col <- c; break }
  if (is.null(gene_col)) gene_col <- names(df)[1]

  # --- p-value column ---
  p_col <- NULL
  cand_p <- c("pvalue","p.value","p_value","P.Value","PValue","p","P","pval","PVAL","padj","adj.pval","adj_pval")
  for (c in cand_p) if (c %in% names(df)) { p_col <- c; break }
  if (is.null(p_col)) p_col <- names(df)[2]

  # --- fc column ---
  # CIE expects column named "fc". If logFC=TRUE in runCIE, "fc" should contain log fold-change.
  fc_col <- NULL
  cand_fc <- c("fc","FC","foldchange","FoldChange","logFC","log2FC","log2FoldChange","log_fold_change","lfc","LFC")
  for (c in cand_fc) if (c %in% names(df)) { fc_col <- c; break }
  if (is.null(fc_col)) fc_col <- if (ncol(df) >= 3) names(df)[3] else NA_character_

  out <- data.frame(
    entrez = safe_int(df[[gene_col]]),
    fc = if (!is.na(fc_col)) suppressWarnings(as.numeric(df[[fc_col]])) else NA_real_,
    pvalue = suppressWarnings(as.numeric(df[[p_col]])),
    stringsAsFactors=FALSE
  )

  out <- out[!is.na(out$entrez) & nzchar(as.character(out$entrez)), , drop=FALSE]
  out <- out[!is.na(out$pvalue), , drop=FALSE]
  out <- out[!duplicated(out$entrez), , drop=FALSE]

  say("Signature rows after cleaning/dedup: ", nrow(out))
  return(out)
}

# Create a temp CIE database directory that matches what CIE expects
make_dbdir <- function(rels_path, ents_path, debug=0) {
  root <- file.path(tempdir(), paste0("cie_db_", paste(sample(c(letters,0:9), 12, TRUE), collapse="")))
  engine <- file.path(root, "engine")
  data_dir <- file.path(engine, "data")
  dir.create(data_dir, recursive=TRUE, showWarnings=FALSE)

  # CIE sometimes looks in engine/, sometimes in engine/data/
  file.copy(rels_path, file.path(engine, "ChIP.rels"), overwrite=TRUE)
  file.copy(rels_path, file.path(data_dir, "ChIP.rels"), overwrite=TRUE)

  file.copy(ents_path, file.path(engine, "ChIP.ents"), overwrite=TRUE)
  file.copy(ents_path, file.path(data_dir, "ChIP.ents"), overwrite=TRUE)

  if (debug == 1) {
    say("[shim] databaseDir: ", paste0(engine, "/"))
    say("[shim] rels file at: ", file.path(data_dir, "ChIP.rels"))
    say("[shim] ents file at: ", file.path(data_dir, "ChIP.ents"))
  }

  return(paste0(engine, "/"))
}

# Filter ents to only nodes present in rels (universe alignment / reproducibility)
filter_ents_to_rels <- function(rels_path, ents_path, debug=0) {
  rels <- read.table(rels_path, header=FALSE, sep="\t", quote="", comment.char="",
                     stringsAsFactors=FALSE)
  if (ncol(rels) < 2) die("rels file has <2 columns: ", rels_path)
  node_ids <- unique(c(rels[[1]], rels[[2]]))

  ents <- read.table(ents_path, header=FALSE, sep="\t", quote="", comment.char="",
                     stringsAsFactors=FALSE)
  if (ncol(ents) < 1) die("ents file has <1 column: ", ents_path)

  before <- nrow(ents)
  ents2 <- ents[ents[[1]] %in% node_ids, , drop=FALSE]
  after <- nrow(ents2)

  if (debug == 1) say("[shim] ents filtered: ", before, " -> ", after)

  out_path <- file.path(tempdir(), paste0("ents_filtered_", basename(ents_path)))
  write.table(ents2, out_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  return(out_path)
}

# Network stats for metadata
network_stats <- function(rels_path, ents_path) {
  rels <- read.table(rels_path, header=FALSE, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)
  edges <- nrow(rels)
  nodes <- unique(c(rels[[1]], rels[[2]]))
  node_count <- length(nodes)

  ents <- read.table(ents_path, header=FALSE, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)
  ents_n <- nrow(ents)

  list(
    rels_edges = edges,
    rels_nodes = node_count,
    ents_rows = ents_n
  )
}

# Signature stats for metadata
signature_stats <- function(sig_df, pthresh, fcthresh, logFC_bool) {
  # NOTE:
  # - sig_df$fc is the value passed to CIE as "fc"
  # - If logFC_bool==TRUE, many pipelines compare abs(fc) >= log(fcthresh)
  #   but some UIs pass fcthresh already in log space.
  # We record BOTH, and count using a conservative interpretation:
  #   if logFC_bool==TRUE, threshold = log(fcthresh)
  #   else threshold = fcthresh
  thr_used <- if (isTRUE(logFC_bool)) log(fcthresh) else fcthresh

  ok_p <- sig_df$pvalue <= pthresh
  ok_fc <- if (all(is.na(sig_df$fc))) rep(FALSE, nrow(sig_df)) else abs(sig_df$fc) >= thr_used

  sig_n <- sum(ok_p & ok_fc, na.rm=TRUE)
  list(
    rows = nrow(sig_df),
    pthresh = pthresh,
    fcthresh_input = fcthresh,
    fcthresh_used_for_count = thr_used,
    logFC = isTRUE(logFC_bool),
    significant_count_est = sig_n
  )
}

write_out <- function(df, out_path) {
  dir.create(dirname(out_path), recursive=TRUE, showWarnings=FALSE)
  ext <- tolower(file_ext(out_path))
  if (ext %in% c("tsv")) {
    write.table(df, out_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  } else {
    write.csv(df, out_path, row.names=FALSE, quote=TRUE)
  }
  say("Wrote: ", out_path)
}

write_meta <- function(meta, out_path) {
  meta_path <- paste0(out_path, ".meta.json")
  dir.create(dirname(meta_path), recursive=TRUE, showWarnings=FALSE)

  # rjson wants a list, not a data.frame
  json_txt <- rjson::toJSON(meta)
  writeLines(json_txt, con=meta_path, useBytes=TRUE)
  say("Wrote: ", meta_path)
}

# Robust call into CIE, adapting to whatever argument names exist in this install
run_cie_safe <- function(sig_df, databaseType, db_dir, pthresh, fcthresh, logFC_bool, methods, cores, debug=0) {
  suppressPackageStartupMessages(library(CIE))

  fmls <- names(formals(CIE::runCIE))
  if (debug == 1) {
    say("CIE::runCIE formals:")
    print(fmls)
  }

  args <- list()

  if ("databaseType" %in% fmls) args$databaseType <- databaseType
  if ("filter" %in% fmls) args$filter <- FALSE

  # signature argument name differs by install
  if ("DGEs" %in% fmls) {
    args$DGEs <- sig_df
  } else if ("DEGs" %in% fmls) {
    args$DEGs <- sig_df
  } else if ("DEG" %in% fmls) {
    args$DEG <- sig_df
  } else {
    args[[length(args)+1]] <- sig_df
  }

  if ("p.thresh" %in% fmls) args$p.thresh <- pthresh
  if ("fc.thresh" %in% fmls) args$fc.thresh <- fcthresh
  if ("logFC" %in% fmls) args$logFC <- logFC_bool

  if ("methods" %in% fmls) args$methods <- methods
  if ("databaseDir" %in% fmls) args$databaseDir <- db_dir

  # IMPORTANT for your install: defaults show useFile=TRUE, so we explicitly set it when available
  if ("useFile" %in% fmls) args$useFile <- TRUE

  if ("expectProgressObject" %in% fmls) args$expectProgressObject <- FALSE
  if ("verbose" %in% fmls) args$verbose <- (debug == 1)

  if ("numCores" %in% fmls) args$numCores <- as.integer(cores)
  if ("progress" %in% fmls) args$progress <- NULL

  res <- do.call(CIE::runCIE, args)
  return(res)
}

# ------------------------------------------------------------
# Resolve tcChIP rels/ents paths (prof default)
# ------------------------------------------------------------

sig_df <- read_signature(opt$sig)

rels_path <- NULL
ents_path <- NULL
tcdir_used <- NULL
tissue_file_used <- NULL

if (tolower(db_type) == "tcchip") {
  tcdir <- opt$tcchip_dir
  if (is.null(tcdir) || !nzchar(tcdir)) {
    guess <- path.expand("~/Downloads/Lab/CIEdata/root/CIEdata/tissueCorrectedChIP")
    if (dir.exists(guess)) tcdir <- guess
  }
  if (is.null(tcdir) || !dir.exists(tcdir)) {
    die("tcChIP selected but tcchip_dir not found. Provide --tcchip_dir PATH.")
  }
  tcdir_used <- normalizePath(tcdir, winslash="/", mustWork=TRUE)

  tissue_file <- if (tolower(tissue) == "all") "all_tissues.rels" else paste0(tissue, ".rels")
  tissue_file_used <- tissue_file
  rels_path <- file.path(tcdir_used, tissue_file)
  if (!file.exists(rels_path)) {
    die("tcChIP rels not found: ", rels_path, " (tissue=", tissue, ")")
  }

  ents_candidate <- file.path(tcdir_used, "ChIPfilter.ents")
  if (!file.exists(ents_candidate)) {
    die("tcChIP ents not found. Expected ChIPfilter.ents at: ", ents_candidate)
  }

  ents_path <- filter_ents_to_rels(rels_path, ents_candidate, debug=opt$debug)
} else {
  die("For now this runner is prof-match focused for tcChIP. (db=", db_type, ")")
}

# ------------------------------------------------------------
# Build CIE databaseDir shim + run
# ------------------------------------------------------------

db_dir <- make_dbdir(rels_path, ents_path, debug=opt$debug)

say(sprintf("Requested (prof-match): db=%s tissue=%s method=%s p.thresh=%.3g fc.thresh=%.3g logFC=%s cores=%d",
            db_type, tissue, opt$method, opt$pval, opt$fc, as.character(opt$use_log_fc==1), opt$cores))

pb <- txtProgressBar(min=0, max=100, style=3)
setTxtProgressBar(pb, 5)

res <- tryCatch({
  setTxtProgressBar(pb, 15)
  run_cie_safe(
    sig_df = sig_df,
    databaseType = "ChIP",     # tcChIP is implemented as ChIP network files on disk
    db_dir = db_dir,
    pthresh = opt$pval,
    fcthresh = opt$fc,
    logFC_bool = (opt$use_log_fc == 1),
    methods = opt$method,
    cores = opt$cores,
    debug = opt$debug
  )
}, error=function(e) {
  close(pb)
  cat("\nCIE ERROR: ", conditionMessage(e), "\n", file=stderr())
  quit(status=1)
})

setTxtProgressBar(pb, 85)

# ------------------------------------------------------------
# Normalize output to match professor file conventions
# ------------------------------------------------------------
tf_tbl <- NULL

if (is.data.frame(res)) {
  tf_tbl <- res
} else if (is.list(res)) {
  if (!is.null(res$regulators) && is.data.frame(res$regulators)) tf_tbl <- res$regulators
  if (is.null(tf_tbl) && !is.null(res$results) && is.data.frame(res$results)) tf_tbl <- res$results
  if (is.null(tf_tbl)) {
    for (x in res) {
      if (is.data.frame(x)) { tf_tbl <- x; break }
    }
  }
}

if (is.null(tf_tbl) || !is.data.frame(tf_tbl)) {
  close(pb)
  die("Could not locate TF table in CIE result object.")
}

if (!("isTF" %in% names(tf_tbl))) tf_tbl$isTF <- TRUE

front <- c("uid","name","id","type","isTF")
front <- front[front %in% names(tf_tbl)]
rest <- setdiff(names(tf_tbl), front)
tf_tbl <- tf_tbl[, c(front, rest), drop=FALSE]

setTxtProgressBar(pb, 100)
close(pb)
say("\"Complete!\"")

write_out(tf_tbl, opt$out)

# ------------------------------------------------------------
# Write reproducibility metadata JSON
# ------------------------------------------------------------
pkg_ver <- function(pkg) {
  v <- tryCatch(as.character(utils::packageVersion(pkg)), error=function(e) NA)
  v
}

net_stats <- network_stats(rels_path, ents_path)
sig_stats <- signature_stats(sig_df, opt$pval, opt$fc, (opt$use_log_fc == 1))

meta <- list(
  engine = "CIE",
  method = opt$method,
  parameters = list(
    db = db_type,
    tissue = tissue,
    p_thresh = opt$pval,
    fc_thresh_input = opt$fc,
    logFC = (opt$use_log_fc == 1),
    cores = opt$cores
  ),
  signature = list(
    path = normalizePath(opt$sig, winslash="/", mustWork=TRUE),
    md5 = md5_file(opt$sig),
    stats = sig_stats
  ),
  network = list(
    tcchip_dir = tcdir_used,
    tissue_rels_file = tissue_file_used,
    rels_path = normalizePath(rels_path, winslash="/", mustWork=TRUE),
    ents_path = normalizePath(ents_path, winslash="/", mustWork=TRUE),
    rels_md5 = md5_file(rels_path),
    ents_md5 = md5_file(ents_path),
    stats = net_stats
  ),
  shim = list(
    databaseDir = db_dir
  ),
  packages = list(
    CIE = pkg_ver("CIE"),
    QuaternaryProd = pkg_ver("QuaternaryProd"),
    fdrtool = pkg_ver("fdrtool"),
    dplyr = pkg_ver("dplyr")
  ),
  created_at = as.character(Sys.time())
)

write_meta(meta, opt$out)

quit(status=0)
