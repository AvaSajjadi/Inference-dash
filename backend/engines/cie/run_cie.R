#!/usr/bin/env Rscript

base_pkgs <- c("datasets","utils","grDevices","graphics","stats","methods","dplyr","magrittr")
Sys.setenv(R_DEFAULT_PACKAGES = paste(base_pkgs, collapse=","))

suppressPackageStartupMessages({
  library(utils)
  library(stats)
  library(dplyr)
  library(magrittr)
  library(CIE)
})

cat("PROGRESS: 3\n")

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
read_any <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)

  df <- tryCatch(
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || ncol(df) <= 1) {
    df <- tryCatch(
      read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  }
  if (is.null(df)) stop("Could not read file: ", path)
  df
}

get_arg <- function(args, flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

pick_first_col <- function(df, candidates) {
  nms <- tolower(colnames(df))
  idx <- match(tolower(candidates), nms, nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx) == 0) return(NULL)
  colnames(df)[idx[1]]
}

write_auto <- function(df, path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    write.csv(df, file = path, row.names = FALSE, quote = FALSE)
  } else {
    write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

make_sidecar_path <- function(base_path, suffix_with_ext) {
  paste0(base_path, suffix_with_ext)
}

normalize_method <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return("Fisher")
  x
}

normalize_database_type <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return("ChIP")
  xl <- tolower(x)
  if (xl %in% c("tcchip", "chip")) return("ChIP")
  x
}

ensure_signature <- function(sig) {
  colnames(sig) <- tolower(colnames(sig))

  entrez_col <- pick_first_col(sig, c("entrez", "entrez_id", "entrezid", "geneid", "gene_id"))
  fc_col     <- pick_first_col(sig, c("fc", "logfc", "log2fc", "fold_change", "foldchange"))
  pval_col   <- pick_first_col(sig, c("pval", "pvalue", "p_value", "padj", "fdr", "qvalue"))

  if (is.null(entrez_col)) stop("signature must contain an entrez column")
  if (is.null(fc_col)) stop("signature must contain a fold-change column")

  if (is.null(pval_col)) {
    sig$pval <- 1
    pval_col <- "pval"
  }

  out <- data.frame(
    entrez = suppressWarnings(as.integer(sig[[entrez_col]])),
    fc     = suppressWarnings(as.numeric(sig[[fc_col]])),
    pval   = suppressWarnings(as.numeric(sig[[pval_col]])),
    stringsAsFactors = FALSE
  )

  out <- out[!is.na(out$entrez) & !is.na(out$fc), , drop = FALSE]
  out$pval[is.na(out$pval)] <- 1
  out
}

ensure_rels <- function(rels) {
  uid_col  <- pick_first_col(rels, c("uid"))
  src_col  <- pick_first_col(rels, c("srcuid", "source", "src"))
  trg_col  <- pick_first_col(rels, c("trguid", "target", "trg"))
  type_col <- pick_first_col(rels, c("type", "mode", "relation"))

  if (is.null(uid_col))  stop("rels must contain uid")
  if (is.null(src_col))  stop("rels must contain srcuid")
  if (is.null(trg_col))  stop("rels must contain trguid")
  if (is.null(type_col)) stop("rels must contain type")

  out <- data.frame(
    uid    = suppressWarnings(as.integer(rels[[uid_col]])),
    srcuid = suppressWarnings(as.integer(rels[[src_col]])),
    trguid = suppressWarnings(as.integer(rels[[trg_col]])),
    type   = as.character(rels[[type_col]]),
    stringsAsFactors = FALSE
  )

  out <- out[!is.na(out$uid) & !is.na(out$srcuid) & !is.na(out$trguid), , drop = FALSE]
  out
}

ensure_ents <- function(ents, rels_clean) {
  uid_col  <- pick_first_col(ents, c("uid"))
  name_col <- pick_first_col(ents, c("name", "symbol", "gene_symbol"))
  id_col   <- pick_first_col(ents, c("id", "entrez", "geneid", "gene_id"))
  type_col <- pick_first_col(ents, c("type"))

  if (!is.null(uid_col)) {
    out <- data.frame(
      uid  = suppressWarnings(as.integer(ents[[uid_col]])),
      name = if (!is.null(name_col)) as.character(ents[[name_col]]) else as.character(ents[[uid_col]]),
      id   = if (!is.null(id_col)) as.character(ents[[id_col]]) else as.character(ents[[uid_col]]),
      type = if (!is.null(type_col)) as.character(ents[[type_col]]) else "unknown",
      stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$uid), , drop = FALSE]
    out <- out[!duplicated(out$uid), , drop = FALSE]
    return(out)
  }

  all_uids <- sort(unique(c(rels_clean$uid, rels_clean$srcuid, rels_clean$trguid)))
  srcs <- unique(rels_clean$srcuid)

  out <- data.frame(
    uid  = all_uids,
    name = as.character(all_uids),
    id   = as.character(all_uids),
    type = ifelse(all_uids %in% srcs, "Protein", "mRNA"),
    stringsAsFactors = FALSE
  )
  out
}

extract_first_df <- function(x) {
  if (is.data.frame(x)) return(x)

  if (is.list(x)) {
    for (nm in names(x)) {
      if (is.data.frame(x[[nm]])) return(x[[nm]])
    }
    for (nm in names(x)) {
      if (is.list(x[[nm]]) && length(x[[nm]]) > 0 && is.data.frame(x[[nm]][[1]])) {
        return(x[[nm]][[1]])
      }
    }
  }

  data.frame()
}

extract_regulators <- function(res) {
  if (is.data.frame(res)) return(res)

  if (is.list(res)) {
    preferred_names <- c(
      "regulators", "regulator", "tfs", "tf", "results", "result",
      "stats", "summary"
    )

    for (nm in preferred_names) {
      if (!is.null(res[[nm]]) && is.data.frame(res[[nm]])) return(res[[nm]])
    }

    for (nm in names(res)) {
      if (grepl("reg", nm, ignore.case = TRUE) && is.data.frame(res[[nm]])) return(res[[nm]])
    }

    return(extract_first_df(res))
  }

  data.frame()
}

extract_pathways <- function(res) {
  if (!is.list(res)) return(data.frame())

  for (nm in names(res)) {
    if (grepl("path", nm, ignore.case = TRUE) && is.data.frame(res[[nm]])) {
      return(res[[nm]])
    }
  }

  for (nm in names(res)) {
    obj <- res[[nm]]
    if (is.list(obj)) {
      for (nm2 in names(obj)) {
        if (grepl("path", nm2, ignore.case = TRUE) && is.data.frame(obj[[nm2]])) {
          return(obj[[nm2]])
        }
      }
    }
  }

  data.frame()
}

make_edges_from_signature <- function(sig_df, rels_clean, ents_clean, tf_df, p_thresh = 0.05, fc_thresh = log2(1.5)) {
  sig_keep <- sig_df[which(sig_df$pval <= p_thresh & abs(sig_df$fc) >= fc_thresh), , drop = FALSE]
  if (nrow(sig_keep) == 0) {
    return(data.frame(srcuid = integer(), trguid = integer(), score = integer()))
  }

  sig_keep$sign <- ifelse(sig_keep$fc > 0, 1L, ifelse(sig_keep$fc < 0, -1L, 0L))

  ents_map <- ents_clean
  ents_map$id_num <- suppressWarnings(as.integer(ents_map$id))

  mapped <- merge(
    sig_keep,
    ents_map[, c("uid", "id_num"), drop = FALSE],
    by.x = "entrez",
    by.y = "id_num",
    all = FALSE
  )

  if (nrow(mapped) == 0) {
    return(data.frame(srcuid = integer(), trguid = integer(), score = integer()))
  }

  trg_sign <- mapped$sign
  names(trg_sign) <- mapped$uid

  rels2 <- rels_clean[rels_clean$trguid %in% names(trg_sign), , drop = FALSE]
  if (nrow(rels2) == 0) {
    return(data.frame(srcuid = integer(), trguid = integer(), score = integer()))
  }

  rel_sign <- ifelse(
    tolower(rels2$type) %in% c("increase", "increases", "activation", "activates", "up"),
    1L,
    ifelse(
      tolower(rels2$type) %in% c("decrease", "decreases", "repression", "represses", "down"),
      -1L,
      0L
    )
  )

  obs_sign <- trg_sign[as.character(rels2$trguid)]
  score <- ifelse(rel_sign == obs_sign, 1L, ifelse(rel_sign == -obs_sign, -1L, 0L))

  edge_df <- data.frame(
    srcuid = rels2$srcuid,
    trguid = rels2$trguid,
    score  = score,
    stringsAsFactors = FALSE
  )

  edge_df <- unique(edge_df)

  if (nrow(tf_df) > 0) {
    uid_candidates <- c("uid", "id")
    uid_col <- uid_candidates[uid_candidates %in% tolower(colnames(tf_df))][1]
    if (!is.na(uid_col)) {
      real_col <- colnames(tf_df)[match(uid_col, tolower(colnames(tf_df)))]
      keep_src <- suppressWarnings(as.integer(tf_df[[real_col]]))
      keep_src <- keep_src[!is.na(keep_src)]
      if (length(keep_src) > 0) {
        edge_df <- edge_df[edge_df$srcuid %in% keep_src, , drop = FALSE]
      }
    }
  }

  edge_df
}

# ------------------------------------------------------------
# Args
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

sig_file  <- get_arg(args, "-s", get_arg(args, "--signature"))
out_path  <- get_arg(args, "-o", "cie_edges.csv")
rels_file <- get_arg(args, "--rels")
ents_file <- get_arg(args, "--ents")

db_arg    <- get_arg(args, "--db", "tcChIP")
tissue    <- get_arg(args, "--tissue", "all")
method    <- normalize_method(get_arg(args, "-m", "Fisher"))

p_thresh <- suppressWarnings(as.numeric(get_arg(args, "-p", "0.05")))
if (is.na(p_thresh)) p_thresh <- 0.05

fc_thresh_user <- suppressWarnings(as.numeric(get_arg(args, "-f", "1.5")))
if (is.na(fc_thresh_user)) fc_thresh_user <- 1.5

if (is.null(sig_file) || is.null(rels_file) || is.null(ents_file)) {
  stop("Missing required args: -s <signature> -o <out> --rels <rels> --ents <ents>")
}

out_dir <- dirname(out_path)
if (identical(out_dir, ".") || is.na(out_dir) || out_dir == "") out_dir <- "."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Load and normalize
# ------------------------------------------------------------
sig_raw  <- read_any(sig_file)
rels_raw <- read_any(rels_file)
ents_raw <- read_any(ents_file)

cat("PROGRESS: 10\n")

sig_df  <- ensure_signature(sig_raw)
rels_df <- ensure_rels(rels_raw)
ents_df <- ensure_ents(ents_raw, rels_df)

cat("PROGRESS: 18\n")

# Professor-style threshold for logFC input:
# user passes fold-change threshold 1.5 -> convert to log2(1.5)
fc_thresh_log2 <- log2(fc_thresh_user)

# ------------------------------------------------------------
# Real CIE call
# ------------------------------------------------------------
res <- runCIE(
  databaseType = normalize_database_type(db_arg),
  filter = TRUE,
  DGEs = sig_df,
  p.thresh = p_thresh,
  fc.thresh = fc_thresh_log2,
  logFC = TRUE,
  methods = method,
  ents = ents_df,
  rels = rels_df,
  useFile = FALSE,
  verbose = FALSE,
  numCores = 1
)

cat("PROGRESS: 25\n")

# ------------------------------------------------------------
# Extract outputs
# ------------------------------------------------------------
tf_df <- extract_regulators(res)
if (is.null(tf_df) || !is.data.frame(tf_df)) tf_df <- data.frame()

pathway_df <- extract_pathways(res)
if (is.null(pathway_df) || !is.data.frame(pathway_df)) pathway_df <- data.frame()

edge_df <- make_edges_from_signature(
  sig_df = sig_df,
  rels_clean = rels_df,
  ents_clean = ents_df,
  tf_df = tf_df,
  p_thresh = p_thresh,
  fc_thresh = fc_thresh_log2
)

# ------------------------------------------------------------
# Write outputs
# Contract used by app1.py:
#   -o <path> is the EDGE file
#   TF table is written as a sidecar:
#       <out>_tfs.tsv
#   Pathway table is written as a sidecar if available:
#       <out>_pathwayEnrichment.tsv
# ------------------------------------------------------------
write_auto(edge_df, out_path)

tf_path <- make_sidecar_path(out_path, "_tfs.tsv")
write.table(
  tf_df,
  file = tf_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

if (nrow(pathway_df) > 0) {
  pathway_path <- make_sidecar_path(out_path, "_pathwayEnrichment.tsv")
  write.table(
    pathway_df,
    file = pathway_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

meta_path <- file.path(out_dir, "cie_run_meta.json")
meta <- list(
  signature = sig_file,
  rels = rels_file,
  ents = ents_file,
  out_edges = out_path,
  out_tfs = tf_path,
  database_arg = db_arg,
  tissue = tissue,
  databaseType_used = normalize_database_type(db_arg),
  method = method,
  p_thresh = p_thresh,
  fc_thresh_user = fc_thresh_user,
  fc_thresh_log2 = fc_thresh_log2,
  logFC = TRUE,
  filter = TRUE,
  rows_signature = nrow(sig_df),
  rows_rels = nrow(rels_df),
  rows_ents = nrow(ents_df),
  rows_edges = nrow(edge_df),
  rows_tfs = nrow(tf_df),
  rows_pathways = nrow(pathway_df)
)
writeLines(jsonlite::toJSON(meta, pretty = TRUE, auto_unbox = TRUE), meta_path)

cat("PROGRESS: 100\n")
cat("CIE run completed\n")
