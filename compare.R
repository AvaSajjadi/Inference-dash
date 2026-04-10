args <- commandArgs(trailingOnly=TRUE)

student_file <- args[1]
prof_file <- args[2]

student <- read.delim(student_file, stringsAsFactors=FALSE)
prof <- read.delim(prof_file, stringsAsFactors=FALSE)

# Compare using Entrez IDs
student_ids <- unique(student$id)
prof_ids <- unique(prof$id)

common <- intersect(student_ids, prof_ids)

cat("\n========== CIE RESULTS COMPARISON ==========\n\n")

cat("Student TF count:", length(student_ids), "\n")
cat("Professor TF count:", length(prof_ids), "\n")
cat("Common TFs:", length(common), "\n\n")

jaccard <- length(common) / length(union(student_ids, prof_ids))

cat("Jaccard similarity:", round(jaccard,3), "\n\n")

cat("Top overlapping TF IDs:\n")
print(head(common,20))
