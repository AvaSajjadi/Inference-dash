# Install R packages needed for CIE analysis

cat("Installing R packages...\n")

# Core packages
packages <- c("dplyr", "magrittr", "data.table")
for (pkg in packages) {
    cat("Installing", pkg, "...\n")
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        if (!install.packages(pkg, repos="https://cloud.r-project.org", quiet = FALSE)) {
            stop(paste("Failed to install", pkg))
        }
    }
}

# CIE package from GitHub
cat("Installing remotes...\n")
if (!require("remotes", character.only = TRUE, quietly = TRUE)) {
    if (!install.packages("remotes", repos="https://cloud.r-project.org", quiet = FALSE)) {
        stop("Failed to install remotes")
    }
}

cat("Installing CIE from GitHub...\n")
if (!require("CIE", character.only = TRUE, quietly = TRUE)) {
    if (!remotes::install_github("cansylab/CIE", upgrade="never")) {
        stop("Failed to install CIE")
    }
}

cat("All packages installed successfully!\n")
