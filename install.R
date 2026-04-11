# Install R packages needed for CIE analysis

cat("Installing R packages...\n")

# Core packages
packages <- c("dplyr", "magrittr", "data.table")
for (pkg in packages) {
    cat("Installing", pkg, "...\n")
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, repos="https://cloud.r-project.org", quiet = FALSE)
    }
}

# CIE package from GitHub
cat("Installing remotes...\n")
if (!require("remotes", character.only = TRUE, quietly = TRUE)) {
    install.packages("remotes", repos="https://cloud.r-project.org", quiet = FALSE)
}

cat("Installing CIE from GitHub...\n")
if (!require("CIE", character.only = TRUE, quietly = TRUE)) {
    remotes::install_github("cansylab/CIE", upgrade="never")
}

cat("All packages installed successfully!\n")
