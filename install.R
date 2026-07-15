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

# CIE package installation — use remotes only (devtools has a heavy dep tree
# that needs fontconfig/freetype dev headers not available on this base image)
cat("Checking remotes...\n")
if (!require("remotes", character.only = TRUE, quietly = TRUE)) {
    install.packages("remotes", repos="https://cloud.r-project.org", quiet = FALSE)
}

cat("Installing CIE...\n")
if (!require("CIE", character.only = TRUE, quietly = TRUE)) {
    github_pat <- Sys.getenv("GITHUB_PAT", "")
    if (nzchar(github_pat)) {
        cat("Using GitHub Personal Access Token for authentication\n")
        Sys.setenv(GITHUB_PAT = github_pat)
    }

    tryCatch({
        cat("Attempting installation from GitHub...\n")
        remotes::install_github("cansylab/CIE", upgrade="never", dependencies=TRUE)
        cat("CIE installed successfully!\n")
    }, error = function(e) {
        cat("Warning: CIE installation failed:", e$message, "\n")
        cat("To install CIE, provide a GitHub token at build time:\n")
        cat("  fly deploy --build-arg GITHUB_TOKEN=<your-token>\n")
        cat("Continuing without CIE. ORNOR analysis will work.\n")
    })
}

cat("All available packages installed successfully!\n")
