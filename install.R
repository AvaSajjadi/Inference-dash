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

# CIE package installation
cat("Installing remotes and devtools...\n")
if (!require("remotes", character.only = TRUE, quietly = TRUE)) {
    install.packages("remotes", repos="https://cloud.r-project.org", quiet = FALSE)
}
if (!require("devtools", character.only = TRUE, quietly = TRUE)) {
    install.packages("devtools", repos="https://cloud.r-project.org", quiet = FALSE)
}

cat("Installing CIE...\n")
if (!require("CIE", character.only = TRUE, quietly = TRUE)) {
    # Set GitHub PAT if available (from Docker build arg)
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
        cat("  docker build --build-arg GITHUB_TOKEN=<your-token> -t inference-dash .\n")
        cat("Or install locally: remotes::install_github('cansylab/CIE')\n")
        cat("Continuing without CIE. ORNOR analysis will work.\n")
    })
}

cat("All available packages installed successfully!\n")
