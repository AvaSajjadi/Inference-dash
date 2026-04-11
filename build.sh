#!/bin/bash
set -e

echo "Setting up R environment..."

# Find Rscript location
RSCRIPT_PATH=$(which Rscript)
if [ -z "$RSCRIPT_PATH" ]; then
    echo "ERROR: Rscript not found in PATH"
    exit 1
fi

echo "Found Rscript at: $RSCRIPT_PATH"

# Save Rscript path for runtime access
echo "$RSCRIPT_PATH" > /app/.rscript_path
chmod 644 /app/.rscript_path

# Install R packages for CIE analysis
echo "Installing R packages..."

"$RSCRIPT_PATH" --vanilla --quiet --slave <<EOF
install.packages("dplyr", repos="https://cloud.r-project.org")
install.packages("magrittr", repos="https://cloud.r-project.org")
install.packages("data.table", repos="https://cloud.r-project.org")
if (!require("remotes")) install.packages("remotes", repos="https://cloud.r-project.org")
remotes::install_github("cansylab/CIE", upgrade="never")
EOF

echo "R packages installed successfully"
echo "Rscript path saved to: /app/.rscript_path"
