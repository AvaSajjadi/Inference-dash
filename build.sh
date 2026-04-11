#!/bin/bash
set -e

# Install R packages for CIE analysis
echo "Installing R packages..."

R --vanilla --quiet --slave <<EOF
install.packages("dplyr", repos="https://cloud.r-project.org")
install.packages("magrittr", repos="https://cloud.r-project.org")
install.packages("data.table", repos="https://cloud.r-project.org")
if (!require("remotes")) install.packages("remotes", repos="https://cloud.r-project.org")
remotes::install_github("cansylab/CIE", upgrade="never")
EOF

echo "R packages installed successfully"
