# Heroku R buildpack requires this file
# Install any R packages needed for CIE

# Base packages needed for statistics/data processing
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("data.table")) install.packages("data.table")
