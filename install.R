# Install R packages needed for CIE analysis
# This runs during Railway/Heroku build

install.packages("dplyr", repos="https://cloud.r-project.org")
install.packages("magrittr", repos="https://cloud.r-project.org")
install.packages("data.table", repos="https://cloud.r-project.org")

# CIE package from GitHub
if (!require("remotes")) install.packages("remotes", repos="https://cloud.r-project.org")
remotes::install_github("cansylab/CIE", upgrade="never")
