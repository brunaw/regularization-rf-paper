#---------------------------------------------------------------------
# Just run the following code to install the used ranger version'
# please install devtools first if needed 
# install.packages("devtools")
devtools::install("code/ranger")

# Test if worked
library(ranger)
sessionInfo() #  [1] ranger_0.12.2 
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Installing other packages

# From CRAN
packages <- c("tidyverse", "furrr", "infotheo", "rminer", 
              "RColorBrewer", "extrafont", "patchwork", 
              "RRF", "data.table", "xtable")
install.packages(packages)
#---------------------------------------------------------------------

