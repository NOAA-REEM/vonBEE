# ----------------------------------------
# setup.R
# subset of Yasumiishi et al. in prep
# EBS salmon code
# 
# updated 2022
# ----------------------------------------

lib_list <- c(
  # these for reshaping and manipulating data:
    "ncdf4",
    "devtools",
    "TMB",
   # "ncfd",
    "magrittr",
    "httr",
    "reshape",
    "dplyr", 
    
  # markdown stuff:
    "knitr",
    "kableExtra",
    
  # These for making plots:
    "RColorBrewer",
    "ggplot2", 
    "mgcv",
    "cowplot",               # 
    "wesanderson",
    "scales",
    "ggforce",
    "grid",
    "processx",
    "plotly",
    "extrafont"
  )

# Install missing libraries:
missing <- setdiff(lib_list, installed.packages()[, 1])
if (length(missing) > 0) install.packages(missing)

# Load libraries:
for(lib in lib_list)
       eval(parse(text=paste("library(",lib,")")))



