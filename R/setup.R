# ----------------------------------------
# setup.R
# subset of Yasumiishi et al. in prep
# EBS salmon code
# 
# updated 2022
# ----------------------------------------
  
    # switches and options:
    #-------------------------------------------
    scaleIN         <-  1      # controls the ratio (relative scaling of window)
    dpiIN           <-  150    # dpi for figures (set to lower res for smaller file size- these will be about 3.5 MB)
    reference_years <-  1970:2000  # for determining the mean for z-scoring
    salmon_years    <-  1976:2020

    # set up directory paths:
    #-------------------------------------------
    ROMS_data_path <- file.path("..","ACLIM2","Data","out")
   
    # Index stuff:
    #-------------------------------------------
    # Juv_summer_strata <-  c(70, 71) 
    # Juv_fall_strata   <- c(70,71, 81, 82, 90)
    # Juv_winter_strata <- c(70,71, 81, 82, 90)
    # Imm_summer_strata <- c(50) #+ offshore (need whole basin is avail )
    # Imm_fall_strata   <- c(61, 62)
    # Imm_winter_strata <- c(61, 62)
 

     varlist <- c("temp_surface5m", 
                      "temp_bottom5m", 
                      "temp_integrated",
                      "EupO_surface5m",
                      "EupO_integrated",
                      "EupS_integrated",
                      "EupS_surface5m",
                      "NCaO_integrated",
                      "NCaO_surface5m",
                      "NCaS_integrated",
                      "NCaS_surface5m")
    # Plotting stuff:
    #-------------------------------------------
    col_ramp <- colorRampPalette(c(brewer.pal(11,"Spectral")[c(3,9)],brewer.pal(11,"RdBu")[11]))
    col      <- col_ramp(7)[c(4,1,6)]
 
    # set up color palettes
    #-------------------------------------------