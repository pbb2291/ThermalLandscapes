# Correlations and Simple Models of Temperature 
# Following Miriam's suggestion for thermal landscapes paper
# PB
# 4/24/23
library(tidyverse)
library(mgcv)
library(lme4)
library(gratia)
library(broom)
library(ggplot2)

# Set output directories
tabled = '/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/RModelFitting/tables/'
figd = '/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/RModelFitting/figs'
datad = '/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/RModelFitting/data'

XY = read_csv(paste0(datad, '/XY.csv'))

# Note: add weather variable that miriam was talking about - vapor pressure
XY_subset = XY %>% select(c('temp_cd', 'pai_cd', 'height_cd', 'AirTemp',
                            'RH', 'Time_num', 'Altitude', "Srad_mJ",
                            "Dewpoint_C", "Specific_RH_gkg", "Leaf_Wet_Perc",
                            "Absolute_RH_gm3"))
# , "Eto_mm"
# https://www.guru99.com/r-pearson-spearman-correlation.html
corrmat = data.frame(cor(XY_subset, method="spearman"))

# as.dist(round(corrmat, 2))

write_csv(corrmat, paste0(tabled, 'TempCohensD_SpearmanCorrMatrix.csv'))


# Fit a mixed mod
mod = glm(temp_cd ~ pai_cd_scaled + AirTemp_scaled + RH_scaled, data=XY)
plot(mod)
summary(mod)

# # Fit a mixed mod
# mod2 = glm(temp_cd ~ pai_cd_scaled + AirTemp_scaled + site_f, data=XY)
# # plot(mod)
# summary(mod2)
# mod2$

# pai_cd_scaled:height_cd_scaled

