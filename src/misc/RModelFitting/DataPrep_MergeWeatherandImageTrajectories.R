# Data Cleaning Script
# Cleans, joins, and saves weather data along with image data
# PB 4/24/23
library(tidyverse)

# Set output directories
tabled = '/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/RModelFitting/tables/TempCohensD_ModelResults'
figd = '/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/RModelFitting/figs'
datad = '/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/RModelFitting/data'

# Load in input datasets
esstats_temp =
  read_csv('/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/stats/TTestandEffectSizeStats_Temperature_byExclosure_byImage_AllSites.csv')

esstats_pai =
  read_csv('/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/stats/TTestandEffectSizeStats_PAI_byExclosure_byImage_AllSites.csv')

esstats_height =
  read_csv('/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/stats/TTestandEffectSizeStats_Height_byExclosure_byImage_AllSites.csv')

ImageTraj = 
  read_csv('/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/data/in/ImageTrajectories/AllSites_ImageTrajectory_withWeather.csv')

# Make ImageTraj site names the same as those in the esfiles
# so that you can join them later
ImageTraj['Site'][ImageTraj['Site'] == 'NkuluEP']  <- "Nkuhlu" 
ImageTraj['Site'][ImageTraj['Site'] == 'LetabaExclosure']  <- "Letaba" 
ImageTraj['Site'][ImageTraj['Site'] == 'MakhohlolaEP']  <- "Makhohlola" 
ImageTraj['Site'][ImageTraj['Site'] == 'BuffaloCamp100m']  <- "BuffaloCamp"
ImageTraj['Site'][ImageTraj['Site'] == 'HlangwineEP']  <- "Hlangwine" 

# Rename site and image f cols to match
ImageTraj = ImageTraj %>% rename("site" = "Site",
                                 "imgf" = "Imagef")

# Last, add .tif to image file names so you can join them
namejoin = function(s) paste0(s,'.tif')
ImageTraj['imgf'] = apply(ImageTraj['imgf'], FUN=namejoin, MARGIN=1)

# Make a smaller XY dataframe
# join on image name
pai = subset(esstats_pai,
             select= c(imgf, cohensd, meandiff))
height = subset(esstats_height,
                select= c(imgf, cohensd, meandiff))

XY =  subset(esstats_temp,
             select= c(site, imgf, cohensd, meandiff)) %>%
  rename("temp_cd" = "cohensd", "temp_md" = "meandiff") %>%
  inner_join(pai, by='imgf') %>%
  rename("pai_cd" = "cohensd", "pai_md" = "meandiff") %>%
  inner_join(height, by='imgf') %>%
  rename("height_cd" = "cohensd", "height_md" = "meandiff") %>%
  inner_join(ImageTraj, by= c('site', 'imgf'))

# Make site a factor
XY = XY %>% mutate(site_f = factor(site))

# Make time as a factor
# Note... this doesn't really seem to work, but idk
# timefactor = function(t) strptime(t, '%H:%M:%S')
timefactor = function(t) as.POSIXct(t, format='%H:%M:%S')
XY['Time_f'] = lapply(XY['Time'], FUN=timefactor)

datetimefactor = function(dt) as.POSIXct(dt, "%Y-%m-%d %H:%M:%S")
XY['DateTime_f'] = apply(XY['DateTime'], FUN=datetimefactor, MARGIN=1)

# Also, make a decimal version of time to use in GAM
# See:
# https://stackoverflow.com/questions/5186972/how-to-convert-time-mmss-to-decimal-form-in-r
XY <- XY %>%
  separate(col = Time, into = c("H", "M", "S"), sep = "\\:", remove = FALSE) %>% 
  mutate(H = as.numeric(H)/24) %>% 
  mutate(M = as.numeric(M)/24/60) %>% 
  mutate(S = as.numeric(S)/24/60/60) %>% 
  mutate(Time_num = H+M+S)

# Scale the important vars - Temp, Pai, Height, Time, AirTemp, and humidity 
XY['AirTemp_scaled'] = XY['AirTemp'] %>% scale()
XY['RH_scaled'] = XY['RH'] %>% scale()
XY['temp_cd_scaled'] = XY['temp_cd'] %>% scale()
XY['pai_cd_scaled'] = XY['pai_cd'] %>% scale()
XY['height_cd_scaled'] = XY['height_cd'] %>% scale()
XY['Time_num_scaled'] = XY['Time_num'] %>% scale()

XY['Srad_mJ_scaled'] = XY['Srad_mJ'] %>% scale()
XY['Dewpoint_C_scaled'] = XY['Dewpoint_C'] %>% scale()
XY['Specific_RH_gkg_scaled'] = XY['Specific_RH_gkg'] %>% scale()
XY['Leaf_Wet_Perc_scaled'] = XY['Leaf_Wet_Perc'] %>% scale()
XY['Absolute_RH_gm3_scaled'] = XY['Absolute_RH_gm3'] %>% scale()
XY['Eto_mm_scaled'] = XY['Eto_mm'] %>% scale()

# save out as csv
write_csv(XY, paste0(datad, '/XY.csv'))
