# Fitting GAMs for Temp differences 
# For Thermal Landscapes
# 3/28/23 
library(tidyverse)
library(mgcv)
library(lme4)
library(gratia)
library(broom)
library(ggplot2)

# Set output directories
tabled = '/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/RModelFitting/tables/TempCohensD_ModelResults'
figd = '/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/RModelFitting/figs/1to1Plots'

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


# # # Plots
# 
# ggplot(data = XY) + 
#   geom_point(aes(x=pai_cd, y=temp_cd, color=site))
# 
# ggplot(data = XY) + 
#   geom_point(aes(x=height_cd, y=temp_cd, color=site))
# 
# ggplot(data = XY) + 
#   geom_point(aes(x=height_cd, y=temp_cd, color=RH))
# 
# ggplot(data = XY) + 
#   geom_point(aes(x=height_cd, y=temp_cd, color=AirTemp))
# 
# ggplot(data = XY) + 
#   geom_point(aes(x=height_cd, y=temp_cd, color=Time_f))
# 
# ggplot(data = XY) + 
#   geom_point(aes(x=height_cd, y=temp_cd, color=DateTime_f))
# 
# ggplot(data = XY) + 
#   geom_point(aes(x=RH_scaled, y=temp_cd, color=site)) + scale_x_time()
# 
# ggplot(data = XY) + 
#   geom_point(aes(x=AirTemp_scaled, y=temp_cd, color=site)) + scale_x_time()
# 
# ggplot(data = XY) + 
#   geom_point(aes(x=Time, y=height_cd, color=site)) + scale_x_time()
# 
# ggplot(data = XY) + 
#   geom_point(aes(x=DateTime_f, y=height_cd, color=site)) + scale_x_time()


# # # Make Models

# make a vector of formulas
# various combinations of Temp, Pai, Height, Time, AirTemp, and humidity

formulas = list(as.formula('temp_cd ~ s(pai_cd_scaled) + site_f'),
             as.formula('temp_cd ~ s(height_cd_scaled) + site_f'),
             as.formula('temp_cd ~ s(pai_cd_scaled) + s(height_cd_scaled) + site_f'),
             as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + site_f'),
             as.formula('temp_cd ~ s(pai_cd_scaled) + s(RH_scaled) + site_f'),
             as.formula('temp_cd ~ s(pai_cd_scaled) + s(Time_num_scaled) + site_f'),
             as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(RH_scaled) + site_f'),
             as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(RH_scaled) + s(Time_num_scaled) + site_f'),
             as.formula('temp_cd ~ s(pai_cd_scaled) + s(RH_scaled) + s(Time_num_scaled) + site_f'),
             as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(Time_num_scaled) + site_f')
             )

# Mixed models - add them in later if you need them
# as.formula('temp_cd_scaled ~ s(pai_cd_scaled, by=site_f) + s(AirTemp_scaled) + s(RH_scaled) + site_f'),
# as.formula('temp_cd_scaled ~ s(pai_cd_scaled) + s(AirTemp_scaled, by=site_f) + s(RH_scaled) + site_f'),
# as.formula('temp_cd_scaled ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(RH_scaled, by=site_f) + site_f'),

library("wesanderson")
# Choose your wes anderson color pallette, v important!
pal = wes_palette("IsleofDogs1", 6, "continuous")

# Loop over formulas, fit them, and save results
for (i in seq_along(formulas)){
  
  mod = gam(formulas[[i]],
            data=XY,
            family=gaussian,
            select=TRUE,
            method="REML")
  
  mod_RMSE = sqrt(mean(mod$residuals**2))
  
  # Save important results
  vartbl_mod = tidy(mod)
  write_csv(vartbl_mod,
            paste0(tabled, '/Mod', i, '_VarTable_TempCohensD_.csv'))

  # use sink to save summary outputs in a txt file
  sink(paste0(tabled, '/Mod', i, '_SummaryTable_TempCohensD_.csv'))
  print(summary(mod))
  print(paste0("AIC = ", mod$aic))
  print(paste0("RMSE = ", mod_RMSE))
  sink()

  # # # Save out a 1-to-1 plot for each model
  
  dfplot = data.frame(temp_cd_fitted = predict(mod)) %>% 
           mutate(temp_cd = XY$temp_cd,
                  Site = XY$site_f,
                  DateTime = XY$DateTime_f,
                  Time = XY$Time_f)
  
  
  p1to1 =  ggplot(data=dfplot) +
    geom_point(mapping = aes(x = temp_cd_fitted,
                             y = temp_cd,
                             colour = Site),
               size=2.2,
               alpha=0.85) +
    geom_abline()  +
    theme_bw() +
    theme(axis.text.y = element_text(colour = "black", size = 12),
          axis.text.x = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 12, colour ="black"),
          legend.position = "right",
          title = element_text(face = "bold", size = 12, colour = "black"),
          legend.title = element_text(size = 12, colour = "black", face = "bold"),
          legend.key=element_blank()) +
    xlab('Fitted') +
    ylab('Observed') +
    coord_fixed(ratio = 1) +
    scale_color_manual(values = pal)
  
  p1to1
  
  ggsave(plot = p1to1,
         filename=paste0(figd,'/Mod', i, '_1to1Plot_TempCohensD_.png'),
         width = 6, height = 5, units = "in", device='png', dpi=400)
  
  # # # TBD - save another set of plots with time added
  # 
  # p1to1 =  ggplot(data=dfplot) +
  #   geom_point(mapping = aes(x = temp_cd_fitted,
  #                            y = temp_cd,
  #                            shape = Site,
  #                            color = Time),
  #              size=3) +
  #   geom_abline()  +
  #   theme(axis.text.y = element_text(colour = "black", size = 12),
  #         axis.text.x = element_text(colour = "black", size = 12),
  #         legend.text = element_text(size = 12, colour ="black"),
  #         legend.position = "right",
  #         title = element_text(face = "bold", size = 12, colour = "black"),
  #         legend.title = element_text(size = 12, colour = "black", face = "bold"),
  #         legend.key=element_blank()) +
  #   xlab('Fitted') +
  #   ylab('Observed') +
  #   labs(shape='Site', colour='Time') +
  #   coord_fixed(ratio = 1)
  # # +
  # # coord_fixed(xlim=c(3, 3.8),
  # #             ylim=c(3, 3.8))
  # 
  # # scale_colour_manual(values = c("grey30", "#D55E00")) +
  # p1to1
  # 
  # ggsave(plot = p1to1,
  #        filename=paste0(figd,'/Mod', i, '_1to1Plot_TempCohensD_withTime.png'),
  #        width = 9, height = 6, units = "in", device='png', dpi=300)
  
}


summary(mod)$dev.expl
