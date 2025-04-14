# Fitting GAMs for Temp differences 
# For Thermal Landscapes
# 3/28/23 
library(tidyverse)
library(mgcv)
library(lme4)
library(gratia)
library(broom)
library(ggplot2)
library(qpcR) # for aikike weights

library("wesanderson")

# Set output directories
# wd = '/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/RModelFitting/tables/initial/TempCohensD_ModelResults'
tabled = '/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/RModelFitting/tables/grapefruit/TempCohensD_ModelResults'
figd = '/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/RModelFitting/figs/grapefruit'
datad = '/n/home02/pbb/scripts/halo-metadata-server/ThermalLandscapes/RModelFitting/data'

XY = read_csv(paste0(datad, '/XY.csv'))

# # # Make Models
# Make site a factor
XY = XY %>% mutate(site_f = factor(site))

# make a vector of formulas
# various combinations of Temp, Pai, Height, Time, AirTemp, and humidity
# grapefruit run model set
formulas = list(as.formula('temp_cd ~ s(pai_cd_scaled) + site_f'),
                as.formula('temp_cd ~ s(height_cd_scaled) + site_f'),
                as.formula('temp_cd ~ s(AirTemp_scaled) + site_f'),
                as.formula('temp_cd ~ s(RH_scaled) + site_f'),
                as.formula('temp_cd ~ s(Time_num_scaled) + site_f'),
                as.formula('temp_cd ~ s(Srad_mJ_scaled) + site_f'),
                as.formula('temp_cd ~ s(Specific_RH_gkg_scaled) + site_f'),
                
                as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + site_f'),
                as.formula('temp_cd ~ s(pai_cd_scaled) + s(RH_scaled) + site_f'),
                as.formula('temp_cd ~ s(pai_cd_scaled) + s(Time_num_scaled) + site_f'),
                as.formula('temp_cd ~ s(pai_cd_scaled) + s(Srad_mJ_scaled) + site_f'),
                as.formula('temp_cd ~ s(pai_cd_scaled) + s(Specific_RH_gkg_scaled) + site_f'),
                
                as.formula('temp_cd ~ s(height_cd_scaled) + s(AirTemp_scaled) + site_f'),
                as.formula('temp_cd ~ s(height_cd_scaled) + s(RH_scaled) + site_f'),
                as.formula('temp_cd ~ s(height_cd_scaled) + s(Time_num_scaled) + site_f'),
                as.formula('temp_cd ~ s(height_cd_scaled) + s(Srad_mJ_scaled) + site_f'),
                as.formula('temp_cd ~ s(height_cd_scaled) + s(Specific_RH_gkg_scaled) + site_f'),
                
                as.formula('temp_cd ~ s(pai_cd_scaled) +
                                      s(RH_scaled) +
                                      s(Srad_mJ_scaled) + site_f'),
                as.formula('temp_cd ~ s(height_cd_scaled) +
                                      s(RH_scaled) +
                                      s(Srad_mJ_scaled) + site_f'),
                as.formula('temp_cd ~ s(RH_scaled) +
                                      s(Srad_mJ_scaled) + site_f'),
                
                as.formula('temp_cd ~ s(pai_cd_scaled) +
                                      s(Time_num_scaled) +
                                      s(Specific_RH_gkg_scaled) + site_f'),
                as.formula('temp_cd ~ s(height_cd_scaled) +
                                      s(Time_num_scaled) +
                                      s(Specific_RH_gkg_scaled) + site_f'),
                as.formula('temp_cd ~ s(Time_num_scaled) +
                                      s(Specific_RH_gkg_scaled) + site_f'),
                
                as.formula('temp_cd ~ s(pai_cd_scaled) +
                                      s(Time_num_scaled) +
                                      s(Specific_RH_gkg_scaled) + site_f'),
                as.formula('temp_cd ~ s(height_cd_scaled) +
                                      s(Time_num_scaled) +
                                      s(Specific_RH_gkg_scaled) + site_f'),
                as.formula('temp_cd ~ s(Time_num_scaled) +
                                      s(Specific_RH_gkg_scaled) + site_f'),
                
                as.formula('temp_cd ~ s(pai_cd_scaled) +
                                      s(AirTemp_scaled) +
                                      s(Srad_mJ_scaled) +
                                      s(Specific_RH_gkg_scaled) + site_f'),
                as.formula('temp_cd ~ s(height_cd_scaled) +
                                      s(AirTemp_scaled) +
                                      s(Srad_mJ_scaled) +
                                      s(Specific_RH_gkg_scaled) + site_f'),
                as.formula('temp_cd ~ s(AirTemp_scaled) +
                                      s(Srad_mJ_scaled) +
                                      s(Specific_RH_gkg_scaled) + site_f')
             )

# Other potential models to add
# as.formula('temp_cd ~ s(pai_cd_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(Specific_RH_gkg_scaled) + site_f'),
# as.formula('temp_cd ~ s(AirTemp_scaled) + site_f'),
# as.formula('temp_cd ~ s(AirTemp_scaled) + s(Specific_RH_gkg_scaled) + site_f'),
# as.formula('temp_cd ~ s(height_cd_scaled) + s(AirTemp_scaled) + site_f'),
# as.formula('temp_cd ~ s(height_cd_scaled) + s(AirTemp_scaled) + s(Specific_RH_gkg_scaled) + site_f'),
#
# -initial model set -
# as.formula('temp_cd ~ s(pai_cd_scaled) + site_f'),
# as.formula('temp_cd ~ s(height_cd_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(height_cd_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(RH_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(Time_num_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(RH_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(RH_scaled) + s(Time_num_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(RH_scaled) + s(Time_num_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(Time_num_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(RH_scaled, by=site_f) + s(Time_num_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(RH_scaled) + s(Time_num_scaled, by=site_f) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled, by=site_f) + s(Time_num_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(Time_num_scaled, by=site_f) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled, RH_scaled) + s(Time_num_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(pai_cd_scaled, AirTemp_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(RH_scaled) + s(pai_cd_scaled, RH_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(pai_cd_scaled, AirTemp_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(pai_cd_scaled, AirTemp_scaled) + s(Time_num_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(RH_scaled) + s(pai_cd_scaled, RH_scaled) + s(Time_num_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(Specific_RH_gkg_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(Srad_mJ_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(Specific_RH_gkg_scaled) + s(Srad_mJ_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(Srad_mJ_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(Specific_RH_gkg_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(Specific_RH_gkg_scaled) + s(Srad_mJ_scaled) + site_f'),
# as.formula('temp_cd ~ s(AirTemp_scaled) + site_f'),
# as.formula('temp_cd ~ s(RH_scaled) + site_f'),
# as.formula('temp_cd ~ s(Specific_RH_gkg_scaled) + site_f'),
# as.formula('temp_cd ~ s(Srad_mJ_scaled) + site_f'),
# as.formula('temp_cd ~ s(Dewpoint_C_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(Dewpoint_C_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled)  + s(AirTemp_scaled) + s(Dewpoint_C_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(Srad_mJ_scaled) + site_f'),
# as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled)')
# 
# Initial Mixed models - add them in later if you need them
# as.formula('temp_cd_scaled ~ s(pai_cd_scaled, by=site_f) + s(AirTemp_scaled) + s(RH_scaled) + site_f'),
# as.formula('temp_cd_scaled ~ s(pai_cd_scaled) + s(AirTemp_scaled, by=site_f) + s(RH_scaled) + site_f'),
# as.formula('temp_cd_scaled ~ s(pai_cd_scaled) + s(AirTemp_scaled) + s(RH_scaled, by=site_f) + site_f'),



# Choose your wes anderson color pallette, v important!
pal = wes_palette("IsleofDogs1", 6, "continuous")

# Initialize List of aic values for weight calculation
aic = c()

# Loop over formulas, fit them, and save results
for (i in seq_along(formulas)){
  
  mod = gam(formulas[[i]],
            data=XY,
            family=gaussian,
            select=FALSE,
            method="REML")
  
  mod_RMSE = sqrt(mean(mod$residuals**2))
  
  # Save important results
  vartbl_mod = tidy(mod)
  write_csv(vartbl_mod,
            paste0(tabled, '/Mod', i, '_VarTable_TempCohensD_.csv'))

  # use sink to save summary outputs in a txt file
  sink(paste0(tabled, '/Mod', i, '_SummaryTable_TempCohensD_.csv'))
  print(paste0("Model: ", i))
  print(summary(mod))
  print(paste0("AIC = ", mod$aic))
  print(paste0("RMSE = ", mod_RMSE))
  sink()
  
  # add aic to list
  aic = aic %>% append(mod$aic)

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
         filename=paste0(figd,'/1to1Plots/Mod', i, '_1to1Plot_TempCohensD_.png'),
         width = 6, height = 5, units = "in", device='png', dpi=400)
  

  partialeffectsplot = draw(mod, residuals = TRUE)
  ggsave(plot = partialeffectsplot,
         filename=paste0(figd,'/PartialEffectsPlots/Mod', i, '_PartialEffectsPlot_TempCohensD.png'),
         width = 6, height = 5, units = "in", device='png', dpi=400)
  appraisalplot = appraise(mod)
  ggsave(plot = appraisalplot,
         filename=paste0(figd,'/AppraisalPlots/Mod', i, '_AppraisalPlot_TempCohensD.png'),
         width = 6, height = 5, units = "in", device='png', dpi=400)
  
}

# Calculate aic weights
aicweights = akaike.weights(aic)

# add model number to df
aicweights = data.frame(aicweights) %>% 
  mutate(modnum=seq_along(formulas))

write_csv(aicweights,
          paste0(tabled, '/AllModels_AkaikeWeights_TempCohensD.csv'))

# # # Fit and plot an important model
# mod = gam(formulas[[26]],
#           data=XY,
#           family=gaussian,
#           select=FALSE,
#           method="REML")
# 
# summary(mod)
# draw(mod, residuals=TRUE)
# 
# mod2 = gam(formula=as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) +
#                     s(Srad_mJ_scaled) + site_f'),
#            data=XY,
#            family=gaussian,
#            select=FALSE,
#            method="REML")
# 
# summary(mod2)
# draw(mod2, residuals=TRUE)
# 
# mod3 = gam(formula=as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled) +
#                                site_f'),
#            data=XY,
#            family=gaussian,
#            select=FALSE,
#            method="REML")
# 
# summary(mod3)
# draw(mod3, residuals=TRUE)
# 
# 
# mod4 = gam(formula=as.formula('temp_cd ~ s(pai_cd_scaled) + s(AirTemp_scaled)'),
#            data=XY,
#            family=gaussian,
#            select=FALSE,
#            method="REML")
# 
# summary(mod4)
# draw(mod4, residuals=TRUE)

# gam.check(mod, rep=500)
# lines(c(-2, 0.5), c(-2, 0.5))

# Make some nice plots of the important variables
p1 =  ggplot(data=XY) + 
    geom_point(aes(x=pai_cd, y=temp_cd, color=site_f))+ 
    theme(axis.text.y = element_text(colour = "black", size = 12),
          axis.text.x = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 12, colour ="black"),
          legend.position = "right",
          title = element_text(face = "bold", size = 12, colour = "black"),
          legend.title = element_text(size = 12, colour = "black", face = "bold"),
          legend.key=element_blank()) +
    theme_bw() +
    xlab('PAI Cohen\'s D') +
    ylab('Temperature Cohen\'s D') +
    coord_fixed(ratio = 1) +
    scale_color_manual(values = pal)
p1
ggsave(plot = p1,
       filename=paste0(figd,'/Scatter_PAICohensD_TempCohensD.png'),
       width = 6, height = 5, units = "in", device='png', dpi=400)

p1 =  ggplot(data=XY) + 
  geom_point(aes(x=AirTemp, y=temp_cd, color=site_f))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right",
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"),
        legend.key=element_blank()) +
  theme_bw() +
  xlab('Air Temperature [C]') +
  ylab('Temperature Cohen\'s D')+
  scale_color_manual(values = pal)
p1
ggsave(plot = p1,
       filename=paste0(figd,'/Scatter_AirTemp_TempCohensD.png'),
       width = 6, height = 5, units = "in", device='png', dpi=400)



p1 =  ggplot(data=XY) + 
  geom_point(aes(x=RH, y=temp_cd, color=site_f))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right",
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"),
        legend.key=element_blank()) +
  theme_bw() +
  xlab('Relative Humidity [%]') +
  ylab('Temperature Cohen\'s D') +
  scale_color_manual(values = pal)
p1
ggsave(plot = p1,
       filename=paste0(figd,'/Scatter_RelHumidity_TempCohensD.png'),
       width = 6, height = 5, units = "in", device='png', dpi=400)


p1 =  ggplot(data=XY) + 
  geom_point(aes(x=Specific_RH_gkg, y=temp_cd, color=site_f))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right",
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"),
        legend.key=element_blank()) +
  theme_bw() +
  xlab('Specific Humidity [g/kg]') +
  ylab('Temperature Cohen\'s D') +
  scale_color_manual(values = pal)
p1
ggsave(plot = p1,
       filename=paste0(figd,'/Scatter_SpecificHumiditygkg_TempCohensD.png'),
       width = 6, height = 5, units = "in", device='png', dpi=400)


p1 =  ggplot(data=XY) + 
  geom_point(aes(x=Srad_mJ, y=temp_cd, color=site_f))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right",
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"),
        legend.key=element_blank()) +
  theme_bw() +
  xlab('Srad_mJ') +
  ylab('Temperature Cohen\'s D') +
  scale_color_manual(values = pal)
p1
ggsave(plot = p1,
       filename=paste0(figd,'/Scatter_SradmJ_TempCohensD.png'),
       width = 6, height = 5, units = "in", device='png', dpi=400)


p1 =  ggplot(data=XY) + 
  geom_point(aes(x=Dewpoint_C, y=temp_cd, color=site_f))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right",
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"),
        legend.key=element_blank()) +
  theme_bw() +
  xlab('Dewpoint_C') +
  ylab('Temperature Cohen\'s D') +
  scale_color_manual(values = pal)
p1
ggsave(plot = p1,
       filename=paste0(figd,'/Scatter_DewpointC_TempCohensD.png'),
       width = 6, height = 5, units = "in", device='png', dpi=400)

p1 =  ggplot(data=XY) + 
  geom_point(aes(x=Time, y=temp_cd, color=site_f))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right",
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"),
        legend.key=element_blank()) +
  theme_bw() +
  xlab('Local Time') +
  ylab('Temperature Cohen\'s D') +
  scale_color_manual(values = pal)
p1
ggsave(plot = p1,
       filename=paste0(figd,'/Scatter_Time_TempCohensD.png'),
       width = 6, height = 5, units = "in", device='png', dpi=400)

