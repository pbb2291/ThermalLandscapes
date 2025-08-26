
# Fit a model of temp cd based on pai, height, and site
tempcd_m1 = gam(temp_cd ~ s(pai_cd) + site_f,
                data=XY,
                family=gaussian,
                select=TRUE,
                method="REML")


summary(tempcd_m1)
draw(tempcd_m1, residuals = TRUE)
appraise(tempcd_m1)
gam.check(tempcd_m1, rep=500)
lines(c(-2, 0.5), c(-2, 0.5))

# Try adding a random slope as well
tempcd_m2 = gam(temp_cd ~ s(pai_cd, by=site_f) + site_f,
                data=XY,
                family=gaussian,
                select=TRUE,
                method="REML")


summary(tempcd_m2)
draw(tempcd_m2, residuals = TRUE)
appraise(tempcd_m2)
gam.check(tempcd_m2, rep=500)
lines(c(-2, 0.5), c(-2, 0.5))

# Incorporate heights as well
tempcd_m3 = gam(temp_cd ~ s(pai_cd, by=site_f) + s(height_cd) + site_f,
                data=XY,
                family=gaussian,
                select=TRUE,
                method="REML")


summary(tempcd_m3)
draw(tempcd_m3, residuals = TRUE)
appraise(tempcd_m3)
gam.check(tempcd_m3, rep=500)
lines(c(-2, 0.5), c(-2, 0.5))

# # # WEATHER MODELS


tempcd_weather1 = gam(temp_cd ~ s(height_cd) + s(pai_cd) + s(AirTemp_scaled) + s(RH_scaled) + site_f,
                      data=XY,
                      family=gaussian,
                      select=TRUE,
                      method="REML")


summary(tempcd_weather1)
draw(tempcd_weather1, residuals = TRUE)
appraise(tempcd_weather1)
gam.check(tempcd_weather1, rep=500)
lines(c(-2, 0.5), c(-2, 0.5))


tempcd_weather1.1 = gam(temp_cd ~ s(pai_cd) + s(RH_scaled) + site_f,
                        data=XY,
                        family=gaussian,
                        select=TRUE,
                        method="REML")


summary(tempcd_weather1.1)
draw(tempcd_weather1.1, residuals = TRUE)
appraise(tempcd_weather1.1)
gam.check(tempcd_weather1.1, rep=500)
lines(c(-2, 0.5), c(-2, 0.5))

# # # try mixed model 

tempcd_weather2.0 = gam(temp_cd ~ s(pai_cd) + s(AirTemp_scaled)+ s(RH_scaled, by=site_f) + site_f,
                        data=XY,
                        family=gaussian,
                        select=TRUE,
                        method="REML")


summary(tempcd_weather.0)
draw(tempcd_weather2.0, residuals = TRUE)
appraise(tempcd_weather2.0)
gam.check(tempcd_weather2.0, rep=500)
lines(c(-2, 0.5), c(-2, 0.5))

tempcd_weather2.1 = gam(temp_cd ~ s(pai_cd) + s(RH_scaled, by=site_f) + site_f,
                        data=XY,
                        family=gaussian,
                        select=TRUE,
                        method="REML")


summary(tempcd_weather2.1)
draw(tempcd_weather2, residuals = TRUE)
appraise(tempcd_weather2)
gam.check(tempcd_weather2, rep=500)
lines(c(-2, 0.5), c(-2, 0.5))

# # # mixed model

tempcd_weather3.0 = gam(temp_cd ~ s(pai_cd) + s(AirTemp_scaled, by=site_f) + site_f,
                        data=XY,
                        family=gaussian,
                        select=TRUE,
                        method="REML")

summary(tempcd_weather3.0)
draw(tempcd_weather3.0, residuals = TRUE)
appraise(tempcd_weather3.0)
gam.check(tempcd_weather3.0, rep=500)
lines(c(-2, 0.5), c(-2, 0.5))

# Time model

tempcd_weather4.0 = gam(temp_cd ~ s(pai_cd) + s(Time_num) + site_f,
                        data=XY,
                        family=gaussian,
                        select=TRUE,
                        method="REML")

summary(tempcd_weather4.0)
draw(tempcd_weather4.0, residuals = TRUE)
appraise(tempcd_weather4.0)
gam.check(tempcd_weather4.0, rep=500)
lines(c(-2, 0.5), c(-2, 0.5))

# Time is significant, but does it stand up to Temp?
tempcd_weather4.1 = gam(temp_cd ~ s(pai_cd) + s(RH_scaled) + s(Time_num) + site_f,
                        data=XY,
                        family=gaussian,
                        select=TRUE,
                        method="REML")

summary(tempcd_weather4.1)
draw(tempcd_weather4.1, residuals = TRUE)
appraise(tempcd_weather4.1)
gam.check(tempcd_weather4.1, rep=500)
lines(c(-2, 0.5), c(-2, 0.5))

# Time is significant, but does it stand up to RH? 
tempcd_weather4.2 = gam(temp_cd ~ s(pai_cd) + s(AirTemp_scaled) +
                          s(RH_scaled) + s(Time_num) + site_f,
                        data=XY,
                        family=gaussian,
                        select=TRUE,
                        method="REML")

summary(tempcd_weather4.2)
draw(tempcd_weather4.2, residuals = TRUE)
appraise(tempcd_weather4.2)
gam.check(tempcd_weather4.2, rep=500)
lines(c(-2, 0.5), c(-2, 0.5))



tidy(tempcd_weather1.1)
tidy(tempcd_weather2.1)
tidy(tempcd_weather3.0)
tidy(tempcd_weather4.0)
tidy(tempcd_weather4.2)
tempcd_weather4.2$aic
tempcd_weather4.0$aic
# best:
tempcd_weather3.0$aic # pai + air temp mixed
tempcd_weather2.1$aic # pai + RH mixed