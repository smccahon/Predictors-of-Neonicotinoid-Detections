#---------------------------------#
#  Predictors of Neonics in Water #
#   Analysis by Shelby McCahon    #
#      Created: 2026-01-13        #
#     Modified: 2026-01-13        #
#---------------------------------#

# load packages
library(tidyverse)
library(MuMIn)
library(glmmTMB)
library(DHARMa)
library(AICcmodavg)

# DETAILS OF ANALYSIS
# STAGE I: Identify the best predictors within each hypothesis
# STAGE II: Identify which hypotheses are supported
# STAGE III: Identify combined effects across hypotheses if multiple hypotheses
#            are supported

# Analysis notes
# all detections in seasonal wetlands, so need to combine with temporary

#------------------------------------------------------------------------------#
#                        load data and organize datasets                    ----                        
#------------------------------------------------------------------------------# 

water <- read.csv("cleaned_data/wetland_data_cleaned_2025-09-30.csv")

# theme for plotting
my_theme <- theme_classic() + theme(
  axis.title.x = element_text(size = 21, margin = margin(t = 12)),
  axis.title.y = element_text(size = 21, margin = margin(r = 12)),
  axis.text.x = element_text(size = 18),
  axis.text.y = element_text(size = 18))

options(tibble.print_max = Inf)

#------------------------------------------------------------------------------#
#                              manipulate data                              ----                        
#------------------------------------------------------------------------------# 
# combine temporary and seasonal permanence classes
water <- water %>% 
  mutate(Permanence = case_when(
    Permanence %in% c("Temporary", "Seasonal") ~ "Temporary/Seasonal",
    Permanence == "Semipermanent" ~ "Semipermanent",
    Permanence == "Permanent" ~ "Permanent"))

# convert characters to factors
water <- water %>% 
  mutate_if(is.character, as.factor)

# filter data (n = 97)
water <- water %>%
  filter(!is.na(WaterNeonicDetection))

# standardize data
water.cs <- water %>%
  mutate(across(where(is.numeric), scale))

# set reference levels
water.cs$DominantCrop <- relevel(water.cs$DominantCrop,
                                 ref = "Grassland")

water.cs$Permanence <- relevel(water.cs$Permanence,
                               ref = "Temporary/Seasonal")

#------------------------------------------------------------------------------#
#                   identify best model for each hypothesis                 ----                        
#------------------------------------------------------------------------------# 

# include all variables that show up in models with support (delta AICc < 2)

# best agricultural model
global.model <- glm(WaterNeonicDetection ~ PercentAg + DominantCrop +
                      NearestCropDistance_m + Buffered,
                    family = "binomial",
                    data = water.cs,
                    na.action = na.fail)

# check for collinearity
car::vif(global.model) # vif < 3

dredge(global.model)

# Model selection table 
#      (Int) Bff DmC    NCD_m    PrA df  logLik  AICc delta weight
# 14 -0.5968   +      1.44400 0.7345  4 -48.093 104.6  0.00  0.413
# 6  -0.3406   +      1.18400         3 -49.871 106.0  1.38  0.207
# 13 -1.2420          0.70910 1.0490  3 -50.107 106.5  1.85  0.164
# 9  -1.1870                  0.5490  2 -51.908 107.9  3.32  0.078
# 2  -0.8675   +                      2 -52.412 109.0  4.33  0.047
# 10 -1.0660   +              0.4101  3 -51.749 109.8  5.14  0.032
# 1  -1.1120                          1 -54.270 110.6  5.96  0.021
# 8  -0.7906   +   +  1.31400         7 -48.339 111.9  7.32  0.011
# 5  -1.1130         -0.05361         2 -54.246 112.6  8.00  0.008
# 16 -0.7463   +   +  1.44300 0.5701  8 -47.620 112.9  8.26  0.007
# 3  -1.6090       +                  5 -51.689 114.0  9.42  0.004
# 15 -1.4570       +  0.67940 0.8578  7 -49.599 114.5  9.84  0.003
# 4  -1.2990   +   +                  6 -51.206 115.3 10.72  0.002
# 11 -1.4290       +          0.3394  6 -51.241 115.4 10.80  0.002
# 7  -1.7280       +  0.23580         6 -51.337 115.6 10.99  0.002
# 12 -1.2920   +   +          0.2080  7 -51.091 117.4 12.82  0.001
# Models ranked by AICc(x) 

# can remove dominant crop
m1 <- glm(WaterNeonicDetection ~ PercentAg + NearestCropDistance_m + Buffered,
          family = "binomial",
          data = water.cs)

summary(m1)
confint(m1)

# model validation --> no big issues
simulationOutput <- simulateResiduals(fittedModel = m1) 
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m1)$PercentAg) #  good
plotResiduals(simulationOutput, form = model.frame(m1)$NearestCropDistance_m) # some pattern
plotResiduals(simulationOutput, form = model.frame(m1)$Buffered)  # good

#------------------------------------------------------------------------------#

# best temporal model
m2 <- glm(WaterNeonicDetection ~ Event,
                    data = water.cs,
                    family = "binomial",
                    na.action = na.fail)

summary(m2)
confint(m2)

# model validation --> no issues
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$Event) #  good

#------------------------------------------------------------------------------#

# best wetland dynamics model
global.model <- glm(WaterNeonicDetection ~ Permanence + Dist_Closest_Wetland_m,
                    family = "binomial",
                    na.action = na.fail,
                    data = water.cs)

# check for collinearity
car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#     (Int) Dst_Cls_Wtl_m Prm df  logLik  AICc delta weight
# 3  0.5465                 +  4 -35.755  79.9  0.00  0.654
# 4  0.5342       -0.2978   +  5 -35.279  81.2  1.27  0.346
# 1 -1.1120                    1 -54.270 110.6 30.64  0.000
# 2 -1.1330       -0.2869      2 -53.724 111.6 31.63  0.000
# Models ranked by AICc(x) 

# retain all predictors
m3 <- glm(WaterNeonicDetection ~ Permanence + Dist_Closest_Wetland_m,
          data = water.cs,
          family = "binomial",
          na.action = na.fail)

summary(m3)
confint(m3)

#------------------------------------------------------------------------------#

# best precipitation model
global.model <- glm(WaterNeonicDetection ~ AnnualSnowfall_in + 
                      AnnualPrecipitation_in + DaysSinceLastPrecipitation_5mm +
                      PrecipitationAmount_7days,
                    family = "binomial",
                    na.action = na.fail,
                    data = water.cs)

# check for collinearity
car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#     (Int)   AnP_in  AnS_in DSL_5mm   PrA_7dy df  logLik  AICc delta weight
# 1  -1.112                                     1 -54.270 110.6  0.00  0.199
# 3  -1.134          -0.2925                    2 -53.505 111.1  0.56  0.151
# 5  -1.129                  -0.2598            2 -53.835 111.8  1.22  0.108
# 7  -1.153          -0.2988 -0.2773            3 -53.045 112.3  1.77  0.082
# 9  -1.113                           0.061350  2 -54.235 112.6  2.02  0.073
# 2  -1.112 -0.01712                            2 -54.268 112.7  2.08  0.070
# 4  -1.136  0.11490 -0.3374                    3 -53.405 113.1  2.48  0.057
# 11 -1.134          -0.2896          0.039110  3 -53.491 113.2  2.66  0.053
# 13 -1.132                  -0.3125 -0.075440  3 -53.799 113.9  3.27  0.039
# 6  -1.129 -0.03454         -0.2640            3 -53.824 113.9  3.32  0.038
# 15 -1.159          -0.3093 -0.3640 -0.112500  4 -52.963 114.4  3.78  0.030
# 8  -1.154  0.09180 -0.3335 -0.2652            4 -52.980 114.4  3.81  0.030
# 10 -1.114 -0.04647                  0.078260  3 -54.219 114.7  4.11  0.025
# 12 -1.136  0.11970 -0.3400         -0.009243  4 -53.404 115.2  4.66  0.019
# 14 -1.132 -0.01186         -0.3102 -0.070060  4 -53.798 116.0  5.45  0.013
# 16 -1.166  0.18150 -0.3875 -0.4126 -0.205000  5 -52.768 116.2  5.61  0.012
# Models ranked by AICc(x) 

# can remove annual precipitation & precipitation amount (last 7 days)

m4 <- glm(WaterNeonicDetection ~ AnnualSnowfall_in + 
                      DaysSinceLastPrecipitation_5mm,
                    family = "binomial",
                    na.action = na.fail,
                    data = water.cs)

summary(m4)
confint(m4)

        