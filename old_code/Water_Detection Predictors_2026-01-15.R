#---------------------------------#
#  Predictors of Neonics in Water #
#   Analysis by Shelby McCahon    #
#      Created: 2026-01-15        #
#     Modified: 2026-01-21        #
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
# STAGE III: Model-average top predictors (95% confidence set)

# Analysis notes
# seasonal wetlands all had detections, so I combined it with temporary wetlands
# event highly influential so I included it in all models and used it as informed null
# results do not change depending on if I include correlated agricultural 
# variables together or not. log(cropland distance) is still the only informative
# metric (no other variables are within 2 delta AIC)

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

# convert factors to binary for correlation testing
 # water <- water %>%
 #   mutate(Season = ifelse(Season == "Spring", 1, 0)) %>% 
 #   mutate(Buffered = ifelse(Buffered == "Y", 1, 0))
 # 
 # water$Year <- as.factor(water$Year)

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

# log transform nearest crop distance due to skew
water <- water %>% 
  mutate(LogCropDistance = NearestCropDistance_m + 0.0001) %>% 
  mutate(LogCropDistance = log(LogCropDistance))

# log transform precipitation amount due to skew
water <- water %>% 
  mutate(LogPrecipitationAmount_7days = PrecipitationAmount_7days + 0.0001) %>% 
  mutate(LogPrecipAmount = log(LogPrecipitationAmount_7days))

# log transform nearest crop distance due to skew
water <- water %>% 
  mutate(LogDaysSinceLastPrecipitation_5mm = DaysSinceLastPrecipitation_5mm + 0.0001) %>% 
  mutate(LogPrecipDays = log(LogDaysSinceLastPrecipitation_5mm))

# log transform nearest crop distance due to skew
water <- water %>% 
   mutate(LogWetlandDistance = log(Dist_Closest_Wetland_m))

# standardize data
water.cs <- water %>%
  mutate(across(where(is.numeric), scale))

# set reference levels
water.cs$DominantCrop <- relevel(water.cs$DominantCrop,
                                 ref = "Grassland")

water.cs$Permanence <- relevel(water.cs$Permanence,
                               ref = "Temporary/Seasonal")

#------------------------------------------------------------------------------#
#                         identify correlations                             ----                        
#------------------------------------------------------------------------------# 

# cor_mat <- cor(water[sapply(water, is.numeric)], 
#                use = "pairwise.complete.obs")
# 
# # Convert matrix to a data frame of pairs
# cor_df <- as.data.frame(as.table(cor_mat))
# 
# # Rename columns
# names(cor_df) <- c("Var1", "Var2", "Correlation")
# 
# # Convert factors to characters
# cor_df$Var1 <- as.character(cor_df$Var1)
# cor_df$Var2 <- as.character(cor_df$Var2)
# 
# # Keep only unique pairs (upper triangle)
# cor_df <- cor_df[cor_df$Var1 < cor_df$Var2, ]
# 
# # Filter by threshold
# high_corr <- subset(cor_df, abs(Correlation) > 0.6)
# 
# # Sort by correlation strength
# high_corr <- high_corr[order(-abs(high_corr$Correlation)), ]
# 
# high_corr

# relevant correlations
  #                      Var1                  Var2 Correlation
 # 21                 Julian                Season  -0.9736649
 # 82                 Season                  SPEI   0.8311283
 # 81                 Julian                  SPEI  -0.7553440
 # 66               Buffered NearestCropDistance_m   0.7409913
 # 46               Buffered             PercentAg  -0.7367155
 # 44  NearestCropDistance_m             PercentAg  -0.7150703
 
#------------------------------------------------------------------------------#
#                   identify best model for each hypothesis                 ----                        
#------------------------------------------------------------------------------# 

# include all variables that show up in models with support (delta AICc < 2)

# log-transformation? YES
# m1 <- glm(WaterNeonicDetection ~ LogCropDistance + Event,
#           family = "binomial",
#           data = water.cs,
#           na.action = na.fail)
# 
# m2 <- glm(WaterNeonicDetection ~ NearestCropDistance_m + Event,
#           family = "binomial",
#           data = water.cs,
#           na.action = na.fail)
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)

# which ag variable should I use? answer: nearest crop distance
# m1 <- glm(WaterNeonicDetection ~ PercentAg + DominantCrop + Event,
#                      family = "binomial",
#                      data = water.cs,
#                      na.action = na.fail)
# 
# m2 <- glm(WaterNeonicDetection ~ Buffered + DominantCrop + Event,
#             family = "binomial",
#             data = water.cs,
#             na.action = na.fail)
#   
# # best model if considering correlations
# m3 <- glm(WaterNeonicDetection ~ LogCropDistance + DominantCrop + Event,
#             family = "binomial",
#             data = water.cs,
#             na.action = na.fail)

# model_names <- paste0("m", 1:3)
# models <- mget(model_names)
# aictab(models, modnames = model_names)


# best agricultural model (considering correlations)
# global.model <- glm(WaterNeonicDetection ~ DominantCrop +
#                       LogCropDistance + Event,
#                     family = "binomial",
#                     data = water.cs,
#                     na.action = na.fail)

# considering all variables (even correlated ones)
global.model <- glm(WaterNeonicDetection ~ LogCropDistance + Event + 
                      PercentAg + Buffered + DominantCrop,
                    family = "binomial",
                    data = water.cs,
                    na.action = na.fail)

# check for collinearity
car::vif(global.model) # vif < 3

dredge(global.model)

# Model selection table 
#    (Intrc) Bffrd DmnnC Event   LgCrD    PrcnA df  logLik  AICc delta weight
# 13 -0.3501                 + -0.6987           5 -37.449  85.6  0.00  0.459
# 14 -0.3539     +           + -0.7872           6 -37.328  87.6  2.03  0.166
# 29 -0.3326                 + -0.7223 -0.03984  6 -37.444  87.8  2.26  0.148
# 21 -0.7844                 +          0.49630  5 -39.601  89.9  4.30  0.053
# 30 -0.3945     +           + -0.7544  0.09225  7 -37.310  89.9  4.32  0.053
# 5  -0.6931                 +                   4 -41.112  90.7  5.10  0.036
# 6  -0.6171     +           +                   5 -40.498  91.7  6.10  0.022
# 22 -0.8041     +           +          0.53500  6 -39.591  92.1  6.56  0.017
# 9  -1.3230                   -1.0200           2 -44.347  92.8  7.26  0.012
# 7  -0.8509           +     +                   8 -38.051  93.7  8.18  0.008
# 15 -0.3503           +     + -0.6973           9 -37.215  94.5  8.94  0.005
# 25 -1.3160                   -1.1090 -0.17630  3 -44.223  94.7  9.15  0.005
# 10 -1.3740     +             -1.0570           3 -44.310  94.9  9.32  0.004
# 23 -0.8471           +     +          0.08525  9 -38.028  96.1 10.57  0.002
# 8  -0.8612     +     +     +                   9 -38.048  96.2 10.61  0.002
# 16 -0.3558     +     +     + -0.7853          10 -37.097  96.8 11.20  0.002
# 26 -1.3130     +             -1.1080 -0.18000  4 -44.223  96.9 11.32  0.002
# 31 -0.3438           +     + -0.7072 -0.03120 10 -37.213  97.0 11.43  0.002
# 11 -0.6509           +       -1.6300           6 -42.593  98.1 12.56  0.001
# 24 -0.8862     +     +     +          0.15790 10 -37.998  98.6 13.00  0.001
# 32 -0.3824     +     +     + -0.7734  0.11730 11 -37.070  99.2 13.69  0.000
# 27 -0.6246           +       -1.6080  0.08891  7 -42.567 100.4 14.84  0.000
# 12 -0.6324     +     +       -1.6190           7 -42.589 100.4 14.88  0.000
# 28 -0.6314     +     +       -1.6110  0.09733  8 -42.567 102.8 17.21  0.000
# 17 -1.1870                            0.54900  2 -51.908 107.9 22.39  0.000
# 2  -0.8675     +                               2 -52.412 109.0 23.40  0.000
# 18 -1.0660     +                      0.41010  3 -51.749 109.8 24.20  0.000
# 1  -1.1120                                     1 -54.270 110.6 25.03  0.000
# 3  -1.6090           +                         5 -51.689 114.0 28.48  0.000
# 4  -1.2990     +     +                         6 -51.206 115.3 29.79  0.000
# 19 -1.4290           +                0.33940  6 -51.241 115.4 29.86  0.000
# 20 -1.2920     +     +                0.20800  7 -51.091 117.4 31.88  0.000
# Models ranked by AICc(x) 

# can remove dominant crop and percent ag
# ** correlated variables didn't end up in there anyways
m1 <- glm(WaterNeonicDetection ~ LogCropDistance + Event,
          family = "binomial",
          data = water.cs)

car::vif(m1) # vif < 2

summary(m1)
confint(m1)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m1, n = 1000) 
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m1)$LogCropDistance) # good
plotResiduals(simulationOutput, form = model.frame(m1)$Event) # good

#------------------------------------------------------------------------------#

# best wetland dynamics model

# log-transformation? NO, although they perform similarily (wt = 0.46)
# does it affect end results? NO (variable not a top predictor anyways)
# m1 <- glm(WaterNeonicDetection ~ LogWetlandDistance + Event,
#           family = "binomial",
#           data = water.cs,
#           na.action = na.fail)
# 
# m2 <- glm(WaterNeonicDetection ~ Dist_Closest_Wetland_m + Event,
#           family = "binomial",
#           data = water.cs,
#           na.action = na.fail)
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)


global.model <- glm(WaterNeonicDetection ~ Permanence + 
                      Dist_Closest_Wetland_m + Event + Porosity,
                    family = "binomial",
                    na.action = na.fail,
                    data = water.cs)

# check for collinearity
car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#     (Int) Dst_Cls_Wtl_m Evn Prm    Prs df  logLik  AICc delta weight
# 7   1.1090                 +   +         6 -38.247  89.4  0.00  0.312
# 3  -0.6931                 +             4 -41.112  90.7  1.23  0.169
# 15  1.0380                 +   + 0.3073  7 -37.987  91.2  1.80  0.126
# 8   1.0010       -0.1591   +   +         7 -38.101  91.5  2.03  0.113
# 11 -0.7026                 +     0.4430  5 -40.402  91.5  2.04  0.113
# 4  -0.7623       -0.2115   +             5 -40.835  92.3  2.90  0.073
# 12 -0.7718       -0.1962   +     0.4339  6 -40.167  93.3  3.84  0.046
# 16  0.9313       -0.1495   +   + 0.2989  8 -37.856  93.3  3.92  0.044
# 13 -0.6453                     + 0.6399  4 -45.780 100.0 10.57  0.002
# 14 -0.6832       -0.3636       + 0.6761  5 -44.945 100.5 11.12  0.001
# 5  -0.3514                     +         3 -47.226 100.7 11.28  0.001
# 6  -0.3663       -0.3415       +         4 -46.522 101.5 12.05  0.001
# 9  -1.3040                       0.8988  2 -49.985 104.1 14.67  0.000
# 10 -1.3350       -0.3728         0.9464  3 -49.043 104.3 14.92  0.000
# 1  -1.1120                               1 -54.270 110.6 21.15  0.000
# 2  -1.1330       -0.2869                 2 -53.724 111.6 22.15  0.000
# Models ranked by AICc(x) 

# distance to nearest wetland can be removed
m2 <- glm(WaterNeonicDetection ~ Permanence + Event + Porosity,
          data = water.cs,
          family = "binomial",
          na.action = na.fail)

summary(m2)
confint(m2)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m2, n = 1000) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$Permanence) #  good
plotResiduals(simulationOutput, form = model.frame(m2)$Event) # good
plotResiduals(simulationOutput, form = model.frame(m2)$Porosity) # pattern


#------------------------------------------------------------------------------#

# log-transformation? no
 # m1 <- glm(WaterNeonicDetection ~ LogPrecipDays + Event,
 #           family = "binomial",
 #           data = water.cs,
 #           na.action = na.fail)
 # 
 # m2 <- glm(WaterNeonicDetection ~ DaysSinceLastPrecipitation_5mm + Event,
 #           family = "binomial",
 #           data = water.cs,
 #           na.action = na.fail)
 # 
 # model_names <- paste0("m", 1:2)
 # models <- mget(model_names)
 # aictab(models, modnames = model_names)
 # 
 # # log-transformation? no
 # m1 <- glm(WaterNeonicDetection ~ LogPrecipAmount + Event,
 #           family = "binomial",
 #           data = water.cs,
 #           na.action = na.fail)
 # 
 # m2 <- glm(WaterNeonicDetection ~ PrecipitationAmount_7days + Event,
 #           family = "binomial",
 #           data = water.cs,
 #           na.action = na.fail)
 # 
 # model_names <- paste0("m", 1:2)
 # models <- mget(model_names)
 # aictab(models, modnames = model_names)

# best precipitation model
global.model <- glm(WaterNeonicDetection ~ AnnualSnowfall_in + 
                      AnnualPrecipitation_in + DaysSinceLastPrecipitation_5mm +
                      PrecipitationAmount_7days + Event,
                    family = "binomial",
                    na.action = na.fail,
                    data = water.cs)

# check for collinearity
car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#      (Int)   AnP_in  AnS_in DSL_5mm Evn   PrA_7dy df  logLik  AICc delta weight
# 9  -0.6931                            +            4 -41.112  90.7  0.00  0.120
# 10 -1.0650 -0.50470                   +            5 -40.130  90.9  0.26  0.106
# 11 -1.6890          -0.5883           +            5 -40.134  90.9  0.27  0.105
# 26 -1.4470 -0.66470                   + -0.506000  6 -39.319  91.6  0.91  0.076
# 27 -2.2100          -0.7427           + -0.470700  6 -39.389  91.7  1.05  0.071
# 30 -1.4950 -0.65150         -0.5013   + -0.813800  7 -38.278  91.8  1.16  0.068
# 31 -2.2450          -0.7332 -0.4669   + -0.761600  7 -38.326  91.9  1.25  0.064
# 13 -0.6369                  -0.2568   +            5 -40.655  92.0  1.31  0.063
# 25 -0.8481                            + -0.287800  5 -40.792  92.2  1.59  0.055
# 29 -0.8984                  -0.5059   + -0.599300  6 -39.663  92.3  1.60  0.054
# 15 -1.5730          -0.5460 -0.2050   +            6 -39.829  92.6  1.93  0.046
# 14 -0.9936 -0.46090         -0.2008   +            6 -39.857  92.6  1.99  0.045
# 12 -1.5680 -0.32280 -0.3769           +            6 -39.866  92.7  2.01  0.044
# 28 -2.1210 -0.44910 -0.4741           + -0.555000  7 -38.924  93.1  2.45  0.035
# 32 -2.1650 -0.43950 -0.4706 -0.4781   + -0.852400  8 -37.890  93.4  2.76  0.030
# 16 -1.4780 -0.28960 -0.3608 -0.1892   +            7 -39.614  94.5  3.83  0.018
# 1  -1.1120                                         1 -54.270 110.6 19.92  0.000
# 3  -1.1340          -0.2925                        2 -53.505 111.1 20.48  0.000
# 5  -1.1290                  -0.2598                2 -53.835 111.8 21.14  0.000
# 7  -1.1530          -0.2988 -0.2773                3 -53.045 112.3 21.69  0.000
# 17 -1.1130                               0.061350  2 -54.235 112.6 21.94  0.000
# 2  -1.1120 -0.01712                                2 -54.268 112.7 22.00  0.000
# 4  -1.1360  0.11490 -0.3374                        3 -53.405 113.1 22.41  0.000
# 19 -1.1340          -0.2896              0.039110  3 -53.491 113.2 22.58  0.000
# 21 -1.1320                  -0.3125     -0.075440  3 -53.799 113.9 23.20  0.000
# 6  -1.1290 -0.03454         -0.2640                3 -53.824 113.9 23.25  0.000
# 23 -1.1590          -0.3093 -0.3640     -0.112500  4 -52.963 114.4 23.70  0.000
# 8  -1.1540  0.09180 -0.3335 -0.2652                4 -52.980 114.4 23.74  0.000
# 18 -1.1140 -0.04647                      0.078260  3 -54.219 114.7 24.04  0.000
# 20 -1.1360  0.11970 -0.3400             -0.009243  4 -53.404 115.2 24.58  0.000
# 22 -1.1320 -0.01186         -0.3102     -0.070060  4 -53.798 116.0 25.37  0.000
# 24 -1.1660  0.18150 -0.3875 -0.4126     -0.205000  5 -52.768 116.2 25.54  0.000
# Models ranked by AICc(x) 

# retain all predictors
m3 <- glm(WaterNeonicDetection ~ AnnualSnowfall_in + 
                      DaysSinceLastPrecipitation_5mm + AnnualPrecipitation_in +
                      PrecipitationAmount_7days + Event,
                    family = "binomial",
                    na.action = na.fail,
                    data = water.cs)

summary(m3)
confint(m3)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m3, n = 1000) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$AnnualSnowfall_in) #  good
plotResiduals(simulationOutput, form = model.frame(m3)$DaysSinceLastPrecipitation_5mm) # pattern but n.s.
plotResiduals(simulationOutput, form = model.frame(m3)$Event) # good
plotResiduals(simulationOutput, form = model.frame(m3)$PrecipitationAmount_7days) # good
plotResiduals(simulationOutput, form = model.frame(m3)$AnnualPrecipitation_in) # good

#------------------------------------------------------------------------------#

# best drought model
m4 <- glm(WaterNeonicDetection ~ SPEI + Event,
          data = water.cs,
          family = "binomial",
          na.action = na.fail)

car::vif(m4)# vif < 3

summary(m4)
confint(m4)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m4)$Event) # good
plotResiduals(simulationOutput, form = model.frame(m4)$SPEI) # good

#------------------------------------------------------------------------------#

# informed null
m5 <- glm(WaterNeonicDetection ~ Event,
          data = water.cs,
          family = "binomial",
          na.action = na.fail)

summary(m5)
confint(m5)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m5) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m5)$Event) #  good

#------------------------------------------------------------------------------#

# STAGE II: Cross-hypothesis

# AIC model selection
model_names <- paste0("m", 1:5)
models <- mget(model_names)
aictab(models, modnames = model_names)

# not considering porosity
# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m1 5 85.56       0.00   0.77   0.77 -37.45
# m2 6 89.43       3.87   0.11   0.88 -38.25
# m5 4 90.66       5.10   0.06   0.94 -41.11
# m4 5 91.03       5.48   0.05   0.98 -40.19
# m3 8 93.42       7.86   0.02   1.00 -37.89

# when i consider porosity
# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m1 5 85.56       0.00   0.82   0.82 -37.45
# m5 4 90.66       5.10   0.06   0.88 -41.11
# m4 5 91.03       5.48   0.05   0.94 -40.19
# m2 7 91.23       5.68   0.05   0.98 -37.99
# m3 8 93.42       7.86   0.02   1.00 -37.89


#------------------------------------------------------------------------------#

# STAGE III: model average
model_avg <- model.avg(m1, m4, m5)

summary(model_avg)
confint(model_avg)

# Model-averaged coefficients:  
#   (full average) 
#                         Estimate Std. Error Adjusted SE z value Pr(>|z|)  
# (Intercept)              -0.1997     0.9173      0.9262   0.216   0.8293  
# LogCropDistance          -0.5716     0.3646      0.3668   1.558   0.1192  
# EventFall 2023           -3.1809     1.2770      1.2934   2.459   0.0139 *
# EventSpring 2022          0.1767     1.1045      1.1181   0.158   0.8744  
# EventSpring 2023         -0.4519     0.9428      0.9540   0.474   0.6358  
# PermanencePermanent      -0.2559     0.7916      0.7940   0.322   0.7472  
# PermanenceSemipermanent  -0.1652     0.5559      0.5585   0.296   0.7674  

#                             2.5 %      97.5 %
# (Intercept)             -2.014970  1.61550948
# LogCropDistance         -1.237916 -0.15949500
# EventFall 2023          -5.715882 -0.64583465
# EventSpring 2022        -2.014809  2.36821352
# EventSpring 2023        -2.321638  1.41792776
# PermanencePermanent     -4.310612 -0.02372223
# PermanenceSemipermanent -3.274323  0.47704909


#------------------------------------------------------------------------------#

# plot results

# extract estimates
b <- coef(model_avg)

# extract CI
ci <- confint(model_avg)

# combine into a dataframe
df <- data.frame(
  term = names(b),
  estimate = b,
  conf_low = ci[, 1],
  conf_high = ci[, 2]
)

# remove intercept
df <- subset(df, term != "(Intercept)")

# order factors from highest to lowest
  # df$term <- factor(
  #   df$term,
  #   levels = df$term[order(df$estimate)]
  # )

# add reference levels to plot
ref_levels <- data.frame(
  term = c("Sampling Event (Fall 2021)"), 
  estimate = 0,
  conf_low = 0,
  conf_high = 0
)

df <- rbind(df, ref_levels)

# relabel terms
df$term <- factor(
  df$term,
  levels = rev(c("EventSpring 2022",
                 "EventSpring 2023",
                 "EventFall 2023",
                 "Sampling Event (Fall 2021)",
                 "LogCropDistance",
                 "SPEI")),
  labels = rev(c("Sampling Event (Spring 2022)",
                 "Sampling Event (Spring 2023)",
                 "Sampling Event (Fall 2023)",
                 "Sampling Event (Fall 2021)",
                 "Log(Distance to Nearest Crop)",
                 "Drought Index (SPEI)")))

# graph results
ggplot(df, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "grey50",
             linewidth = 0.74) +
  geom_errorbar(aes(xmin = conf_low, 
                     xmax = conf_high), 
                 width = 0.3, 
                 linewidth = 1.15) +
  geom_point(size = 3.5) + my_theme +
  xlim(-10,5) +
  labs(
    x = "Standardized Parameter Estimate",
    y = "")

ggsave(filename = "Model-averaged Water Neonic Results_2026-01-20.tiff",
       path = "C:/Users/shelb/OneDrive - University of Idaho/Masters Project/Analysis/Predictors of Neonics/figures", 
       width = 13, height = 8, units = "in", device='tiff', dpi=300)







ggplot(water, aes(x = Event, y = SPEI)) + geom_boxplot()
