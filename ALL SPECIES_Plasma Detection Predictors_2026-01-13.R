#---------------------------------#
# Predictors of Neonics in Plasma #
#   Analysis by Shelby McCahon    #
#      Created: 2026-01-13        #
#     Modified: 2026-01-13        #
#---------------------------------#

# load packages
library(tidyverse)
library(MuMIn)
library(glmmTMB)
library(DHARMa)

# DETAILS OF ANALYSIS
# STAGE I: Identify the best predictors within each hypothesis
# STAGE II: Identify which hypotheses are supported
# STAGE III: Identify combined effects across hypotheses if multiple hypotheses
#            are supported


# ANALYSIS NOTES: 
# I included sampling event in all models due to VERY strong temporal effect
# Event became informed null
# I considered including site as a random effect but it created model
# convergence issues and singular fit warnings
# All models had VIF < 3

# UNRESOLVED:
# Should distance to nearest wetland be included? If not, model averaging not
# needed because wetland dynamics model not informative. I'm getting opposite 
# result so I'm not sure what to do. Let's wait to see what we learn from

#------------------------------------------------------------------------------#
#                        load data and organize datasets                    ----                        
#------------------------------------------------------------------------------# 

birds <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")

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
birds <- birds %>% 
  mutate(Permanence = case_when(
    Permanence %in% c("Temporary", "Seasonal") ~ "Temporary/Seasonal",
    Permanence == "Semipermanent" ~ "Semipermanent",
    Permanence == "Permanent" ~ "Permanent"))

# convert characters to factors
birds <- birds %>% 
  mutate_if(is.character, as.factor)

# filter data (n = 161)
birds <- birds %>%
  filter(!is.na(PlasmaDetection)) %>% 
  filter(!is.na(Sex)) %>% 
  group_by(Site) %>% 
  filter(n() >= 3) %>%
  group_by(Species) %>% 
  filter(n() >= 3)

# standardize data
birds.cs <- birds %>%
  mutate(across(where(is.numeric), scale))

#------------------------------------------------------------------------------#
#                        random effects structure                           ----                        
#------------------------------------------------------------------------------# 

# should I have site as a random effect? no
# even when species is removed, model without random effect has higher
# explanatory power, although it is close. most models in hypotheses have
# singular fit warnings or model convergence issues, so its best to remove
# site as a random effect

# m1 <- glmmTMB(PlasmaDetection ~ AnnualSnowfall_in + 
#                DaysSinceLastPrecipitation_5mm + PrecipitationAmount_7days +
#                Event + PercentAg + DominantCrop + NearestCropDistance_m + 
#                Permanence + Percent_Exposed_Shoreline + SPEI + 
#                seconds_since_midnight + Sex,
#              family = "binomial",
#              data = birds.cs,
#              REML = TRUE)
# 
# # model not stable when species is included
# m2 <- glmmTMB(PlasmaDetection ~ AnnualSnowfall_in + 
#                 DaysSinceLastPrecipitation_5mm + PrecipitationAmount_7days +
#                 Event + PercentAg + DominantCrop + NearestCropDistance_m + 
#                 Permanence + Percent_Exposed_Shoreline + SPEI + 
#                 seconds_since_midnight + Sex + (1|Site),
#               family = "binomial",
#               data = birds.cs,
#               REML = TRUE)
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)

#------------------------------------------------------------------------------#
#                   identify best model for each hypothesis                 ----                        
#------------------------------------------------------------------------------# 

# include all variables that show up in models with support (delta AICc < 2)

# best agricultural model
global.model <- glm(PlasmaDetection ~ PercentAg + DominantCrop +
                      NearestCropDistance_m + Event,
                    family = "binomial",
                    data = birds.cs,
                    na.action = na.fail)

# check for collinearity
car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#      (Int) DmC Evn   NCD_m      PrA df   logLik  AICc delta weight
# 8   0.9597   +   + -0.6339           8  -69.205 155.4  0.00  0.306
# 7  -0.8143       + -0.8036           5  -72.508 155.4  0.05  0.299
# 15 -0.7893       + -0.8095 -0.03780  6  -72.499 157.5  2.19  0.102
# 16  1.0620   +   + -0.6516 -0.08771  9  -69.181 157.6  2.20  0.102
# 4   0.8038   +   +                   7  -71.417 157.6  2.21  0.101
# 12  0.6516   +   +          0.20660  8  -71.261 159.5  4.11  0.039
# 3  -0.8109       +                   4  -75.724 159.7  4.35  0.035
# 11 -0.8865       +          0.18290  5  -75.411 161.2  5.85  0.016
# 10 -0.1985   +              0.55470  5  -88.889 188.2 32.81  0.000
# 6   0.1982   +     -0.3888           5  -88.989 188.4 33.01  0.000
# 14 -0.1608   +     -0.2539  0.37550  6  -88.283 189.1 33.75  0.000
# 2   0.4700   +                       4  -91.050 190.4 35.00  0.000
# 1  -0.6286                           1 -104.021 210.1 54.71  0.000
# 5  -0.6344         -0.1869           2 -103.496 211.1 55.71  0.000
# 13 -0.6401         -0.3283 -0.25980  3 -102.721 211.6 56.24  0.000
# 9  -0.6295                 -0.07986  2 -103.909 211.9 56.54  0.000
# Models ranked by AICc(x) 

# can remove % ag
m1 <- glm(PlasmaDetection ~ DominantCrop + NearestCropDistance_m + Event,
          family = "binomial",
          data = birds.cs)

summary(m1)
confint(m1)

# model validation --> no severe issues
simulationOutput <- simulateResiduals(fittedModel = m1) 
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m1)$DominantCrop) #  good
plotResiduals(simulationOutput, form = model.frame(m1)$NearestCropDistance_m) # some pattern
plotResiduals(simulationOutput, form = model.frame(m1)$Event)  # good

#------------------------------------------------------------------------------#

# best precipitation model
global.model <- glm(PlasmaDetection ~ DaysSinceLastPrecipitation_5mm +
           PrecipitationAmount_7days + AnnualSnowfall_in + Event,
         family = "binomial",
         data = birds.cs,
         na.action = na.fail)

car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#      (Int)   AnS_in DSL_5mm Evn PrA_7dy df   logLik  AICc delta weight
# 15 -0.8320           0.6834   +  0.5851  6  -73.137 158.8  0.00  0.245
# 16 -0.3443  0.37320  0.6554   +  0.5873  7  -72.433 159.6  0.78  0.166
# 5  -0.8109                    +          4  -75.724 159.7  0.88  0.158
# 6  -0.2546  0.42160           +          5  -74.871 160.1  1.31  0.128
# 7  -0.7717           0.2596   +          5  -75.106 160.6  1.78  0.101
# 13 -0.8456                    +  0.1511  5  -75.461 161.3  2.49  0.071
# 8  -0.2976  0.36170  0.2219   +          6  -74.453 161.5  2.63  0.066
# 14 -0.2367  0.47100           +  0.1912  6  -74.458 161.5  2.64  0.066
# 9  -0.6623                       0.5028  2  -99.882 203.8 45.02  0.000
# 11 -0.6634           0.1234      0.5913  3  -99.765 205.7 46.86  0.000
# 10 -0.6632 -0.06288              0.4899  3  -99.821 205.8 46.98  0.000
# 12 -0.6644 -0.07113  0.1322      0.5832  4  -99.688 207.6 48.81  0.000
# 3  -0.6409          -0.2943              2 -102.565 209.2 50.39  0.000
# 1  -0.6286                               1 -104.021 210.1 51.25  0.000
# 4  -0.6424 -0.10580 -0.2714              3 -102.387 210.9 52.11  0.000
# 2  -0.6322 -0.16080                      2 -103.571 211.2 52.40  0.000
# Models ranked by AICc(x) 

# retain all predictors
m2 <- glm(PlasmaDetection ~ PrecipitationAmount_7days +
            DaysSinceLastPrecipitation_5mm + AnnualSnowfall_in + Event,
          family = "binomial",
          data = birds.cs)

summary(m2)
confint(m2)

# model validation --> no issues
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$PrecipitationAmout_7days) #  good
plotResiduals(simulationOutput, form = model.frame(m2)$AnnualSnowfall_in) # good
plotResiduals(simulationOutput, form = model.frame(m2)$Event)  # good

#------------------------------------------------------------------------------#

# best wetland dynamics model
global.model <- glm(PlasmaDetection ~ Permanence + Percent_Exposed_Shoreline +
                      Event + Dist_Closest_Wetland_m,
                    family = "binomial",
                    data = birds.cs,
                    na.action = na.fail)

car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#      (Int) Dst_Cls_Wtl_m Evn Prc_Exp_Shr Prm df   logLik  AICc delta weight
# 4  -0.7397        0.5827   +                  5  -72.198 154.8  0.00  0.431
# 12 -0.8509        0.8543   +               +  7  -70.532 155.8  1.01  0.260
# 8  -0.7402        0.5883   +    -0.05280      6  -72.167 156.9  2.10  0.151
# 16 -0.8595        0.8409   +    -0.04803   +  8  -70.513 158.0  3.19  0.087
# 3  -0.8109                 +                  4  -75.724 159.7  4.92  0.037
# 11 -0.9326                 +               +  6  -74.516 161.6  6.79  0.014
# 7  -0.8112                 +    -0.01160      5  -75.723 161.8  7.05  0.013
# 15 -0.9638                 +    -0.16560   +  7  -74.286 163.3  8.52  0.006
# 9  -0.5328                                 +  3 -100.820 207.8 53.01  0.000
# 13 -0.5854                      -0.13650   +  4 -100.528 209.3 54.53  0.000
# 10 -0.5276        0.1162                   +  4 -100.615 209.5 54.70  0.000
# 1  -0.6286                                    1 -104.021 210.1 55.28  0.000
# 2  -0.6360        0.2346                      2 -103.057 210.2 55.41  0.000
# 14 -0.5889        0.1455        -0.16170   +  5 -100.218 210.8 56.04  0.000
# 6  -0.6412        0.2701        -0.20460      3 -102.345 210.8 56.06  0.000
# 5  -0.6319                      -0.15970      2 -103.568 211.2 56.43  0.000
# Models ranked by AICc(x) 

# can remove % exposed shoreline
m3 <- glm(PlasmaDetection ~ Permanence + Event + Dist_Closest_Wetland_m,
                    family = "binomial",
                    data = birds.cs)

# model validation --> no severe issues
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$Permanence) #  good
plotResiduals(simulationOutput, form = model.frame(m3)$Event)  # good
plotResiduals(simulationOutput, form = model.frame(m3)$Dist_Closest_Wetland_m)  # some pattern

#------------------------------------------------------------------------------#

# best temporal model
global.model <- glm(PlasmaDetection ~ Event + Julian +
                    seconds_since_midnight,
                    data = birds.cs,
                    family = "binomial",
                    na.action = na.fail)

car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#     (Int) Evn     Jln scn_snc_mdn df   logLik  AICc delta weight
# 2 -0.8109   +                      4  -75.724 159.7  0.00  0.330
# 6 -0.9302   +            -0.31460  5  -74.673 159.7  0.03  0.325
# 8 -0.9254   + -0.3480    -0.35170  6  -74.143 160.8  1.13  0.188
# 4 -0.7919   + -0.2652              5  -75.399 161.2  1.48  0.157
# 3 -0.7336     -0.9720              2  -90.701 185.5 25.77  0.000
# 7 -0.7396     -0.9953    -0.09872  3  -90.576 187.3 27.60  0.000
# 1 -0.6286                          1 -104.021 210.1 50.36  0.000
# 5 -0.6296                 0.08497  2 -103.896 211.9 52.16  0.000
# Models ranked by AICc(x) 

# no evidence to drop any variables

m4 <- glm(PlasmaDetection ~ Event + seconds_since_midnight + Julian,
          family = "binomial",
          data = birds.cs)

# model validation --> no issues
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m4)$seconds_since_midnight) #  good
plotResiduals(simulationOutput, form = model.frame(m4)$Event)  # good
plotResiduals(simulationOutput, form = model.frame(m4)$Julian)  # good

#------------------------------------------------------------------------------#

# drought model
m5 <- glm(PlasmaDetection ~ SPEI + Event,
          data = birds.cs,
          family = "binomial")

car::vif(m5) # vif < 2

summary(m5)
confint(m5)

# model validation --> no issues
simulationOutput <- simulateResiduals(fittedModel = m5) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m5)$SPEI) #  good
plotResiduals(simulationOutput, form = model.frame(m5)$Event)  # good

#------------------------------------------------------------------------------#

# life history model
m6 <- glm(PlasmaDetection ~ Sex + Species + Event,
          data = birds.cs,
          family = "binomial")

car::vif(m6) # vif < 2

summary(m6)
confint(m6) 
# LBDO estimates are unstable because all were captured in spring 2023

# model validation --> no issues
simulationOutput <- simulateResiduals(fittedModel = m6) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m6)$Species) #  good
plotResiduals(simulationOutput, form = model.frame(m6)$Event)  # good
plotResiduals(simulationOutput, form = model.frame(m6)$Sex)  # good

#------------------------------------------------------------------------------#

# best pesticide detection model
m7 <- glm(PlasmaDetection ~ EnvDetection + Event,
          data = birds.cs,
          family = "binomial")

car::vif(m7) # vif < 2

summary(m7)
confint(m7) 

# model validation --> no issues
simulationOutput <- simulateResiduals(fittedModel = m7) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m7)$EnvDetection) #  good
plotResiduals(simulationOutput, form = model.frame(m7)$Event)  # good

#------------------------------------------------------------------------------#

# informed null
m8 <- glm(PlasmaDetection ~ Event,
          data = birds.cs,
          family = "binomial")

#------------------------------------------------------------------------------#

# STAGE II: Cross-hypothesis

# AIC model selection
model_names <- paste0("m", 1:8)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#     K   AICc Delta_AICc AICcWt Cum.Wt     LL
# m1  8 155.36       0.00   0.45   0.45 -69.20 **agriculture
# m3  7 155.80       0.44   0.36   0.82 -70.53 **wetland dynamics
# m2  7 159.60       4.24   0.05   0.87 -72.43
# m8  4 159.70       4.35   0.05   0.92 -75.72
# m5  5 160.79       5.43   0.03   0.95 -75.20
# m4  6 160.83       5.48   0.03   0.98 -74.14
# m7  5 161.83       6.47   0.02   1.00 -75.72
# m6 13 168.05      12.70   0.00   1.00 -69.79


#------------------------------------------------------------------------------#

# STAGE III: model average
model_avg <- model.avg(m1, m3)

summary(model_avg)

# (conditional average) 
# Estimate Std. Error Adjusted SE z value Pr(>|z|)   
# (Intercept)                    0.1537     1.2332      1.2378   0.124  0.90121   
# DominantCropGrassland         -1.0832     0.8759      0.8829   1.227  0.21989   
# DominantCropSoybean           -1.7074     0.9965      1.0045   1.700  0.08918 . 
# DominantCropWheat             -2.1667     0.9212      0.9286   2.333  0.01963 * 
#   NearestCropDistance_m         -0.6339     0.3498      0.3526   1.798  0.07225 . 
# EventFall 2023                -0.9800     0.6205      0.6254   1.567  0.11709   
# EventSpring 2022               4.3036     1.4650      1.4745   2.919  0.00352 **
#   EventSpring 2023              -0.7840     0.8202      0.8265   0.949  0.34282   
# PermanenceSemipermanent        0.5745     0.5042      0.5082   1.130  0.25828   
# PermanenceTemporary/Seasonal  -0.8531     0.9773      0.9850   0.866  0.38644   
# Dist_Closest_Wetland_m         0.8543     0.3132      0.3157   2.706  0.00681 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## relationship of detection ~ distance to wetland doesn't make sense...
## if i remove this variable, the agriculture model is the only one that has
## support, so model averaging wouldn't be needed. let's see if we get this
## result with the other analyses

#------------------------------------------------------------------------------#

# examine relationships
ggplot(birds, aes(y = Dist_Closest_Wetland_m, x = PlasmaDetection)) + 
  geom_boxplot() + my_theme

#------------------------------------------------------------------------------#

# STAGE III: NEW GLOBAL MODEL

global.model <- glm(PlasmaDetection ~ Permanence + Event + 
                      Dist_Closest_Wetland_m + DominantCrop +
                      NearestCropDistance_m,
                    family = "binomial",
                    data = birds.cs,
                    na.action = na.fail)

car::vif(global.model) # vif < 3

dredge(global.model)

# Model selection table 
#         (Int) Dst_Cls_Wtl_m DmC Evn   NCD_m Prm df   logLik  AICc delta weight
# 30 -0.569300       1.29200       + -1.1520   +  8  -65.858 148.7  0.00  0.556
# 14 -0.663500       0.57200       + -0.7108      6  -69.193 150.9  2.27  0.179
# 32  0.270800       1.24100   +   + -1.1020   + 11  -64.315 152.4  3.74  0.086
# 16  0.994000       0.52000   +   + -0.7256      9  -66.843 152.9  4.21  0.068
# 6  -0.739700       0.58270       +              5  -72.198 154.8  6.12  0.026
# 15  0.959700                 +   + -0.6339      8  -69.205 155.4  6.69  0.020
# 13 -0.814300                     + -0.8036      5  -72.508 155.4  6.74  0.019
# 22 -0.850900       0.85430       +           +  7  -70.532 155.8  7.13  0.016
# 8   0.666300       0.43250   +   +              8  -69.727 156.4  7.74  0.012
# 7   0.803800                 +   +              7  -71.417 157.6  8.90  0.006
# 24  0.301400       0.82480   +   +           + 10  -68.635 158.7 10.07  0.004
# 29 -0.864900                     + -0.7290   +  7  -72.239 159.2 10.55  0.003
# 5  -0.810900                     +              4  -75.724 159.7 11.04  0.002
# 31  0.969900                 +   + -0.6277   + 10  -69.200 159.9 11.20  0.002
# 21 -0.932600                     +           +  6  -74.516 161.6 12.91  0.001
# 23  0.912900                 +   +           +  9  -71.368 161.9 13.26  0.001
# 11  0.198200                 +     -0.3888      5  -88.989 188.4 39.70  0.000
# 27  0.521200                 +     -0.3214   +  7  -87.440 189.6 40.95  0.000
# 12  0.006562       0.14980   +     -0.4227      6  -88.716 190.0 41.31  0.000
# 19  0.714700                 +               +  6  -88.802 190.1 41.49  0.000
# 3   0.470000                 +                  4  -91.050 190.4 41.69  0.000
# 28  0.583300      -0.03206   +     -0.3119   +  8  -87.430 191.8 43.14  0.000
# 20  0.970100      -0.14750   +               +  7  -88.572 191.9 43.21  0.000
# 4   0.416700       0.04758   +                  5  -91.021 192.4 43.76  0.000
# 17 -0.532800                                 +  3 -100.820 207.8 59.13  0.000
# 25 -0.463200                       -0.2238   +  4 -100.143 208.5 59.88  0.000
# 18 -0.527600       0.11620                   +  4 -100.615 209.5 60.82  0.000
# 26 -0.446400       0.15010         -0.2437   +  5  -99.810 210.0 61.34  0.000
# 1  -0.628600                                    1 -104.021 210.1 61.40  0.000
# 2  -0.636000       0.23460                      2 -103.057 210.2 61.53  0.000
# 10 -0.641200       0.24970         -0.2018      3 -102.415 211.0 62.32  0.000
# 9  -0.634400                       -0.1869      2 -103.496 211.1 62.40  0.000
# Models ranked by AICc(x) 
