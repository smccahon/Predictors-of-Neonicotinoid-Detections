#---------------------------------------------------#
# Predictors of Neonics in Lesser Yellowlegs Plasma #
#           Analysis by Shelby McCahon              #
#              Created: 2026-01-15                  #
#             Modified: 2026-01-15                  #
#---------------------------------------------------#

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


# ANALYSIS NOTES: 
# I included sampling event in all models due to VERY strong temporal effect
# Event became informed null
# I considered including site as a random effect but it created model
# convergence issues and singular fit warnings
# All models had VIF < 3

# *** EVENT IS GOING TO CAUSE ISSUES -> ALL 11 BIRDS IN SPRING 2022 HAVE DETECTIONS!!!
# How do i handle event? I can't just do spring/fall because no birds in spring 2023

#------------------------------------------------------------------------------#
#                        load data and organize datasets                    ----                        
#------------------------------------------------------------------------------# 

leye <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")

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
leye <- leye %>% 
  mutate(Permanence = case_when(
    Permanence %in% c("Temporary", "Seasonal") ~ "Temporary/Seasonal",
    Permanence == "Semipermanent" ~ "Semipermanent",
    Permanence == "Permanent" ~ "Permanent"))

# filter to Lesser Yellowlegs (n = 54)
leye <- subset(leye, Species == "Lesser Yellowlegs")

# convert characters to factors
leye <- leye %>% 
  mutate_if(is.character, as.factor)

# standardize data
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# set reference levels
leye.cs$DominantCrop <- relevel(leye.cs$DominantCrop,
                                 ref = "Grassland")

leye.cs$Permanence <- relevel(leye.cs$Permanence,
                               ref = "Temporary/Seasonal")

#------------------------------------------------------------------------------#
#                   identify best model for each hypothesis                 ----                        
#------------------------------------------------------------------------------# 

# include all variables that show up in models with support (delta AICc < 2)

# best agricultural model
global.model <- glm(PlasmaDetection ~ PercentAg + DominantCrop +
                      NearestCropDistance_m + Event,
                    family = "binomial",
                    data = leye.cs,
                    na.action = na.fail)

# check for collinearity
car::vif(global.model) # vif < 3

dredge(global.model)

# Model selection table 
#         (Int) DmC Evn   NCD_m      PrA df  logLik AICc delta weight
# 7   1.322e-01       + -0.6554           4 -22.575 54.0  0.00  0.350
# 3  -1.542e-01       +                   3 -23.984 54.4  0.48  0.275
# 11 -2.657e-01       +          0.39690  4 -23.303 55.4  1.45  0.169
# 15  3.243e-01       + -0.8707 -0.24960  5 -22.498 56.2  2.28  0.112
# 4  -7.984e-01   +   +                   5 -23.504 58.3  4.29  0.041
# 8   2.702e-01   +   + -0.6952           6 -22.558 58.9  4.94  0.030
# 12 -3.010e-01   +   +          0.38120  6 -23.302 60.4  6.42  0.014
# 16  1.933e-01   +   + -0.9107 -0.37190  7 -22.459 61.4  7.39  0.009
# 5  -3.532e-01         -0.6000           2 -35.623 75.5 21.51  0.000
# 1  -2.985e-01                           1 -36.835 75.7 21.78  0.000
# 13 -3.821e-01         -0.9862 -0.42780  3 -35.001 76.5 22.52  0.000
# 6   2.402e-01   +     -0.8633           4 -34.202 77.2 23.25  0.000
# 9  -2.987e-01                  0.04690  2 -36.821 77.9 23.91  0.000
# 2   6.483e-16   +                       3 -36.272 79.0 25.06  0.000
# 14  2.086e-01   +     -0.8881 -0.04065  5 -34.200 79.6 25.68  0.000
# 10  4.743e-01   +              0.49230  4 -35.599 80.0 26.05  0.000
# Models ranked by AICc(x) 

# can remove dominant crop
m1 <- glm(PlasmaDetection ~ PercentAg + NearestCropDistance_m + Event,
          family = "binomial",
          data = leye.cs)

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
                    data = leye.cs,
                    na.action = na.fail)

car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#      (Int)   AnS_in DSL_5mm Evn PrA_7dy df   logLik  AICc delta weight
# 15 -0.8221           0.8863   + 0.67500  6  -74.536 161.6  0.00  0.344
# 16 -0.2883  0.40050  0.8952   + 0.69120  7  -73.655 162.0  0.42  0.280
# 7  -0.7588           0.3740   +          5  -76.856 164.1  2.49  0.099
# 5  -0.8109                    +          4  -78.152 164.6  2.95  0.079
# 8  -0.2550  0.37590  0.3605   +          6  -76.084 164.7  3.09  0.073
# 6  -0.2362  0.42960           +          5  -77.216 164.8  3.21  0.069
# 13 -0.8273                    + 0.07713  5  -78.085 166.5  4.94  0.029
# 14 -0.2342  0.44880           + 0.10490  6  -77.093 166.7  5.11  0.027
# 9  -0.6340                      0.47670  2 -104.523 213.1 51.52  0.000
# 11 -0.6348           0.1687     0.60320  3 -104.330 214.8 53.21  0.000
# 10 -0.6345 -0.05537             0.46670  3 -104.473 215.1 53.50  0.000
# 12 -0.6353 -0.06031  0.1726     0.59530  4 -104.271 216.8 55.19  0.000
# 3  -0.6170          -0.2945              2 -106.965 218.0 56.41  0.000
# 1  -0.6046                               1 -108.460 218.9 57.35  0.000
# 4  -0.6180 -0.09226 -0.2781              3 -106.821 219.8 58.19  0.000
# 2  -0.6071 -0.13720                      2 -108.118 220.3 58.71  0.000
# Models ranked by AICc(x) 

# retain all predictors
m2 <- glm(PlasmaDetection ~ PrecipitationAmount_7days +
            DaysSinceLastPrecipitation_5mm + AnnualSnowfall_in + Event,
          family = "binomial",
          data = leye.cs)

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
                    data = leye.cs,
                    na.action = na.fail)

car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#      (Int) Dst_Cls_Wtl_m Evn Prc_Exp_Shr Prm df   logLik  AICc delta weight
# 4  -0.7393       0.58570   +                  5  -74.436 159.2  0.00  0.441
# 12 -0.8345       0.87550   +               +  7  -72.793 160.3  1.04  0.261
# 8  -0.7389       0.59300   +   -0.056660      6  -74.401 161.3  2.08  0.155
# 16 -0.8407       0.86630   +   -0.035730   +  8  -72.783 162.5  3.23  0.088
# 3  -0.8109                 +                  4  -78.152 164.6  5.31  0.031
# 7  -0.8109                 +   -0.001477      5  -78.152 166.7  7.43  0.011
# 11 -0.9119                 +               +  6  -77.201 166.9  7.68  0.009
# 15 -0.9411                 +   -0.162000   +  7  -76.982 168.7  9.42  0.004
# 9  -0.5596                                 +  3 -104.443 215.0 55.79  0.000
# 10 -0.5547       0.09329                   +  4 -104.303 216.9 57.61  0.000
# 13 -0.5989                     -0.094060   +  4 -104.305 216.9 57.61  0.000
# 14 -0.6025       0.11580       -0.116600   +  5 -104.099 218.6 59.33  0.000
# 1  -0.6046                                    1 -108.460 218.9 59.70  0.000
# 2  -0.6109       0.21980                      2 -107.575 219.2 59.98  0.000
# 6  -0.6141       0.24870       -0.152700      3 -107.165 220.5 61.23  0.000
# 5  -0.6061                     -0.106900      2 -108.250 220.6 61.33  0.000
# Models ranked by AICc(x) 

# can remove % exposed shoreline
m3 <- glm(PlasmaDetection ~ Permanence + Event + Dist_Closest_Wetland_m,
          family = "binomial",
          data = leye.cs)

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
                    data = leye.cs,
                    family = "binomial",
                    na.action = na.fail)

car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#     (Int) Evn      Jln scn_snc_mdn df   logLik  AICc delta weight
# 2 -0.8109   +                       4  -78.152 164.6  0.00  0.404
# 6 -0.9191   +             -0.27580  5  -77.295 165.0  0.41  0.329
# 4 -0.8039   + -0.08583              5  -78.117 166.6  2.06  0.144
# 8 -0.9134   + -0.13980    -0.28720  6  -77.207 166.9  2.39  0.122
# 3 -0.6982     -0.91650              2  -95.647 195.4 30.82  0.000
# 7 -0.7013     -0.92840    -0.05596  3  -95.604 197.4 32.80  0.000
# 1 -0.6046                           1 -108.460 218.9 54.39  0.000
# 5 -0.6059                  0.09691  2 -108.291 220.7 56.10  0.000
# Models ranked by AICc(x) 

# no evidence to drop any variables

m4 <- glm(PlasmaDetection ~ Event + seconds_since_midnight + Julian,
          family = "binomial",
          data = leye.cs)

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
          data = leye.cs,
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
          data = leye.cs,
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
          data = leye.cs,
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
          data = leye.cs,
          family = "binomial")

#------------------------------------------------------------------------------#

# STAGE II: Cross-hypothesis

# AIC model selection
model_names <- paste0("m", 1:8)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#   K   AICc Delta_AICc AICcWt Cum.Wt     LL
# m3  7 160.29       0.00   0.47   0.47 -72.79 ** wetland dynamics
# m1  8 161.89       1.60   0.21   0.69 -72.49 ** agriculture
# m2  7 162.01       1.72   0.20   0.88 -73.66 ** precipitation
# m8  4 164.55       4.26   0.06   0.94 -78.15
# m5  5 166.42       6.13   0.02   0.96 -78.02
# m7  5 166.67       6.38   0.02   0.98 -78.15
# m4  6 166.94       6.65   0.02   1.00 -77.21
# m6 13 174.29      14.00   0.00   1.00 -72.95

#------------------------------------------------------------------------------#

# STAGE III: model average
model_avg <- model.avg(m1, m2, m3, m8)

summary(model_avg)
confint(model_avg)

# Model-averaged coefficients:  
#   (full average) 
#                                Estimate Std. Error Adjusted SE z value Pr(>|z|)   
# (Intercept)                    -1.06234    1.17965     1.18458   0.897  0.36982   
# PermanencePermanent             0.50440    0.84145     0.84559   0.597  0.55084   
# PermanenceSemipermanent         0.74758    0.98817     0.99144   0.754  0.45083   
# EventFall 2023                 -1.07436    0.74587     0.75033   1.432  0.15219   
# EventSpring 2022                4.59082    1.46372     1.47269   3.117  0.00183 **
# EventSpring 2023               -0.53188    0.75624     0.76166   0.698  0.48497   
# Dist_Closest_Wetland_m          0.43989    0.48890     0.48964   0.898  0.36898   
# DominantCropCanola              0.12321    0.43212     0.43452   0.284  0.77676   
# DominantCropSoybean            -0.13634    0.40805     0.40999   0.333  0.73948   
# DominantCropWheat              -0.26177    0.57963     0.58097   0.451  0.65230   
# NearestCropDistance_m          -0.13258    0.29162     0.29228   0.454  0.65010   
# PrecipitationAmount_7days       0.14660    0.31921     0.31974   0.458  0.64660   
# DaysSinceLastPrecipitation_5mm  0.18985    0.40022     0.40072   0.474  0.63566   
# AnnualSnowfall_in               0.08493    0.21894     0.21967   0.387  0.69904  

#                                     2.5 %    97.5 %
# (Intercept)                    -3.3840801 1.2593923
# PermanencePermanent            -0.8776680 2.8855475
# PermanenceSemipermanent        -0.3240878 3.3000422
# EventFall 2023                 -2.5449887 0.3962606
# EventSpring 2022                1.7044028 7.4772410
# EventSpring 2023               -2.0247041 0.9609348
# Dist_Closest_Wetland_m          0.2690524 1.4820445
# DominantCropCanola             -0.9797089 2.0707654
# DominantCropSoybean            -1.9360813 0.7287588
# DominantCropWheat              -2.4803879 0.1622846
# NearestCropDistance_m          -1.2413917 0.0673040
# PrecipitationAmount_7days       0.0544267 1.3280695
# DaysSinceLastPrecipitation_5mm  0.2002021 1.5902052
# AnnualSnowfall_in              -0.2229812 1.0239035



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
  term = c("Sampling Event (Fall 2021)", 
           "Wetland Permanence (Temporary/Seasonal)",
           "Dominant Land Use (Grassland)"), 
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
                 "PermanenceSemipermanent",
                 "PermanencePermanent",
                 "Wetland Permanence (Temporary/Seasonal)",
                 "DaysSinceLastPrecipitation_5mm",
                 "Dist_Closest_Wetland_m",
                 "PrecipitationAmount_7days",
                 "AnnualSnowfall_in",
                 "DominantCropCanola",
                 "DominantCropSoybean",
                 "DominantCropWheat",
                 "Dominant Land Use (Grassland)",
                 "NearestCropDistance_m")),
  labels = rev(c("Sampling Event (Spring 2022)",
                 "Sampling Event (Spring 2023)",
                 "Sampling Event (Fall 2023)",
                 "Sampling Event (Fall 2021)",
                 "Wetland Permanence (Semi-permanent)",
                 "Wetland Permanence (Permanent)",
                 "Wetland Permanence (Temporary/Seasonal)",
                 "# Days Since Precipitation Event",
                 "Nearest Wetland Distance",
                 "Precipitation Amount (Last 7 days)",
                 "Annual Snowfall",
                 "Dominant Land Use (Canola)",
                 "Dominant Land Use (Soybean)",
                 "Dominant Land Use (Wheat)",
                 "Dominant Land Use (Grassland)",
                 "Nearest Crop Distance")))

# graph results
ggplot(df, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "grey50",
             linewidth = 0.74) +
  geom_errorbarh(aes(xmin = conf_low, 
                     xmax = conf_high), 
                 width = 0.3, 
                 linewidth = 1.15) +
  geom_point(size = 3.5) + my_theme +
  labs(
    x = "Standardized Parameter Estimate",
    y = "")








# standardize data
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))