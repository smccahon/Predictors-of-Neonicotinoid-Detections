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
library(AICcmodavg)
library(lattice)
library(dplyr)
library(tidyr)

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

# convert integers to numeric
birds <- birds %>% 
  mutate_if(is.integer, as.numeric)

# filter data (n = 167)
birds <- birds %>%
  filter(!is.na(PlasmaDetection)) %>% 
  filter(!is.na(Sex)) %>% 
  group_by(Species) %>% 
  filter(n() >= 3)

#------------------------------------------------------------------------------#
#                         identify correlations                             ----                        
#------------------------------------------------------------------------------# 

cor_mat <- cor(birds[sapply(birds, is.numeric)], 
               use = "pairwise.complete.obs")

# Convert matrix to a data frame of pairs
cor_df <- as.data.frame(as.table(cor_mat))

# Rename columns
names(cor_df) <- c("Var1", "Var2", "Correlation")

# Convert factors to characters
cor_df$Var1 <- as.character(cor_df$Var1)
cor_df$Var2 <- as.character(cor_df$Var2)

# Keep only unique pairs (upper triangle)
cor_df <- cor_df[cor_df$Var1 < cor_df$Var2, ]

# Filter by threshold
high_corr <- subset(cor_df, abs(Correlation) > 0.6)

# Sort by correlation strength
high_corr <- high_corr[order(-abs(high_corr$Correlation)), ]

high_corr

# relevant correlations

#                               Var1                      Var2 Correlation
# 623 DaysSinceLastPrecipitation_5mm PrecipitationAmount_7days  -0.7599430
# 347                         Julian                      SPEI  -0.6599927
# 233 DaysSinceLastPrecipitation_5mm                    Julian   0.6191898

#------------------------------------------------------------------------------#
#                             standardize data                              ----                        
#------------------------------------------------------------------------------# 

# standardize data
birds.cs <- birds %>%
  mutate(across(where(is.numeric), scale))

# set reference levels
birds.cs$DominantCrop <- relevel(birds.cs$DominantCrop,
                                 ref = "Grassland")

birds.cs$Permanence <- relevel(birds.cs$Permanence,
                               ref = "Temporary/Seasonal")

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
# 7  -0.7958       + -0.7010           5  -75.391 161.2  0.00  0.343
# 8   0.4163   +   + -0.5870           8  -72.489 161.9  0.73  0.237
# 15 -0.7267       + -0.7231 -0.10480  6  -75.319 163.2  2.01  0.126
# 4   0.2897   +   +                   7  -74.557 163.8  2.66  0.090
# 16  0.4298   +   + -0.5896 -0.01346  9  -72.488 164.1  2.97  0.078
# 3  -0.8109       +                   4  -78.152 164.6  3.40  0.063
# 12  0.1144   +   +          0.24890  8  -74.325 165.6  4.41  0.038
# 11 -0.8652       +          0.12790  5  -77.992 166.4  5.20  0.025
# 10 -0.5307   +              0.61620  5  -92.411 195.2 34.04  0.000
# 14 -0.5196   +     -0.2249  0.46670  6  -91.903 196.3 35.18  0.000
# 6  -0.1257   +     -0.3888           5  -93.113 196.6 35.45  0.000
# 2   0.1335   +                       4  -95.283 198.8 37.66  0.000
# 1  -0.6046                           1 -108.460 218.9 57.79  0.000
# 5  -0.6085         -0.1564           2 -108.065 220.2 59.05  0.000
# 13 -0.6148         -0.3021 -0.26880  3 -107.177 220.5 59.35  0.000
# 9  -0.6060                 -0.10400  2 -108.263 220.6 59.44  0.000
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
# days since last precipitation event & precipitation amount are correlated
# which to use?

# this has a slightly higher AIC score (wt = 0.58)
# m1 <- glm(PlasmaDetection ~ DaysSinceLastPrecipitation_5mm +
#             AnnualSnowfall_in + AnnualPrecipitation_in + Event,
#           family = "binomial",
#           data = birds.cs,
#           na.action = na.fail)
# 
# m2 <- glm(PlasmaDetection ~ PrecipitationAmount_7days +
#             AnnualSnowfall_in + AnnualPrecipitation_in + Event,
#           family = "binomial",
#           data = birds.cs,
#           na.action = na.fail)

# AIC model selection
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)

global.model <- glm(PlasmaDetection ~ DaysSinceLastPrecipitation_5mm +
           AnnualSnowfall_in + 
             AnnualPrecipitation_in + Event,
         family = "binomial",
         data = birds.cs,
         na.action = na.fail)

car::vif(global.model) # vif < 2

dredge(global.model)
# 
# Model selection table 
#      (Int) AnP_in   AnS_in DSL_5mm Evn df   logLik  AICc delta weight
# 13 -0.7588                  0.3740   +  5  -76.856 164.1  0.00  0.190
# 10 -0.4871 0.3679                    +  5  -77.013 164.4  0.31  0.163
# 9  -0.8109                           +  4  -78.152 164.6  0.47  0.151
# 15 -0.2550         0.37590  0.3605   +  6  -76.084 164.7  0.61  0.140
# 11 -0.2362         0.42960           +  5  -77.216 164.8  0.72  0.133
# 14 -0.5577 0.2422           0.2878   +  6  -76.421 165.4  1.28  0.100
# 12 -0.2128 0.2715  0.26850           +  6  -76.731 166.0  1.90  0.073
# 16 -0.2459 0.1096  0.31510  0.3220   +  7  -76.017 166.7  2.65  0.050
# 8  -0.7067 1.1180 -0.69200 -0.5428      4  -95.065 198.4 34.29  0.000
# 4  -0.6889 0.9464 -0.70110              3  -98.663 203.5 39.39  0.000
# 6  -0.6759 0.7157          -0.5631      3  -99.727 205.6 41.52  0.000
# 2  -0.6386 0.5169                       2 -103.765 211.6 47.52  0.000
# 5  -0.6170                 -0.2945      2 -106.965 218.0 53.92  0.000
# 1  -0.6046                              1 -108.460 218.9 54.86  0.000
# 7  -0.6180        -0.09226 -0.2781      3 -106.821 219.8 55.70  0.000
# 3  -0.6071        -0.13720              2 -108.118 220.3 56.22  0.000
# Models ranked by AICc(x) 

# retain all predictors
m2 <- glm(PlasmaDetection ~ DaysSinceLastPrecipitation_5mm + AnnualSnowfall_in + 
            AnnualPrecipitation_in + Event,
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
plotResiduals(simulationOutput, form = model.frame(m2)$AnnualPrecipitation_in)  # good

#------------------------------------------------------------------------------#

# best wetland dynamics model
global.model <- glm(PlasmaDetection ~ Permanence +
                      Event + Dist_Closest_Wetland_m,
                    family = "binomial",
                    data = birds.cs,
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

m3 <- glm(PlasmaDetection ~ Permanence + Event + Dist_Closest_Wetland_m,
                    family = "binomial",
                    data = birds.cs)

car::vif(m3)

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

# can drop Julian

m4 <- glm(PlasmaDetection ~ Event + seconds_since_midnight,
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
# m3  7 160.29       0.00   0.56   0.56 -72.79 ** wetland dynamics
# m1  8 161.89       1.60   0.25   0.81 -72.49 ** agriculture
# m8  4 164.55       4.26   0.07   0.87 -78.15 ** informed null
# m4  5 164.96       4.67   0.05   0.93 -77.29 ** temporal
# m5  5 166.42       6.13   0.03   0.95 -78.02
# m7  5 166.67       6.38   0.02   0.98 -78.15
# m2  7 166.74       6.45   0.02   1.00 -76.02
# m6 13 174.29      14.00   0.00   1.00 -72.95

#------------------------------------------------------------------------------#

# STAGE III: model average
model_avg <- model.avg(m1, m4, m3, m8)

summary(model_avg)
confint(model_avg)

# Model-averaged coefficients:  
#   (full average) 
#                         Estimate Std. Error Adjusted SE z value Pr(>|z|)   
# (Intercept)             -1.25028    1.18075     1.18610   1.054  0.29183   
# PermanencePermanent      0.60298    0.88712     0.89181   0.676  0.49895   
# PermanenceSemipermanent  0.89370    1.01821     1.02201   0.874  0.38187   
# EventFall 2023          -0.88832    0.61405     0.61866   1.436  0.15104   
# EventSpring 2022         4.70095    1.48725     1.49637   3.142  0.00168 **
# EventSpring 2023        -0.50820    0.74461     0.74953   0.678  0.49775   
# Dist_Closest_Wetland_m   0.52587    0.49044     0.49132   1.070  0.28448   
# DominantCropCanola       0.14729    0.46870     0.47134   0.312  0.75467   
# DominantCropSoybean     -0.16298    0.44126     0.44340   0.368  0.71319   
# DominantCropWheat       -0.31293    0.62099     0.62249   0.503  0.61516   
# NearestCropDistance_m   -0.15850    0.31235     0.31308   0.506  0.61268   
# seconds_since_midnight  -0.01601    0.08264     0.08289   0.193  0.84679  
  
#                              2.5 %    97.5 %
# (Intercept)             -3.5750031 1.0744401
# PermanencePermanent     -0.8776680 2.8855475
# PermanenceSemipermanent -0.3240878 3.3000422
# EventFall 2023          -2.1008839 0.3242355
# EventSpring 2022         1.7681076 7.6337895
# EventSpring 2023        -1.9772582 0.9608492
# Dist_Closest_Wetland_m   0.2690524 1.4820445
# DominantCropCanola      -0.9797089 2.0707654
# DominantCropSoybean     -1.9360813 0.7287588
# DominantCropWheat       -2.4803879 0.1622846
# NearestCropDistance_m   -1.2413917 0.0673040
# seconds_since_midnight  -0.6991770 0.1476504



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
             "Dist_Closest_Wetland_m",
             "DaysSinceLastPrecipitation_5mm",
             "AnnualSnowfall_in",
             "AnnualPrecipitation_in",
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
              "Nearest Wetland Distance",
             "# Days Since Precipitation Event",
             "Annual Snowfall",
             "Annual Precipitation",
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
  geom_errorbar(aes(xmin = conf_low, 
                     xmax = conf_high), 
                 width = 0.3, 
                 linewidth = 1.15) +
  geom_point(size = 3.5) + my_theme +
  labs(
    x = "Standardized Parameter Estimate",
    y = "")
  





