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
library(lme4)
library(tidyr)

# DETAILS OF ANALYSIS
# STAGE I: Identify the best predictors within each hypothesis (models < 2 delta AICc)
# STAGE II: Identify which hypotheses are supported (Identify 95% Confidence Set)
# STAGE III: Model average


# ANALYSIS NOTES: 
# site was initially considered as a nested random effect within event to
# account for multiple captures at a site (max = 25 individuals), but
# the estimated variance was negligible and resulted in a singular fit warning. 
# Therefore, only event was
# included as a random effect in the final models. 
# Season was a more informative predictor than Julian day

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

# remove species with no variation in plasma detection (n = 161)
birds <- birds %>% 
  filter(!Species == "Longbilled Dowitcher")

#------------------------------------------------------------------------------#
#                         identify correlations                             ----                        
#------------------------------------------------------------------------------# 
# 
# cor_mat <- cor(birds[sapply(birds, is.numeric)], 
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

#                               Var1                      Var2 Correlation
# 347                         Julian                      SPEI  -0.6599927

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

# should I have event as a random effect? YES (aic wt = 0.99)

#  m1 <- glmmTMB(PlasmaDetection ~ PrecipitationAmount_7days + 
#                  AnnualPrecipitation_in + AnnualSnowfall_in +
#                  DaysSinceLastPrecipitation_5mm +
#                  PercentAg + DominantCrop +  
#                  Permanence + SPEI + Julian +
#                  seconds_since_midnight + Sex + (1|Event),
#                family = "binomial",
#                data = birds.cs,
#                REML = TRUE)
#  
#  m2 <- glmmTMB(PlasmaDetection ~ PrecipitationAmount_7days + 
#                  AnnualPrecipitation_in + AnnualSnowfall_in +
#                  DaysSinceLastPrecipitation_5mm +
#                  PercentAg + DominantCrop  + 
#                  Permanence + SPEI + Julian +
#                  seconds_since_midnight + Sex,
#                family = "binomial",
#                data = birds.cs,
#                REML = TRUE)
#  
#  model_names <- paste0("m", 1:2)
#  models <- mget(model_names)
#  aictab(models, modnames = model_names)

#------------------------------------------------------------------------------#
#                   identify best model for each hypothesis                 ----                        
#------------------------------------------------------------------------------# 

# include all variables that show up in models with support (delta AICc < 2)

# best agricultural model
global.model <- glmer(PlasmaDetection ~ PercentAg + DominantCrop +
                      (1|Event),
                    family = "binomial",
                    data = birds.cs,
                    na.action = na.fail)

# check for collinearity
car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#    (Intrc) DmnnC   PrcnA df  logLik  AICc delta weight
# 2 0.35400     +          5 -79.303 169.0  0.00  0.415
# 1 0.00688                2 -82.929 169.9  0.94  0.259
# 4 0.73120     + 0.34590  6 -78.798 170.1  1.15  0.234
# 3 0.00523       0.01998  3 -82.925 172.0  3.01  0.092
# Models ranked by AICc(x) 
# Random terms (all models): 
#   1 | Event

# retain both predictors
m1 <- glmer(PlasmaDetection ~ PercentAg + DominantCrop + (1|Event),
          family = "binomial",
          data = birds.cs)

summary(m1)
confint(m1, method = "Wald") 

# model validation --> patterns not n.s.; likely to uneven sample sizes in
# dominant crop variable
simulationOutput <- simulateResiduals(fittedModel = m1, n = 1000) 
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m1)$DominantCrop) #  significant, low n
plotResiduals(simulationOutput, form = model.frame(m1)$PercentAg) # good

#------------------------------------------------------------------------------#

# best wetland dynamics model
global.model <- glmer(PlasmaDetection ~ Permanence + Dist_Closest_Wetland_m +
                      (1|Event),
                    family = "binomial",
                    data = birds.cs,
                    na.action = na.fail)

car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#      (Int) Dst_Cls_Wtl_m Prm df  logLik  AICc delta weight
# 2 -0.08022        0.5602      3 -82.424 171.0  0.00  0.679
# 4 -0.66130        0.7573   +  5 -81.352 173.1  2.08  0.240
# 1 -0.07147                    2 -85.866 175.8  4.81  0.061
# 3  0.47830                 +  4 -84.929 178.1  7.11  0.019
# Models ranked by AICc(x) 
# Random terms (all models): 
#   1 | Event

# remove wetland permanence
m2 <- glmer(PlasmaDetection ~ (1|Event) + Dist_Closest_Wetland_m,
                    family = "binomial",
                    data = birds.cs)

# model validation --> pattern but not severe
simulationOutput <- simulateResiduals(fittedModel = m2, n = 1000) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput) # pattern but almost n.s.
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$Dist_Closest_Wetland_m)  # some pattern but n.s.

#------------------------------------------------------------------------------#

# best temporal model
# m1 <- glmer(PlasmaDetection ~ (1|Event) + Season + 
#                     seconds_since_midnight,
#                     data = birds.cs,
#                     family = "binomial",
#                     na.action = na.fail)
# 
# m2 <- glmer(PlasmaDetection ~ (1|Event) + Julian + 
#                         seconds_since_midnight,
#                       data = birds.cs,
#                       family = "binomial",
#                       na.action = na.fail)
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names) # season is a better predictor

global.model <- glmer(PlasmaDetection ~ (1|Event) + Season + 
              seconds_since_midnight,
            data = birds.cs,
            family = "binomial",
            na.action = na.fail)

car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#      (Int) Ssn scn_snc_mdn df  logLik  AICc delta weight
# 2 -1.26000   +              3 -84.689 175.5  0.00  0.312
# 1 -0.07147                  2 -85.866 175.8  0.28  0.271
# 4 -1.32400   +     -0.2512  4 -83.971 176.2  0.66  0.224
# 3 -0.06776         -0.2467  3 -85.172 176.5  0.96  0.193
# Models ranked by AICc(x) 
# Random terms (all models): 
#   1 | Event

# retain all predictors
m3 <- glmer(PlasmaDetection ~ (1|Event) + Season + seconds_since_midnight,
          family = "binomial",
          data = birds.cs)

# model validation --> no issues
simulationOutput <- simulateResiduals(fittedModel = m3, n = 1000) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$seconds_since_midnight) #  good
plotResiduals(simulationOutput, form = model.frame(m3)$Julian)  # good

#------------------------------------------------------------------------------#

# drought model
m4 <- glmer(PlasmaDetection ~ SPEI + (1|Event),
          data = birds.cs,
          family = "binomial")

summary(m4)
confint(m4)

# model validation --> issues
simulationOutput <- simulateResiduals(fittedModel = m4, n = 1000) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m5)$SPEI) #  good
plotResiduals(simulationOutput, form = model.frame(m5)$Event)  # good

#------------------------------------------------------------------------------#

# life history model
# issues with long-billed dowitchers ()
m5 <- glmer(PlasmaDetection ~ Sex + Species + (1|Event),
          data = birds.cs,
          family = "binomial")

car::vif(m5) # vif < 2

summary(m5)
confint(m5)

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
m6 <- glmer(PlasmaDetection ~ EnvDetection + (1|Event),
          data = birds.cs,
          family = "binomial")

car::vif(m7) # vif < 2

summary(m7)
confint(m7) 

# model validation --> no issues
simulationOutput <- simulateResiduals(fittedModel = m7, n = 1000) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m7)$EnvDetection) #  good
plotResiduals(simulationOutput, form = model.frame(m7)$Event)  # good

#------------------------------------------------------------------------------#

# informed null
m7 <- glmer(PlasmaDetection ~ (1|Event),
          data = birds.cs,
          family = "binomial")

#------------------------------------------------------------------------------#

# STAGE II: Cross-hypothesis

# AIC model selection
model_names <- paste0("m", 1:7)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#     K   AICc Delta_AICc AICcWt Cum.Wt     LL
# m2  3 171.00       0.00   0.60   0.60 -82.42
# m1  6 172.52       1.52   0.28   0.88 -80.00
# m7  2 175.80       4.81   0.05   0.93 -85.87
# m4  3 177.10       6.10   0.03   0.96 -85.48
# m6  3 177.88       6.88   0.02   0.98 -85.87
# m3  4 177.95       6.96   0.02   1.00 -84.85
# m5 11 184.52      13.53   0.00   1.00 -80.41

#------------------------------------------------------------------------------#

# STAGE III: model average
model_avg <- model.avg(m1, m2, m7, m3)

summary(model_avg)
confint(model_avg)

# (full average) 
#                         Estimate Std. Error Adjusted SE z value Pr(>|z|)
# (Intercept)            -0.007747   1.034001    1.040837   0.007    0.994
# Dist_Closest_Wetland_m  0.343111   0.320346    0.321001   1.069    0.285
# DominantCropCanola      0.122681   0.457418    0.460263   0.267    0.790
# DominantCropSoybean    -0.207676   0.472591    0.474452   0.438    0.662
# DominantCropWheat      -0.382575   0.699536    0.700891   0.546    0.585
# NearestCropDistance_m  -0.152125   0.285479    0.286113   0.532    0.595
# SeasonSpring            0.112623   0.595809    0.596942   0.189    0.850
# seconds_since_midnight -0.011472   0.069219    0.069440   0.165    0.869
  
#                             2.5 %      97.5 %
# (Intercept)            -2.0477499  2.03225628
# Dist_Closest_Wetland_m  0.1368806  0.98344817 *
# DominantCropCanola     -1.1008091  1.95729991 
# DominantCropSoybean    -1.9811394  0.53126221
# DominantCropWheat      -2.6389017 -0.03201642 *
# NearestCropDistance_m  -1.1009259  0.03887639
# SeasonSpring           -0.3051902  5.23711145
# seconds_since_midnight -0.6686544  0.16629401

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
  term = c("Season (Fall)", 
           "Dominant Land Use (Grassland)"), 
  estimate = 0,
  conf_low = 0,
  conf_high = 0
)

df <- rbind(df, ref_levels)

# relabel terms
df$term <- factor(
  df$term,
  levels = rev(c("SeasonSpring",
                 "Season (Fall)",
                 "seconds_since_midnight",
                 "Dist_Closest_Wetland_m",
                 "DominantCropCanola",
                 "DominantCropSoybean",
                 "DominantCropWheat",
                 "Dominant Land Use (Grassland)",
                 "NearestCropDistance_m")),
  labels = rev(c("Season (Spring)",
                 "Season (Fall)",
                 "Capture Time",
                 "Distance to Nearest Wetland",
                 "Dominant Land Use (Canola)",
                 "Dominant Land Use (Soybean)",
                 "Dominant Land Use (Wheat)",
                 "Dominant Land Use (Grassland)",
                 "Distance to Nearest Crop")))

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
  




