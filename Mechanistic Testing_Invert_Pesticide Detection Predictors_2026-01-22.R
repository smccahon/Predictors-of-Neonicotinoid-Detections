#-----------------------------------#
#  Predictors of Neonics in Inverts #
#     Analysis by Shelby McCahon    #
#        Created: 2026-01-22        #
#       Modified: 2026-01-22        #
#-----------------------------------#

# load packages
library(tidyverse)
library(MuMIn)
library(glmmTMB)
library(DHARMa)
library(AICcmodavg)

# DETAILS OF ANALYSIS
# STAGE I: Identify the best candidate model for each hypothesis using AICc
# STAGE II: Identify which hypotheses are supported using model selection
# STAGE III: Model-average top predictors (95% confidence set)

# ANALYSIS NOTES



#------------------------------------------------------------------------------#
#                        load data and organize datasets                    ----                        
#------------------------------------------------------------------------------# 

invert <- read.csv("cleaned_data/invert_data_cleaned_2025-08-11.csv")

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

# convert for correlation matrix
 # invert <- invert %>%
 #   mutate(Buffered = ifelse(Buffered == "Y", 1, 0))
 # 
 # invert <- invert %>%
 #   mutate(Season = ifelse(Season == "Spring", 1, 0))

# convert characters to factors
invert <- invert %>% 
  mutate_if(is.character, as.factor)

# filter data (n = 51)
invert <- invert %>%
  filter(!is.na(InvertPesticideDetection))

# log transform nearest crop distance due to skew
invert <- invert %>% 
  mutate(LogCropDistance = NearestCropDistance_m + 0.0001) %>% 
  mutate(LogCropDistance = log(LogCropDistance))

# log transform nearest wetland distance due to skew
invert <- invert %>% 
  mutate(LogWetlandDistance = Dist_Closest_Wetland_m) %>% 
  mutate(LogWetlandDistance = log(LogWetlandDistance))

# log transform precipitation amount due to skew
invert <- invert %>% 
  mutate(LogPrecipAmount = PrecipitationAmount_7days) %>% 
  mutate(LogPrecipAmount = log(LogPrecipAmount))

# log transform precipitation events due to skew
invert <- invert %>% 
  mutate(LogPrecipDays = DaysSinceLastPrecipitation_5mm + 0.0001) %>% 
  mutate(LogPrecipDays = log(LogPrecipDays))

# standardize data
invert.cs <- invert %>%
  mutate(across(where(is.numeric), scale))

# set reference levels
invert.cs$DominantCrop <- relevel(invert.cs$DominantCrop,
                                 ref = "Grassland")

invert.cs$Permanence <- relevel(invert.cs$Permanence,
                               ref = "Temporary")

#------------------------------------------------------------------------------#
#                         identify correlations                             ----                        
#------------------------------------------------------------------------------# 

# code from a previous analysis (only relevant ones are shown)

#                               Var1                      Var2 Correlation
#    8                         Julian                    Season  -0.9387924
# 156                       Buffered                 PercentAg  -0.7439067
# 151          NearestCropDistance_m                 PercentAg  -0.7261301
# 72                        Buffered     NearestCropDistance_m   0.6990801
# 288 DaysSinceLastPrecipitation_5mm PrecipitationAmount_7days  -0.6370290

#------------------------------------------------------------------------------#
#                   identify best model for each hypothesis                 ----                        
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# ***think more about other hypotheses...check to see if correlations cause
#  issues
#------------------------------------------------------------------------------#

# log transformation for crop distance improved fit (wt = 0.55)

# *agriculture candidate models ----

###... LANDSCAPE-SCALE EXPOSURE HYPOTHESIS ----
# (detections are influenced by overall agricultural intensity)
# Mechanism: inverts in wetlands with more cropland cover are more likely to
# encounter neonics
m1 <- glm(InvertPesticideDetection ~ PercentAg, 
          data = invert.cs,
          family = "binomial")

###... LOCAL CROP PROXIMITY HYPOTHESIS ----
# (detections are influenced by local proximity to crops)
# Mechanism: Inverts in wetlands close to crop fields are more likely to be exposed
# log transformation improves fit (wt. = 0.95)

m2 <- glm(InvertPesticideDetection ~ LogCropDistance, 
          data = invert.cs,
          family = "binomial")

###... CROP-TYPE SPECIFICITY HYPOTHESIS ----
# (detections are influenced by different crop-type exposure)
# Mechanism: Certain crop types are more heavily treated with neonics which
# affects exposure and detection
m3 <- glm(InvertPesticideDetection ~ DominantCrop, 
          data = invert.cs,
          family = "binomial")

###... VEGETATION BUFFER HYPOTHESIS ----
# (detections are reduced in wetlands with shoreline vegetation buffers)
# Mechanism: Exposure is lower in wetlands with buffers (50m)

m4 <- glm(InvertPesticideDetection ~ Buffered, 
          data = invert.cs,
          family = "binomial")


###... LOCAL AND LANDSCAPE EXPOSURE HYPOTHESIS ----
# (detections are influenced by local and landscape-scale cropland cover)
# Mechanism: Exposure depends on both proximity and landscape-level
# cropland cover intensity

### these variables are moderately correlated! 
m5 <- glm(InvertPesticideDetection ~ PercentAg + LogCropDistance, 
          data = invert.cs,
          family = "binomial")




model_names <- paste0("m", 1:4)
models <- mget(model_names)
aictab(models, modnames = model_names)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m1, n = 1000) 
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m1)$CropDistance) # good

#------------------------------------------------------------------------------#

# best wetland dynamics model
# log-transformation? YES (although similar wt. = 0.55)
#  m1 <- glm(InvertPesticideDetection ~ LogWetlandDistance,
#            family = "binomial",
#            data = invert.cs,
#            na.action = na.fail)
# #
#  m2 <- glm(InvertPesticideDetection ~ Dist_Closest_Wetland_m,
#            family = "binomial",
#            data = invert.cs,
#            na.action = na.fail)
#  
#  model_names <- paste0("m", 1:2)
#  models <- mget(model_names)
#  aictab(models, modnames = model_names)


global.model <- glm(InvertPesticideDetection ~ Permanence + 
                      LogWetlandDistance,
                    family = "binomial",
                    na.action = na.fail,
                    data = invert.cs)

# check for collinearity
car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#   (Intrc)  LgWtD Prmnn df  logLik AICc delta weight
# 3 -1.6090            +  4 -31.563 72.0  0.00  0.357
# 1 -0.1967               1 -35.105 72.3  0.30  0.308
# 4 -1.5330 0.3447     +  5 -30.940 73.2  1.22  0.194
# 2 -0.1992 0.2241        2 -34.799 73.8  1.85  0.141
# Models ranked by AICc(x) 

# keep both
m2 <- glm(InvertPesticideDetection ~ Permanence + LogWetlandDistance,
          data = invert.cs,
          family = "binomial",
          na.action = na.fail)

summary(m2)
confint(m2)

# model validation --> some issues in the residuals but log transformation helped
simulationOutput <- simulateResiduals(fittedModel = m2, n = 1000) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$Permanence) #  good
plotResiduals(simulationOutput, form = model.frame(m2)$LogWetlandDistance) #  some pattern

#------------------------------------------------------------------------------#
# best precipitation model

# log-transformation? no 
# m1 <- glm(InvertPesticideDetection ~ LogPrecipAmount,
#           family = "binomial",
#           data = invert.cs,
#           na.action = na.fail)
# 
# m2 <- glm(InvertPesticideDetection ~ PrecipitationAmount_7days,
#           family = "binomial",
#           data = invert.cs,
#           na.action = na.fail)
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)
  
# log-transformation? YES
# m1 <- glm(InvertPesticideDetection ~ LogPrecipDays,
#           family = "binomial",
#           data = invert.cs,
#           na.action = na.fail)
# 
# m2 <- glm(InvertPesticideDetection ~ DaysSinceLastPrecipitation_5mm,
#           family = "binomial",
#           data = invert.cs,
#           na.action = na.fail)
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)

#### taking into consideration of correlations

# precipitation amount performed better than days since event (AICc wt = 0.80)
# m1 <- glm(InvertPesticideDetection ~ LogPrecipDays,
#           family = "binomial",
#           data = invert.cs,
#           na.action = na.fail) 
# 
# m2 <- glm(InvertPesticideDetection ~ PrecipitationAmount_7days,
#           family = "binomial",
#           data = invert.cs,
#           na.action = na.fail)
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)
 

# including all correlated variables
 global.model <- glm(InvertPesticideDetection ~ AnnualSnowfall_in + 
                       AnnualPrecipitation_in + PrecipitationAmount_7days +
                       LogPrecipDays,
                     family = "binomial",
                     na.action = na.fail,
                     data = invert.cs)
 
# removing correlated variables
# retain all variables...
 global.model <- glm(InvertPesticideDetection ~ AnnualSnowfall_in + 
                       AnnualPrecipitation_in + PrecipitationAmount_7days,
                     family = "binomial",
                     na.action = na.fail,
                     data = invert.cs)
 
# check for collinearity
car::vif(global.model) # vif < 2

dredge(global.model)

# considering correlations
# retain all variables...
m3 <- glm(InvertPesticideDetection ~ PrecipitationAmount_7days +
            AnnualSnowfall_in + 
            AnnualPrecipitation_in + LogPrecipDays,
          family = "binomial",
          na.action = na.fail,
          data = invert.cs)


# model validation --> no severe issues
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$AnnualSnowfall_in) # good
plotResiduals(simulationOutput, form = model.frame(m3)$PrecipitationAmount_7days) # pattern
plotResiduals(simulationOutput, form = model.frame(m3)$AnnualPrecipitation_in) # good

#------------------------------------------------------------------------------#

# temporal model
m4 <- glm(InvertPesticideDetection ~ Season,
          data = invert.cs,
          family = "binomial",
          na.action = na.fail)

summary(m4)
confint(m4)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m4)$Season) #  good

#------------------------------------------------------------------------------#

# null model

m5 <- glm(InvertPesticideDetection ~ 1,
          family = "binomial",
          data = invert.cs)

# STAGE II: Cross-hypothesis

# AIC model selection
model_names <- paste0("m", 1:5)
models <- mget(model_names)
aictab(models, modnames = model_names)

# NOT CONSIDERING CORRELATIONS
# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m3 5 69.48       0.00   0.65   0.65 -29.07
# m5 1 72.29       2.81   0.16   0.81 -35.11
# m2 5 73.21       3.73   0.10   0.91 -30.94
# m4 2 74.33       4.85   0.06   0.97 -35.04
# m1 4 75.51       6.03   0.03   1.00 -33.32

# CONSIDERING CORRELATIONS
# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m3 5 69.48       0.00   0.61   0.61 -29.07
# m5 1 72.29       2.81   0.15   0.76 -35.11
# m2 5 73.21       3.73   0.09   0.86 -30.94
# m1 2 73.32       3.83   0.09   0.95 -34.53
# m4 2 74.33       4.85   0.05   1.00 -35.04


# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m3 2 65.53       0.00   0.90   0.90 -30.64
# m2 4 72.00       6.47   0.04   0.94 -31.56
# m5 1 72.29       6.77   0.03   0.97 -35.11
# m1 2 73.32       7.79   0.02   0.99 -34.53
# m4 2 74.33       8.80   0.01   1.00 -35.04

ggplot(invert,
       aes(x = InvertPesticideDetection,
           y = PrecipitationAmount_7days)) + geom_boxplot() + my_theme

#------------------------------------------------------------------------------#

# STAGE III: model average
model_avg <- model.avg(m3, m2, m5)

summary(model_avg)
confint(model_avg)

# Model-averaged coefficients:  
#   (full average) 
#                            Estimate Std. Error Adjusted SE z value Pr(>|z|)  
# (Intercept)               -0.376024   0.453721    0.462008   0.814   0.4157  
# PrecipitationAmount_7days -1.081054   0.494948    0.505149   2.140   0.0323 *
# PermanencePermanent        0.074846   0.441503    0.444677   0.168   0.8663  
# PermanenceSeasonal         0.008447   0.266269    0.273126   0.031   0.9753  
# PermanenceSemipermanent    0.066761   0.412826    0.416514   0.160   0.8727  

#                                2.5 %     97.5 %
# (Intercept)               -1.2815427  0.5294941
# PrecipitationAmount_7days -2.0374959 -0.2096799
# PermanencePermanent       -0.3929666  4.3472920
# PermanenceSeasonal        -2.4945680  2.9408551
# PermanenceSemipermanent   -0.7080897  4.2352669

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
df$term <- factor(
  df$term,
  levels = df$term[order(df$estimate)]
)

# add reference levels to plot
ref_levels <- data.frame(
  term = c("Wetland Permanence (Temporary)"), 
  estimate = 0,
  conf_low = 0,
  conf_high = 0
)

df <- rbind(df, ref_levels)

# relabel terms
df$term <- factor(
  df$term,
  levels = rev(c("PermanencePermanent",
                 "PermanenceSemipermanent",
                 "PermanenceSeasonal",
                 "Wetland Permanence (Temporary)",
                 "BufferedY",
                 "No Buffer Presence",
                 "AnnualPrecipitation_in",
                 "AnnualSnowfall_in",
                 "PrecipitationAmount_7days",
                 "PercentAg",
                 "NearestCropDistance_m")),
  labels = rev(c("Wetland Permanence (Permanent)",
                 "Wetland Permanence (Semi-permanent)",
                 "Wetland Permanence (Temporary/Seasonal)",
                 "Buffer Presence",
                 "No Buffer Presence",
                 "Annual Precipitation",
                 "Annual Snowfall",
                 "Precipitation Amount (Last 7 Days)",
                 "% Surrounding Cropland Cover",
                 "Nearest Crop Distance"
                 )))

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


# plot relationships
ggplot(invert, aes(y = LogPrecipAmount,
                  x = InvertPesticideDetection)) +
  geom_point()
