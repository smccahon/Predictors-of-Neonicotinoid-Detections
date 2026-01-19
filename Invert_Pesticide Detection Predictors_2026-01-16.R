#-----------------------------------#
#  Predictors of Neonics in Inverts #
#     Analysis by Shelby McCahon    #
#        Created: 2026-01-16        #
#       Modified: 2026-01-16        #
#-----------------------------------#

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

invert <- invert %>%
  mutate(Buffered = ifelse(Buffered == "Y", 1, 0))

# convert characters to factors
invert <- invert %>% 
  mutate_if(is.character, as.factor)

# filter data (n = 51)
invert <- invert %>%
  filter(!is.na(InvertPesticideDetection))

# log transform nearest crop distance due to skew
invert <- invert %>% 
  mutate(CropDistance = NearestCropDistance_m + 0.0001) %>% 
  mutate(CropDistance = log(CropDistance))

# log transform precipitation amount due to skew
invert <- invert %>% 
  mutate(Precip = PrecipitationAmount_7days + 0.0001) %>% 
  mutate(Precip = log(Precip))

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
# 
 cor_mat <- cor(invert[sapply(invert, is.numeric)], 
                use = "pairwise.complete.obs")
# 
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

 #                               Var1                      Var2 Correlation
 # 156                       Buffered                 PercentAg  -0.7439067
 # 151          NearestCropDistance_m                 PercentAg  -0.7261301
 # 72                        Buffered     NearestCropDistance_m   0.6990801
 # 288 DaysSinceLastPrecipitation_5mm PrecipitationAmount_7days  -0.6370290

#------------------------------------------------------------------------------#
#                   identify best model for each hypothesis                 ----                        
#------------------------------------------------------------------------------# 

# crop distance has the best fit (AICc wt. = 0.57)
# m1 <- glm(InvertPesticideDetection ~ PercentAg + DominantCrop,
#                   family = "binomial",
#                   data = invert.cs,
#                   na.action = na.fail)
#  
# m2 <- glm(InvertPesticideDetection ~ Buffered + DominantCrop,
#          family = "binomial",
#          data = invert.cs,
#          na.action = na.fail)
# 
# m3 <- glm(InvertPesticideDetection ~ CropDistance + DominantCrop,
#          family = "binomial",
#          data = invert.cs,
#          na.action = na.fail)
#  
# model_names <- paste0("m", 1:3)
# models <- mget(model_names)
# aictab(models, modnames = model_names)


global.model <- glm(InvertPesticideDetection ~ CropDistance + DominantCrop,
                family = "binomial",
                data = invert.cs,
                na.action = na.fail)
 
# check for collinearity
car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#    (Intrc)  CrpDs DmnnC df  logLik AICc delta weight
# 1 -0.1967               1 -35.105 72.3  0.00  0.585
# 2 -0.2067 0.3249        2 -34.533 73.3  1.02  0.351
# 3 -0.2231            +  5 -33.330 78.0  5.70  0.034
# 4 -0.5269 0.5829     +  6 -32.149 78.2  5.92  0.030
# Models ranked by AICc(x) 

# can remove dominant crop
m1 <- glm(InvertPesticideDetection ~ CropDistance,
          family = "binomial",
          data = invert.cs,
          na.action = na.fail)

summary(m1)
confint(m1)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m1) 
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m1)$CropDistance) # good

#------------------------------------------------------------------------------#

# best wetland dynamics model
global.model <- glm(InvertPesticideDetection ~ Permanence + 
                      Dist_Closest_Wetland_m,
                    family = "binomial",
                    na.action = na.fail,
                    data = invert.cs)

# check for collinearity
car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#     (Int) Dst_Cls_Wtl_m Prm df  logLik AICc delta weight
# 3 -1.6090                 +  4 -31.563 72.0  0.00  0.402
# 1 -0.1967                    1 -35.105 72.3  0.30  0.347
# 2 -0.1977      -0.12230      2 -35.015 74.3  2.28  0.128
# 4 -1.6450      -0.09826   +  5 -31.514 74.4  2.37  0.123
# Models ranked by AICc(x) 

# distance to closest wetland can be removed
m2 <- glm(InvertPesticideDetection ~ Permanence,
          data = invert.cs,
          family = "binomial",
          na.action = na.fail)

summary(m2)
confint(m2)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$Permanence) #  good

#------------------------------------------------------------------------------#

# precipitation amount performed better than days since event (AICc wt = 0.81)
#  m1 <- glm(InvertPesticideDetection ~ Precip,
#            family = "binomial",
#            data = invert.cs,
#            na.action = na.fail) 
#  
#  m2 <- glm(InvertPesticideDetection ~ DaysSinceLastPrecipitation_5mm,
#            family = "binomial",
#            data = invert.cs,
#            na.action = na.fail)
# 
#  model_names <- paste0("m", 1:2)
#  models <- mget(model_names)
#  aictab(models, modnames = model_names)
# 
# # best precipitation model
# global.model <- glm(InvertPesticideDetection ~ AnnualSnowfall_in + 
#                       AnnualPrecipitation_in + PrecipitationAmount_7days,
#                     family = "binomial",
#                     na.action = na.fail,
#                     data = invert.cs)

m3 <- glm(InvertPesticideDetection ~ PrecipitationAmount_7days,
          family = "binomial",
          na.action = na.fail,
          data = invert.cs)

summary(m3)
confint(m3)

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

# is log-transformation better? no
# m1 <- glm(InvertPesticideDetection ~ PrecipitationAmount_7days,
#           family = "binomial",
#           na.action = na.fail,
#           data = invert.cs)
# 
# m2 <- glm(InvertPesticideDetection ~ Precip,
#           family = "binomial",
#           na.action = na.fail,
#           data = invert.cs)
# 
# 
#   model_names <- paste0("m", 1:2)
#   models <- mget(model_names)
#   aictab(models, modnames = model_names)

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
model_avg <- model.avg(m3, m2)

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
ggplot(invert, aes(y = PrecipitationAmount_7days,
                  x = InvertPesticideDetection)) +
  geom_boxplot() + geom_point()

cor(invert$NearestCropDistance_m, invert$PercentAg) # -0.73
