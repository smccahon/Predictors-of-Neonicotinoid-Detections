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
# combine temporary and seasonal permanence classes
# not necessary here, but this keeps it consistent across other analyses
invert <- invert %>% 
  mutate(Permanence = case_when(
    Permanence %in% c("Temporary", "Seasonal") ~ "Temporary/Seasonal",
    Permanence == "Semipermanent" ~ "Semipermanent",
    Permanence == "Permanent" ~ "Permanent"))

# convert characters to factors
invert <- invert %>% 
  mutate_if(is.character, as.factor)

# filter data (n = 51)
invert <- invert %>%
  filter(!is.na(InvertPesticideDetection))

# standardize data
invert.cs <- invert %>%
  mutate(across(where(is.numeric), scale))

# set reference levels
invert.cs$DominantCrop <- relevel(invert.cs$DominantCrop,
                                 ref = "Grassland")

invert.cs$Permanence <- relevel(invert.cs$Permanence,
                               ref = "Temporary/Seasonal")

#------------------------------------------------------------------------------#
#                   identify best model for each hypothesis                 ----                        
#------------------------------------------------------------------------------# 

# include all variables that show up in models with support (delta AICc < 2)

# best agricultural model
global.model <- glm(InvertPesticideDetection ~ PercentAg + DominantCrop +
                      NearestCropDistance_m + Buffered,
                    family = "binomial",
                    data = invert.cs,
                    na.action = na.fail)

# check for collinearity
car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#         (Int) Bff DmC  NCD_m     PrA df  logLik AICc delta weight
# 6  -9.667e-01   +     -1.913          3 -31.504 69.5  0.00  0.276
# 5  -2.851e-01         -0.809          2 -32.829 69.9  0.39  0.227
# 13 -3.381e-01         -1.416 -0.6087  3 -31.836 70.2  0.66  0.198
# 14 -9.223e-01   +     -2.239 -0.4024  4 -31.151 71.2  1.65  0.121
# 1  -1.967e-01                         1 -35.105 72.3  2.77  0.069
# 9  -1.981e-01                 0.1622  2 -34.944 74.1  4.62  0.027
# 2  -1.112e-01   +                     2 -34.993 74.2  4.72  0.026
# 7  -5.236e-03       + -1.061          6 -30.547 75.0  5.48  0.018
# 8  -7.927e-01   +   + -1.994          7 -29.651 75.9  6.39  0.011
# 10 -1.769e-01   +             0.1376  3 -34.941 76.4  6.87  0.009
# 15 -3.283e-01       + -1.455 -0.6323  7 -29.975 76.6  7.04  0.008
# 3  -2.231e-01       +                 5 -33.330 78.0  8.47  0.004
# 16 -9.768e-01   +   + -2.267 -0.4958  8 -29.322 78.1  8.55  0.004
# 11 -1.291e-02       +         0.3049  6 -33.086 80.1 10.56  0.001
# 4  -7.573e-16   +   +                 6 -33.195 80.3 10.78  0.001
# 12  3.796e-02   +   +         0.2588  7 -33.073 82.7 13.23  0.000
# Models ranked by AICc(x) 

# can remove dominant crop
m1 <- glm(InvertPesticideDetection ~ PercentAg + 
            NearestCropDistance_m + 
            Buffered,
          family = "binomial",
          data = invert.cs)

summary(m1)
confint(m1)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m1) 
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m1)$PercentAg) #  pattern
plotResiduals(simulationOutput, form = model.frame(m1)$NearestCropDistance_m) # good
plotResiduals(simulationOutput, form = model.frame(m1)$Buffered)  # good

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
#      (Int) Dst_Cls_Wtl_m Prm df  logLik AICc delta weight
# 3 -1.4660                 +  3 -31.577 69.7  0.00  0.592
# 4 -1.4630      -0.08867   +  4 -31.537 71.9  2.28  0.190
# 1 -0.1967                    1 -35.105 72.3  2.63  0.159
# 2 -0.1977      -0.12230      2 -35.015 74.3  4.62  0.059
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

# best precipitation model
global.model <- glm(InvertPesticideDetection ~ AnnualSnowfall_in + 
                      AnnualPrecipitation_in + DaysSinceLastPrecipitation_5mm +
                      PrecipitationAmount_7days,
                    family = "binomial",
                    na.action = na.fail,
                    data = invert.cs)

# check for collinearity
car::vif(global.model) # vif < 3

dredge(global.model)

# Model selection table 
#      (Int)   AnP_in AnS_in  DSL_5mm PrA_7dy df  logLik AICc delta weight
# 9  -0.3275                           -1.124  2 -30.638 65.5  0.00  0.269
# 10 -0.3070  0.46270                  -1.342  3 -29.855 66.2  0.70  0.190
# 11 -0.3587          0.4131           -1.322  3 -29.908 66.3  0.80  0.180
# 13 -0.3380                 -0.06824  -1.196  3 -30.626 67.8  2.24  0.088
# 12 -0.3291  0.30240 0.2388           -1.382  4 -29.707 68.3  2.76  0.068
# 14 -0.3325  0.48330        -0.16680  -1.527  4 -29.791 68.5  2.93  0.062
# 15 -0.3703          0.4137 -0.06936  -1.399  4 -29.896 68.7  3.14  0.056
# 5  -0.1946                  0.63200          2 -33.017 70.3  4.76  0.025
# 16 -0.3493  0.33070 0.2232 -0.13610  -1.534  5 -29.665 70.7  5.14  0.021
# 7  -0.1917          0.2842  0.70460          3 -32.574 71.7  6.13  0.013
# 1  -0.1967                                   1 -35.105 72.3  6.77  0.009
# 6  -0.1915  0.09328         0.65760          3 -32.971 72.5  6.93  0.008
# 8  -0.1966 -0.08974 0.3328  0.69250          4 -32.546 74.0  8.44  0.004
# 3  -0.1978          0.1462                   2 -34.974 74.2  8.67  0.004
# 2  -0.1967 -0.02760                          2 -35.100 74.5  8.93  0.003
# 4  -0.1990 -0.15020 0.2275                   3 -34.877 76.3 10.74  0.001
# Models ranked by AICc(x) 

# can remove days since last precipitation event
m3 <- glm(InvertPesticideDetection ~  AnnualSnowfall_in + 
            AnnualPrecipitation_in +
            PrecipitationAmount_7days,
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
# m3 4 68.28       0.00   0.52   0.52 -29.71
# m2 3 69.67       1.38   0.26   0.78 -31.58
# m1 4 71.17       2.89   0.12   0.90 -31.15
# m5 1 72.29       4.01   0.07   0.97 -35.11
# m4 2 74.33       6.05   0.03   1.00 -35.04

#------------------------------------------------------------------------------#

# STAGE III: model average
model_avg <- model.avg(m3, m2, m1)

summary(model_avg)
confint(model_avg)

# Model-averaged coefficients:  
#   (full average) 
#                           Estimate Std. Error Adjusted SE z value Pr(>|z|)
# (Intercept)               -0.73769    0.71761     0.72719   1.014    0.310
# AnnualSnowfall_in          0.13749    0.35550     0.36387   0.378    0.706
# AnnualPrecipitation_in     0.17411    0.39761     0.40665   0.428    0.669
# PrecipitationAmount_7days -0.79545    0.79111     0.79649   0.999    0.318
# PermanencePermanent        0.52903    0.92897     0.93382   0.567    0.571
# PermanenceSemipermanent    0.46743    0.86405     0.87032   0.537    0.591
# PercentAg                 -0.05463    0.22480     0.22853   0.239    0.811
# NearestCropDistance_m     -0.30396    0.96688     0.97643   0.311    0.756
# BufferedY                  0.19189    0.69102     0.70038   0.274    0.784

#                                 2.5 %     97.5 %
# (Intercept)               -2.16295095  0.6875656
# AnnualSnowfall_in         -0.65027064  1.1278519
# AnnualPrecipitation_in    -0.67447914  1.2792486
# PrecipitationAmount_7days -2.44088167 -0.3221639
# PermanencePermanent        0.27885110  3.3892726
# PermanenceSemipermanent   -0.08533193  3.3263074
# PercentAg                 -1.37194266  0.5671668
# NearestCropDistance_m     -5.45372384  0.9762370
# BufferedY                 -1.27870197  4.1053824

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
  term = c("Wetland Permanence (Temporary/Seasonal)",
           "No Buffer Presence"), 
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
                 "Wetland Permanence (Temporary/Seasonal)",
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
