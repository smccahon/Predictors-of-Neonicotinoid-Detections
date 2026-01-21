#---------------------------------#
# Predictors of Neonics in Plasma #
#   Analysis by Shelby McCahon    #
#      Created: 2026-01-13        #
#     Modified: 2026-01-20        #
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
library(lme4)

# DETAILS OF ANALYSIS
# STAGE I: Identify the best predictors within each hypothesis (< 2 delta AICc)
# STAGE II: Identify which hypotheses are supported
# STAGE III: Model average the top 95% confidence set

# ANALYSIS NOTES: 
# I included sampling event in all models due to VERY strong temporal effect.
# Event became informed null.
# I considered including site as a random effect but it created model
# convergence issues and singular fit warnings (~0 variation).
# I used "Migratory Status" rather than "Species" because it had a much stronger
# effect according to AICc (AICc wt = 0.99).
# I used Sampling Event rather than Season + Year because it had much stronger
# support (AICc wt = 1.00).
# I log-transformed nearest cropland distance to
# account for heavy right-skew in the data. There was stronger support for
# log-transformation of cropland distance. I did the same for wetland distance
# but model without log-transformation had stronger support so I did not
# include the transformation. 
# There is no correlation in any of the ag variables.

# All models had VIF < 3

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

# convert factors to binary for correlation testing
 birds <- birds %>%
   mutate(Season = ifelse(Season == "Spring", 1, 0)) %>% 
   mutate(MigStatus = ifelse(MigStatus == "Migrant", 1, 0)) %>% 
   mutate(Sex = ifelse(Sex == "F", 1, 0))
# 
# birds$Year <- as.factor(birds$Year)

# combine temporary and seasonal permanence classes
birds <- birds %>% 
  mutate(Permanence = case_when(
    Permanence %in% c("Temporary", "Seasonal") ~ "Temporary/Seasonal",
    Permanence == "Semipermanent" ~ "Semipermanent",
    Permanence == "Permanent" ~ "Permanent"))

# log transform nearest crop distance due to skew
birds <- birds %>% 
  mutate(CropDistance = NearestCropDistance_m + 0.0001) %>% 
  mutate(LogCropDistance = log(CropDistance))

# log transform nearest wetland distance due to skew
birds <- birds %>% 
  mutate(WetlandDistance = Dist_Closest_Wetland_m + 0.0001) %>% 
  mutate(LogWetlandDistance = log(WetlandDistance))

# convert characters to factors
birds <- birds %>% 
  mutate_if(is.character, as.factor)

# convert integers to numeric
birds <- birds %>% 
  mutate_if(is.integer, as.numeric)

# filter data (n = 171)
birds <- birds %>%
  filter(!is.na(PlasmaDetection)) %>% 
  filter(!is.na(Sex))

#------------------------------------------------------------------------------#
#                         identify correlations                             ----                        
#------------------------------------------------------------------------------# 

#  cor_mat <- cor(birds[sapply(birds, is.numeric)], 
#                 use = "pairwise.complete.obs")
#  
#  # Convert matrix to a data frame of pairs
#  cor_df <- as.data.frame(as.table(cor_mat))
#  
# # # Rename columns
#  names(cor_df) <- c("Var1", "Var2", "Correlation")
# # 
# # # Convert factors to characters
#  cor_df$Var1 <- as.character(cor_df$Var1)
#  cor_df$Var2 <- as.character(cor_df$Var2)
# # 
#  # Keep only unique pairs (upper triangle)
#  cor_df <- cor_df[cor_df$Var1 < cor_df$Var2, ]
# # 
#  # Filter by threshold
#  high_corr <- subset(cor_df, abs(Correlation) > 0.6)
#  
#  # Sort by correlation strength
#  high_corr <- high_corr[order(-abs(high_corr$Correlation)), ]
#  high_corr

# relevant correlations to check more for collinearity

#                               Var1                      Var2 Correlation
# 10                          Julian                    Season  -0.9593386
# 436                         Season                      SPEI   0.7964145
# 770 DaysSinceLastPrecipitation_5mm PrecipitationAmount_7days  -0.7586376
# 445                         Julian                      SPEI  -0.6736094
# 305 DaysSinceLastPrecipitation_5mm                    Julian   0.6277669

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
# model without site as a random effect has higher explanatory power (AICc wt = 0.77)
# including site also caused singular fit model warnings

# m1 <- glmmTMB(PlasmaDetection ~ Event + PercentAg + DominantCrop + 
#               NearestCropDistance_m + Dist_Closest_Wetland_m +
#               Permanence + SPEI + 
#               time_hours + Sex + MigStatus,
#               family = "binomial",
#               data = birds.cs,
#               REML = TRUE)
#  
# 
# m2 <- glmmTMB(PlasmaDetection ~ Event + PercentAg + DominantCrop + 
#                NearestCropDistance_m + Dist_Closest_Wetland_m +
#                Permanence + SPEI + 
#                time_hours + Sex + MigStatus + (1|Site),
#                family = "binomial",
#                data = birds.cs,
#                REML = TRUE)
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)

#------------------------------------------------------------------------------#
#                   identify best model for each hypothesis                 ----                        
#------------------------------------------------------------------------------# 

# include all variables that show up in models with support (delta AICc < 2)
# log-transforming crop distance improved model fit

# best agricultural model
global.model <- glm(PlasmaDetection ~ LogCropDistance + Event + PercentAg +
                    DominantCrop,
                    family = "binomial",
                    data = birds.cs,
                    na.action = na.fail)


# check for collinearity
car::vif(global.model) # vif < 4

dredge(global.model)

# Model selection table 
#    (Intrc) DmnnC Event   LgCrD    PrcnA df   logLik  AICc delta weight
# 7  -0.5162           + -0.5918           5  -75.327 161.0  0.00  0.520
# 15 -0.3873           + -0.7086 -0.24750  6  -74.878 162.3  1.25  0.278
# 8   0.3153     +     + -0.8985           8  -74.020 164.9  3.91  0.073
# 4  -0.7119     +     +                   7  -75.761 166.2  5.19  0.039
# 3  -0.8109           +                   4  -79.036 166.3  5.30  0.037
# 16  0.5825     +     + -0.9227  0.17870  9  -73.941 167.0  5.98  0.026
# 11 -0.8355           +          0.09026  5  -78.951 168.3  7.25  0.014
# 12 -0.5834     +     +          0.08506  8  -75.743 168.4  7.36  0.013
# 14  1.4560     +       -0.9379  0.60220  6  -91.038 194.6 33.57  0.000
# 6   0.7804     +       -0.9091           5  -93.142 196.6 35.63  0.000
# 10  0.9559     +                0.56870  5  -95.877 202.1 41.10  0.000
# 2   0.3302     +                         4  -97.967 204.2 43.16  0.000
# 13 -0.6270             -0.5668 -0.48390  3 -105.662 217.5 56.45  0.000
# 5  -0.6001             -0.3385           2 -109.070 222.2 61.19  0.000
# 1  -0.5896                               1 -111.408 224.8 63.82  0.000
# 9  -0.5954                     -0.20840  2 -110.553 225.2 64.16  0.000
# Models ranked by AICc(x) 

# remove dominant crop
m1 <- glm(PlasmaDetection ~ LogCropDistance + Event + PercentAg,
          family = "binomial",
          data = birds.cs)

summary(m1)
confint(m1)

# model validation --> no issues
simulationOutput <- simulateResiduals(fittedModel = m1, n = 1000) 
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m1)$DominantCrop) #  good
plotResiduals(simulationOutput, form = model.frame(m1)$LogCropDistance) # good
plotResiduals(simulationOutput, form = model.frame(m1)$Event)  # good
plotResiduals(simulationOutput, form = model.frame(m1)$PercentAg)  # good

#------------------------------------------------------------------------------#

# best wetland dynamics model
# model without log-transformation is better
global.model <- glm(PlasmaDetection ~ Event + Dist_Closest_Wetland_m +
                      Permanence + Porosity,
                    family = "binomial",
                    data = birds.cs,
                    na.action = na.fail)

car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#     (Int) Dst_Cls_Wtl_m Evn Prm df   logLik  AICc delta weight
# 4 -0.7169       0.62340   +      5  -75.328 161.0  0.00  0.555
# 8 -2.1330       1.01700   +   +  7  -73.510 161.7  0.69  0.394
# 3 -0.8109                 +      4  -79.036 166.3  5.29  0.039
# 7 -0.1073                 +   +  6  -78.119 168.7  7.73  0.012
# 5  0.3254                     +  3 -106.308 218.8 57.74  0.000
# 6  0.2934       0.07419       +  4 -106.210 220.7 59.64  0.000
# 1 -0.5896                        1 -111.408 224.8 63.82  0.000
# 2 -0.5940       0.18020          2 -110.766 225.6 64.58  0.000
# Models ranked by AICc(x) 

# keep all predictors
m2 <- glm(PlasmaDetection ~ Permanence + Event + Dist_Closest_Wetland_m,
                    family = "binomial",
                    data = birds.cs)

car::vif(m2)

# model validation --> no severe issues
simulationOutput <- simulateResiduals(fittedModel = m2, n =1000) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$Permanence) #  good
plotResiduals(simulationOutput, form = model.frame(m2)$Event)  # good
plotResiduals(simulationOutput, form = model.frame(m2)$Dist_Closest_Wetland_m)  # good

#------------------------------------------------------------------------------#

# best temporal model
# Julian was considered but had high VIF (4.46) and high correlation
global.model <- glm(PlasmaDetection ~ Event + time_hours,
                    data = birds.cs,
                    family = "binomial",
                    na.action = na.fail)

car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#     (Int) Evn scn_snc_mdn df   logLik  AICc delta weight
# 2 -0.8109   +              4  -79.036 166.3  0.00  0.645
# 4 -0.8386   +     -0.1953  5  -78.572 167.5  1.19  0.355
# 1 -0.5896                  1 -111.408 224.8 58.53  0.000
# 3 -0.5954          0.2016  2 -110.623 225.3 59.00  0.000
# Models ranked by AICc(x) 

# keep both predictors
m3 <- glm(PlasmaDetection ~ Event + time_hours,
          family = "binomial",
          data = birds.cs)

# model validation --> no issues
simulationOutput <- simulateResiduals(fittedModel = m3, n = 1000) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$time_hours) #  good
plotResiduals(simulationOutput, form = model.frame(m3)$Event)  # good
plotResiduals(simulationOutput, form = model.frame(m3)$Julian)  # good

#------------------------------------------------------------------------------#

# drought model
m4 <- glm(PlasmaDetection ~ SPEI + Event,
          data = birds.cs,
          family = "binomial",
          na.action = na.fail)

car::vif(m4) # vif < 3

summary(m4)
confint(m4)

# model validation --> no issues
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m4)$Event)  # good
plotResiduals(simulationOutput, form = model.frame(m4)$SPEI)  # good

#------------------------------------------------------------------------------#

# life history model
global.model <- glm(PlasmaDetection ~ Sex + MigStatus + Event,
          data = birds.cs,
          family = "binomial",
          na.action = na.fail)

car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#   (Intrc) Event MgStt Sex df   logLik  AICc delta weight
# 2 -0.8109     +            4  -79.036 166.3  0.00  0.350
# 6 -0.5332     +         +  5  -77.978 166.3  0.01  0.349
# 4 -0.8891     +     +      5  -78.795 168.0  1.64  0.154
# 8 -0.6081     +     +   +  6  -77.766 168.0  1.73  0.147
# 7 -0.4319           +   +  3 -105.934 218.0 51.70  0.000
# 3 -0.8550           +      2 -108.144 220.4 54.05  0.000
# 5 -0.1452               +  2 -108.538 221.1 54.83  0.000
# 1 -0.5896                  1 -111.408 224.8 58.53  0.000
# Models ranked by AICc(x) 

summary(m5)
confint(m5) 

# keep all predictors
m5 <- glm(PlasmaDetection ~ Sex + MigStatus + Event,
                    data = birds.cs,
                    family = "binomial",
                    na.action = na.fail)

# model validation --> no issues
simulationOutput <- simulateResiduals(fittedModel = m5) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m5)$MigStatus) #  good
plotResiduals(simulationOutput, form = model.frame(m5)$Event)  # good
plotResiduals(simulationOutput, form = model.frame(m5)$Sex)  # good

#------------------------------------------------------------------------------#

# best pesticide detection model
m6 <- glm(PlasmaDetection ~ EnvDetection + Event,
          data = birds.cs,
          family = "binomial")

car::vif(m6) # vif < 2

summary(m6)
confint(m6) 

# model validation --> no issues
simulationOutput <- simulateResiduals(fittedModel = m6) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m6)$EnvDetection) #  good
plotResiduals(simulationOutput, form = model.frame(m6)$Event)  # good

#------------------------------------------------------------------------------#

# informed null
m7 <- glm(PlasmaDetection ~ Event,
          data = birds.cs,
          family = "binomial")

#------------------------------------------------------------------------------#

# precipitation?
global.model <- glm(PlasmaDetection ~ DaysSinceLastPrecipitation_5mm + 
                      PrecipitationAmount_7days + AnnualSnowfall_in +
                      AnnualPrecipitation_in + Event,
                    data = birds.cs,
                    family = "binomial",
                    na.action = na.fail)

car::vif(global.model) # vif < 2

dredge(global.model)

# retain all

m8 <- glm(PlasmaDetection ~ DaysSinceLastPrecipitation_5mm + 
                      PrecipitationAmount_7days + AnnualSnowfall_in +
                      AnnualPrecipitation_in + Event,
                    data = birds.cs,
                    family = "binomial",
                    na.action = na.fail)


#------------------------------------------------------------------------------#

# STAGE II: Cross-hypothesis

# AIC model selection
model_names <- paste0("m", 1:8)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K   AICc Delta_AICc AICcWt Cum.Wt     LL
# m2 7 161.71       0.00   0.49   0.49 -73.51
# m1 6 162.27       0.56   0.37   0.87 -74.88
# m7 4 166.31       4.61   0.05   0.92 -79.04
# m3 5 167.51       5.80   0.03   0.94 -78.57
# m5 6 168.04       6.34   0.02   0.97 -77.77
# m6 5 168.42       6.71   0.02   0.98 -79.03
# m4 5 168.43       6.72   0.02   1.00 -79.03

#------------------------------------------------------------------------------#

# STAGE III: model average
model_avg <- model.avg(m1, m2, m3, m7, m8)

summary(model_avg)
confint(model_avg)


# Model-averaged coefficients:  
#   (full average) 
#                          Estimate Std. Error Adjusted SE z value Pr(>|z|)   
# (Intercept)             -1.336180   1.197225    1.201722   1.112  0.26619   
# PermanencePermanent      0.719052   1.014281    1.018378   0.706  0.48014   
# PermanenceSemipermanent  0.899308   1.106704    1.109993   0.810  0.41783   
# EventFall 2023          -0.885278   0.653444    0.657592   1.346  0.17822   
# EventSpring 2022         4.854519   1.679124    1.687329   2.877  0.00401 **
# EventSpring 2023        -0.903842   0.751638    0.756321   1.195  0.23207   
# Dist_Closest_Wetland_m   0.532546   0.566288    0.567111   0.939  0.34770   
# LogCropDistance         -0.280216   0.384098    0.384629   0.729  0.46629   
# PercentAg               -0.097865   0.202756    0.203722   0.480  0.63095   
# time_hours              -0.005622   0.047696    0.047882   0.117  0.90653 

#                              2.5 %     97.5 %
# (Intercept)             -3.6915111  1.0191509
# PermanencePermanent     -0.6651340  3.4124209
# PermanenceSemipermanent -0.1895411  3.6255313
# EventFall 2023          -2.1741350  0.4035795
# EventSpring 2022         1.5474149  8.1616228
# EventSpring 2023        -2.3862027  0.5785197
# Dist_Closest_Wetland_m   0.3350684  1.6996345
# LogCropDistance         -1.2291730 -0.1880566
# PercentAg               -0.7582946  0.2633313
# time_hours              -0.5997960  0.2092773



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
  #  df$term,
  #   levels = df$term[order(df$estimate)]
  # )

# add reference levels to plot
ref_levels <- data.frame(
  term = c("Sampling Event (Fall 2021)", 
           "Wetland Permanence (Temporary/Seasonal)"),
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
             "time_hours",
             "PercentAg",
             "LogCropDistance")),
  labels = rev(c("Sampling Event (Spring 2022)",
             "Sampling Event (Spring 2023)",
             "Sampling Event (Fall 2023)",
             "Sampling Event (Fall 2021)",
             "Wetland Permanence (Semi-permanent)",
             "Wetland Permanence (Permanent)",
             "Wetland Permanence (Temporary/Seasonal)",
             "Distance to Nearest Wetland",
             "Capture Time",
             "% Surrounding Cropland Cover",
             "Log(Distance to Nearest Crop)")))

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
  
ggsave(filename = "Model-averaged Plasma Neonic Results_2026-01-20.tiff",
       path = "C:/Users/shelb/OneDrive - University of Idaho/Masters Project/Analysis/Predictors of Neonics/figures", 
       width = 13, height = 8, units = "in", device='tiff', dpi=300)

#------------------------------------------------------------------------------#


# interpreting results
# what's so interesting about isolated wetlands?

# not related to surrounding cropland
ggplot(birds,
       aes(x = Dist_Closest_Wetland_m,
           y = PercentAg)) + geom_point() + my_theme

ggplot(birds,
       aes(x = Event,
           y = Dist_Closest_Wetland_m)) + geom_boxplot() + my_theme

ggplot(birds,
       aes(x = PlasmaDetection,
           y = Dist_Closest_Wetland_m)) + geom_boxplot() + my_theme


