#-----------------------------------#
#  Predictors of Neonics in Inverts #
#     Analysis by Shelby McCahon    #
#        Created: 2026-01-22        #
#       Modified: 2026-01-23        #
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

# variables with correlation > 0.80 not considered in the same model

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

# *agriculture candidate models ----

m1 <- glm(InvertPesticideDetection ~ PercentAg, 
          data = invert.cs,
          family = "binomial")

# log transformation improves fit (wt. = 0.55)
m2 <- glm(InvertPesticideDetection ~ LogCropDistance, 
          data = invert.cs,
          family = "binomial")

m3 <- glm(InvertPesticideDetection ~ DominantCrop, 
          data = invert.cs,
          family = "binomial")

m4 <- glm(InvertPesticideDetection ~ Buffered, 
          data = invert.cs,
          family = "binomial")

m5 <- glm(InvertPesticideDetection ~ PercentAg + LogCropDistance, 
          data = invert.cs,
          family = "binomial")

m6 <- glm(InvertPesticideDetection ~ Buffered + PercentAg, 
          data = invert.cs,
          family = "binomial")

# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m5 3 73.29       0.00   0.28   0.28 -33.39
# m2 2 73.32       0.03   0.28   0.56 -34.53
# m1 2 74.14       0.85   0.18   0.74 -34.94
# m4 2 74.24       0.95   0.17   0.91 -34.99
# m6 3 76.39       3.10   0.06   0.97 -34.94
# m3 5 77.99       4.70   0.03   1.00 -33.33


model_names <- paste0("m", 1:6)
models <- mget(model_names)
aictab(models, modnames = model_names)

m.ag <- glm(InvertPesticideDetection ~ LogCropDistance, 
          data = invert.cs,
          family = "binomial")


# model validation --> some pattern, not extreme
simulationOutput <- simulateResiduals(fittedModel = m.ag, n = 2000) 
plot(simulationOutput)
testDispersion(m.ag) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, 
              form = model.frame(m.ag)$LogCropDistance) # pattern

#------------------------------------------------------------------------------#

# *hydrology candidate models ----

m1 <- glm(InvertPesticideDetection ~ Permanence,
          family = "binomial",
          data = invert.cs)

m2 <- glm(InvertPesticideDetection ~ Porosity,
          family = "binomial",
          data = invert.cs)

m3 <- glm(InvertPesticideDetection ~ LogWetlandDistance,
          family = "binomial",
          data = invert.cs)

m4 <- glm(InvertPesticideDetection ~ Permanence + LogWetlandDistance,
          family = "binomial",
          data= invert.cs)

model_names <- paste0("m", 1:4)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m1 4 72.00       0.00   0.34   0.34 -31.56
# m2 2 72.06       0.06   0.33   0.68 -33.91
# m4 5 73.21       1.22   0.19   0.86 -30.94
# m3 2 73.85       1.85   0.14   1.00 -34.80

m.hydrology <- glm(InvertPesticideDetection ~ Permanence,
                   family = "binomial",
                   data = invert.cs)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m.hydrology, 
                                      n = 1000) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.hydrology)$Permanence) #  good


#------------------------------------------------------------------------------#

# *water quality candidate models ----

m1 <- glm(InvertPesticideDetection ~ pH_probe,
            family = "binomial",
            data= invert.cs)

m2 <- glm(InvertPesticideDetection ~ WaterTemp,
          family = "binomial",
          data= invert.cs)

m3 <- glm(InvertPesticideDetection ~ pH_probe + WaterTemp,
          family = "binomial",
          data= invert.cs)

model_names <- paste0("m", 1:3)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m1 2 72.80       0.00   0.57   0.57 -34.28
# m2 2 74.45       1.64   0.25   0.81 -35.10
# m3 3 75.03       2.23   0.19   1.00 -34.26

m.waterquality <- glm(InvertPesticideDetection ~ pH_probe,
                      family = "binomial",
                      data= invert.cs)


#------------------------------------------------------------------------------#

# *precipitation candidate models ----

# log transformation supported (wt = 0.80)
m1 <- glm(InvertPesticideDetection ~ LogPrecipDays,
          family = "binomial",
          data = invert.cs)

# log transformation was not supported (wt = 0.28)
m2 <- glm(InvertPesticideDetection ~ PrecipitationAmount_7days,
          family = "binomial",
          data = invert.cs)

m3 <- glm(InvertPesticideDetection ~ AnnualSnowfall_in,
          family = "binomial",
          data = invert.cs)

m4 <- glm(InvertPesticideDetection ~ AnnualPrecipitation_in,
          family = "binomial",
          data = invert.cs)

model_names <- paste0("m", 1:4)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m2 2 65.53       0.00   0.78   0.78 -30.64
# m1 2 68.29       2.77   0.20   0.98 -32.02
# m3 2 74.20       8.67   0.01   0.99 -34.97
# m4 2 74.45       8.93   0.01   1.00 -35.10

m.precip <- glm(InvertPesticideDetection ~ PrecipitationAmount_7days,
                family = "binomial",
                data = invert.cs)

# model validation --> no severe issues
simulationOutput <- simulateResiduals(fittedModel = m.precip) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$PrecipitationAmount_7days) # good

#------------------------------------------------------------------------------#

# *temporal candidate models ----
m.time <- glm(InvertPesticideDetection ~ Season,
              data = invert.cs,
              family = "binomial",
              na.action = na.fail)


# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m.time) 
plot(simulationOutput)
testDispersion(m.time) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.time)$Season) #  good

#------------------------------------------------------------------------------#

# null model

m.null <- glm(InvertPesticideDetection ~ 1,
          family = "binomial",
          data = invert.cs)


#------------------------------------------------------------------------------#
#                 STAGE II: evaluate support for hypotheses                 ----                        
#------------------------------------------------------------------------------#

models <- list(
  Time        = m.time,
  Agriculture = m.ag,
  Precip      = m.precip,
  Hydrology   = m.hydrology,
  Null        = m.null,
  Water       = m.waterquality
)

aictab(cand.set = models)

# Model selection based on AICc:
#   
#             K  AICc Delta_AICc AICcWt Cum.Wt     LL
# Precip      2 65.53       0.00   0.88   0.88 -30.64
# Hydrology   4 72.00       6.47   0.03   0.92 -31.56
# Null        1 72.29       6.77   0.03   0.95 -35.11
# Water       2 72.80       7.28   0.02   0.97 -34.28
# Agriculture 2 73.32       7.79   0.02   0.99 -34.53
# Time        2 74.31       8.79   0.01   1.00 -35.03

confint(m.precip)
summary(m.precip)

# view results
ggplot(invert,
       aes(y = PrecipitationAmount_7days,
           x = InvertPesticideDetection)) +
  geom_boxplot() + my_theme

#------------------------------------------------------------------------------#
#                             PLOT TOP MODEL(s)                             ----                        
#------------------------------------------------------------------------------#


