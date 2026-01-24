#---------------------------------------------------#
# Predictors of Neonics in Lesser Yellowlegs Plasma #
#           Analysis by Shelby McCahon              #
#              Created: 2026-01-23                  #
#             Modified: 2026-01-23                  #
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


# ANALYSIS NOTES: 
# I included sampling event in all models due to VERY strong temporal effect
# Event became informed null
# I considered including site as a random effect but it created model
# convergence issues and singular fit warnings
# All models had VIF < 3

# ALL 11 birds IN SPRING 2022 HAVE DETECTIONS! will need to use
# Julian instead of sampling event (or year or season) because there is not 
# enough variation for model estimation

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

# log transform nearest crop distance due to right skew
leye <- leye %>% 
  mutate(CropDistance = NearestCropDistance_m + 0.0001) %>% 
  mutate(LogCropDistance = log(CropDistance))

# log transform nearest wetland distance due to right skew
leye <- leye %>% 
  mutate(WetlandDistance = Dist_Closest_Wetland_m) %>% 
  mutate(LogWetlandDistance = log(WetlandDistance))

# log transform precipitation amount due to right skew
leye <- leye %>% 
  mutate(LogPrecipAmount = PrecipitationAmount_7days + 0.0001) %>% 
  mutate(LogPrecipAmount = log(LogPrecipAmount))

# log transform precipitation event due to right skew
leye <- leye %>% 
  mutate(LogPrecipDays = DaysSinceLastPrecipitation_5mm) %>% 
  mutate(LogPrecipDays = log(LogPrecipDays))

# standardize data
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# set reference levels
leye.cs$DominantCrop <- relevel(leye.cs$DominantCrop,
                                 ref = "Grassland")

leye.cs$Permanence <- relevel(leye.cs$Permanence,
                               ref = "Temporary/Seasonal")

#------------------------------------------------------------------------------#
#                               correlations                                ----                        
#------------------------------------------------------------------------------# 
 
 #                               Var1                      Var2 Correlation
 # 11                          Julian                    Season  -0.9383480
 # 464         AnnualPrecipitation_in                      SPEI   0.8821235
 # 782 DaysSinceLastPrecipitation_5mm PrecipitationAmount_7days  -0.8493425
 # 436                         Season                      SPEI   0.7339089
 # 29          AnnualPrecipitation_in                    Season   0.7277824
 # 461              AnnualSnowfall_in                      SPEI   0.7013480
 # 608 DaysSinceLastPrecipitation_5mm    Dist_Closest_Wetland_m   0.6669624
 # 238          NearestCropDistance_m                 PercentAg  -0.6621564
 # 754         AnnualPrecipitation_in         AnnualSnowfall_in   0.6449298
 # 775         Dist_Closest_Wetland_m PrecipitationAmount_7days  -0.6136009

#------------------------------------------------------------------------------#
#               STAGE I: identify best model for each hypothesis            ----                        
#------------------------------------------------------------------------------#

# *agriculture candidate models ----

m1 <- glm(PlasmaDetection ~ PercentAg + Julian, 
          data = leye.cs,
          family = "binomial")

# log transformation does not improve fit (wt. = 0.48)
m2 <- glm(PlasmaDetection ~ NearestCropDistance_m + Julian, 
          data = leye.cs,
          family = "binomial")

m3 <- glm(PlasmaDetection ~ DominantCrop + Julian, 
          data = leye.cs,
          family = "binomial")

m4 <- glm(PlasmaDetection ~ PercentAg + NearestCropDistance_m + Julian, 
          data = leye.cs,
          family = "binomial")

model_names <- paste0("m", 1:4)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m2 3 53.71       0.00   0.33   0.33 -23.62
# m1 3 53.74       0.03   0.33   0.66 -23.63
# m3 4 54.79       1.07   0.20   0.86 -22.99
# m4 4 55.46       1.75   0.14   1.00 -23.32


# final candidate model
m.ag <- glm(PlasmaDetection ~ NearestCropDistance_m + Julian,
            family = "binomial",
            data = leye.cs)

car::vif(m.ag) # vif < 2


# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m.ag) 
plot(simulationOutput)
testDispersion(m.ag) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.ag)$Julian) #  good
plotResiduals(simulationOutput, form = model.frame(m.ag)$NearestCropDistance_m) # good

#------------------------------------------------------------------------------#
# *contaminant candidate model ----

m.pesticide <- glm(PlasmaDetection ~ EnvDetection + Julian,
                   family = "binomial",
                   data = leye.cs)

car::vif(m.pesticide) # vif < 2

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m.pesticide, 
                                      n = 2000) 
plot(simulationOutput)
testDispersion(m.pesticide) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.pesticide)$EnvDetection) # good
plotResiduals(simulationOutput, form = model.frame(m.pesticide)$Julian) # pattern but n.s.

#------------------------------------------------------------------------------#

# *precipitation candidate models ----

# log transformation not needed (wt = 0.49)
m1 <- glm(PlasmaDetection ~ DaysSinceLastPrecipitation_5mm + Julian,
          family = "binomial",
          data = leye.cs)

# log transformation improves fit (wt = 0.55)
m2 <- glm(PlasmaDetection ~ LogPrecipAmount + Julian,
          family = "binomial",
          data = leye.cs)

m3 <- glm(PlasmaDetection ~ AnnualSnowfall_in + Julian,
          family = "binomial",
          data = leye.cs)

m4 <- glm(PlasmaDetection ~ AnnualPrecipitation_in + Julian,
          family = "binomial",
          data = leye.cs)

model_names <- paste0("m", 1:4)
models <- mget(model_names)
aictab(models, modnames = model_names)

m.precip <- glm(PlasmaDetection ~ AnnualPrecipitation_in + Julian,
          family = "binomial",
          data = leye.cs)

car::vif(m.precip) # VIF ~ 3

# model validation --> some pattern but no severe issues
simulationOutput <- simulateResiduals(fittedModel = m.precip) 
plot(simulationOutput)
testDispersion(m.precip) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.precip)$Julian)  # some pattern
plotResiduals(simulationOutput, form = model.frame(m.precip)$AnnualPreciptiation_in)  # some pattern

#------------------------------------------------------------------------------#

# *drought candidate model ----
m.drought <- glm(PlasmaDetection ~ SPEI + Julian,
                 family = "binomial",
                 data = leye.cs)

car::vif(m.drought) # vif < 3

# model validation --> some pattern but no severe issues
simulationOutput <- simulateResiduals(fittedModel = m.drought) 
plot(simulationOutput)
testDispersion(m.drought) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.drought)$SPEI)  # some pattern
plotResiduals(simulationOutput, form = model.frame(m.drought)$Julian)  # some pattern

#------------------------------------------------------------------------------#
# *life history candidate models ----

m1 <- glm(PlasmaDetection ~ Age + Julian,
          family = "binomial",
          data = leye.cs)

m2 <- glm(PlasmaDetection ~ Sex + Julian,
          family = "binomial",
          data = leye.cs)



model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m1 3 56.38       0.00   0.53   0.53 -24.95
# m2 3 56.58       0.21   0.47   1.00 -25.05

m.lifehistory <- glm(PlasmaDetection ~ Age + Julian,
                     family = "binomial",
                     data = leye.cs)

car::vif(m.lifehistory) # vif < 2

# model validation --> some pattern but no severe issues
simulationOutput <- simulateResiduals(fittedModel = m.lifehistory) 
plot(simulationOutput)
testDispersion(m.lifehistory) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.lifehistory)$Age)  # good
plotResiduals(simulationOutput, form = model.frame(m.lifehistory)$Julian)  # some pattern


#------------------------------------------------------------------------------#
# *hydrology candidate models ----

m1 <- glm(PlasmaDetection ~ Permanence + Julian,
          family = "binomial",
          data = leye.cs)

m2 <- glm(PlasmaDetection ~ Porosity + Julian,
          family = "binomial",
          data = leye.cs)

# log transformation not needed (wt = 0.48)
m3 <- glm(PlasmaDetection ~ Dist_Closest_Wetland_m + Julian,
          family = "binomial",
          data=leye.cs)

m4 <- glm(PlasmaDetection ~ Permanence + Dist_Closest_Wetland_m + Julian,
          family = "binomial",
          data=leye.cs)

model_names <- paste0("m", 1:4)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m3 3 55.46       0.00   0.42   0.42 -24.49
# m1 4 56.26       0.79   0.28   0.70 -23.72
# m2 3 57.06       1.60   0.19   0.89 -25.29
# m4 5 58.22       2.76   0.11   1.00 -23.49

m.hydrology <- glm(PlasmaDetection ~ Dist_Closest_Wetland_m + Julian,
                  family = "binomial",
                  data=leye.cs)


# model validation --> some pattern but not severe
simulationOutput <- simulateResiduals(fittedModel = m.hydrology) 
plot(simulationOutput)
testDispersion(m.hydrology) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.hydrology)$Julian)  # pattern but n.s.
plotResiduals(simulationOutput, 
              form = model.frame(m.hydrology)$Dist_Closest_Wetland_m)  # pattern but n.s.

#------------------------------------------------------------------------------#
# *temporal candidate models ----

m1 <- glm(PlasmaDetection ~ time_hours + Julian,
          family = "binomial",
          data = leye.cs)

m2 <- glm(PlasmaDetection ~ Julian,
          family = "binomial",
          data = leye.cs)

model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m2 2 54.83       0.00   0.74   0.74 -25.30
# m1 3 56.92       2.09   0.26   1.00 -25.22

m.time <- glm(PlasmaDetection ~ Julian,
              family = "binomial",
              data = leye.cs)

# model validation --> some pattern but not severe
simulationOutput <- simulateResiduals(fittedModel = m.time) 
plot(simulationOutput)
testDispersion(m.time) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.time)$Julian) # some pattern

#------------------------------------------------------------------------------#
#                 STAGE II: evaluate support for hypotheses                 ----                        
#------------------------------------------------------------------------------#

models <- list(
  Time        = m.time,
  Agriculture = m.ag,
  Precip      = m.precip,
  LifeHistory = m.lifehistory,
  Hydrology   = m.hydrology,
  Pesticide   = m.pesticide,
  Drought     = m.drought
)

aictab(cand.set = models)

# Model selection based on AICc:
#   
#             K  AICc Delta_AICc AICcWt Cum.Wt     LL
# Agriculture 3 53.71       0.00   0.26   0.26 -23.62
# Pesticide   3 53.90       0.19   0.24   0.49 -23.71
# Time        2 54.83       1.12   0.15   0.64 -25.30
# Hydrology   3 55.46       1.75   0.11   0.75 -24.49
# Precip      3 55.55       1.84   0.10   0.85 -24.54
# Drought     3 56.13       2.41   0.08   0.93 -24.82
# LifeHistory 3 56.38       2.66   0.07   1.00 -24.95

confint(m.ag)
summary(m.ag) # marginal
