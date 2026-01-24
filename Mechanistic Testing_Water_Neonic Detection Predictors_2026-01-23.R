#-----------------------------------#
#  Predictors of Neonics in Water   #
#     Analysis by Shelby McCahon    #
#        Created: 2026-01-23        #
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

water <- read.csv("cleaned_data/wetland_data_cleaned_2025-09-30.csv")

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
#  water <- water %>%
#    mutate(Buffered = ifelse(Buffered == "Y", 1, 0))
# # 
#  water <- water %>%
#    mutate(Season = ifelse(Season == "Spring", 1, 0))
#  
# combine temporary and seasonal permanence classes to increase sample size
water <- water %>% 
 mutate(Permanence = case_when(
   Permanence %in% c("Temporary", "Seasonal") ~ "Temporary/Seasonal",
   Permanence == "Semipermanent" ~ "Semipermanent",
   Permanence == "Permanent" ~ "Permanent"))

# convert characters to factors
water <- water %>% 
  mutate_if(is.character, as.factor)

# filter data (n = 97)
water <- water %>%
  filter(!is.na(WaterNeonicDetection))

# remove wetlands with corn to see if dominant crop has effect (n = 94)
# it doesn't, so include the three wetlands back in
# water <- water %>%
#   filter(DominantCrop != "Corn")

# log transform nearest crop distance due to skew
water <- water %>% 
  mutate(LogCropDistance = NearestCropDistance_m + 0.0001) %>% 
  mutate(LogCropDistance = log(LogCropDistance))

# log transform nearest wetland distance due to skew
water <- water %>% 
  mutate(LogWetlandDistance = Dist_Closest_Wetland_m) %>% 
  mutate(LogWetlandDistance = log(LogWetlandDistance))

# log transform precipitation amount due to skew
water <- water %>% 
  mutate(LogPrecipAmount = PrecipitationAmount_7days + 0.0001) %>% 
  mutate(LogPrecipAmount = log(LogPrecipAmount))

# log transform precipitation events due to skew
water <- water %>% 
  mutate(LogPrecipDays = DaysSinceLastPrecipitation_5mm + 0.0001) %>% 
  mutate(LogPrecipDays = log(LogPrecipDays))

# standardize data
water.cs <- water %>%
  mutate(across(where(is.numeric), scale))

# set reference levels
water.cs$DominantCrop <- relevel(water.cs$DominantCrop,
                                  ref = "Grassland")

water.cs$Permanence <- relevel(water.cs$Permanence,
                                ref = "Temporary/Seasonal")

#------------------------------------------------------------------------------#
#                         identify correlations                             ----                        
#------------------------------------------------------------------------------# 

 #                               Var1                      Var2 Correlation
 # 27                          Julian                    Season  -0.9740772
 # 132                         Season                      SPEI   0.8481201
 # 131                         Julian                      SPEI  -0.7800471
 # 111                       Buffered     NearestCropDistance_m   0.7386995
 # 85                        Buffered                 PercentAg  -0.7344918
 # 83           NearestCropDistance_m                 PercentAg  -0.7133819
 # 443                         Julian                 WaterTemp   0.7117079
 # 444                         Season                 WaterTemp  -0.6756220
 # 638 DaysSinceLastPrecipitation_5mm           LogPrecipAmount  -0.6693172
 # 42                       Diversity                    Season  -0.6197850
 # 337                LogPrecipAmount PrecipitationAmount_7days   0.6137498
 # 101                LogCropDistance                 PercentAg  -0.6007487

#------------------------------------------------------------------------------#
#                   identify best model for each hypothesis                 ----                        
#------------------------------------------------------------------------------#

# *agriculture candidate models ----

m1 <- glm(WaterNeonicDetection ~ PercentAg + Event, 
          data = water.cs,
          family = "binomial")

# log transformation improves fit (wt. = 0.97)

m2 <- glm(WaterNeonicDetection ~ LogCropDistance + Event, 
          data = water.cs,
          family = "binomial")

m3 <- glm(WaterNeonicDetection ~ DominantCrop + Event, 
          data = water.cs,
          family = "binomial")

m4 <- glm(WaterNeonicDetection ~ Buffered + Event, 
          data = water.cs,
          family = "binomial")

m5 <- glm(WaterNeonicDetection ~ PercentAg + LogCropDistance +
            Event, 
          data = water.cs,
          family = "binomial")

m6 <- glm(WaterNeonicDetection ~ Buffered + PercentAg + Event, 
          data = water.cs,
          family = "binomial")

# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m2 5 85.56       0.00   0.65   0.65 -37.45
# m5 6 87.82       2.26   0.21   0.86 -37.44
# m1 5 89.86       4.30   0.08   0.93 -39.60
# m4 5 91.66       6.10   0.03   0.96 -40.50
# m6 6 92.12       6.56   0.02   0.99 -39.59
# m3 8 93.74       8.18   0.01   1.00 -38.05


model_names <- paste0("m", 1:6)
models <- mget(model_names)
aictab(models, modnames = model_names)

m.ag <- glm(WaterNeonicDetection ~ LogCropDistance + Event, 
            data = water.cs,
            family = "binomial")


# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m.ag, 
                                      n = 2000) 
plot(simulationOutput)
testDispersion(m.ag) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, 
              form = model.frame(m.ag)$LogCropDistance) # good
plotResiduals(simulationOutput, 
              form = model.frame(m.ag)$Event) # good

#------------------------------------------------------------------------------#

# *hydrology candidate models ----

m1 <- glm(WaterNeonicDetection ~ Permanence + Event,
          family = "binomial",
          data = water.cs)

m2 <- glm(WaterNeonicDetection ~ Porosity + Event,
          family = "binomial",
          data = water.cs)

# log transformation not supported (wt = 0.54)
m3 <- glm(WaterNeonicDetection ~ Dist_Closest_Wetland_m + Event,
          family = "binomial",
          data = water.cs)

m4 <- glm(WaterNeonicDetection ~ Permanence +
            LogWetlandDistance + Event,
          family = "binomial",
          data= water.cs)

model_names <- paste0("m", 1:4)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m1 6 89.43       0.00   0.52   0.52 -38.25
# m2 5 91.46       2.04   0.19   0.71 -40.40
# m4 7 91.71       2.28   0.17   0.88 -38.23
# m3 5 92.33       2.90   0.12   1.00 -40.84

m.hydrology <- glm(WaterNeonicDetection ~ Permanence + Event,
                   family = "binomial",
                   data = water.cs)

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

# *precipitation candidate models ----

# log transformation not supported (wt = 0.39)
m1 <- glm(WaterNeonicDetection ~ DaysSinceLastPrecipitation_5mm + Event,
          family = "binomial",
          data = water.cs)

# log transformation was not supported (wt = 0.46)
m2 <- glm(WaterNeonicDetection ~ PrecipitationAmount_7days + Event,
          family = "binomial",
          data = water.cs)


m3 <- glm(WaterNeonicDetection ~ AnnualSnowfall_in + Event,
          family = "binomial",
          data = water.cs)

m4 <- glm(WaterNeonicDetection ~ AnnualPrecipitation_in + Event,
          family = "binomial",
          data = water.cs)

model_names <- paste0("m", 1:4)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m4 5 90.92       0.00   0.32   0.32 -40.13
# m3 5 90.93       0.01   0.32   0.64 -40.13
# m1 5 91.97       1.05   0.19   0.83 -40.65
# m2 5 92.24       1.32   0.17   1.00 -40.79

m.precip <- glm(WaterNeonicDetection ~ AnnualPrecipitation_in + Event,
                family = "binomial",
                data = water.cs)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m.precip) 
plot(simulationOutput)
testDispersion(m.precip) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = 
                model.frame(m.precip)$AnnualPrecipitation_in) # good

plotResiduals(simulationOutput, form = 
                model.frame(m.precip)$Event) # good
#------------------------------------------------------------------------------#

# *temporal candidate model ----

m.time <- glm(WaterNeonicDetection ~ Event,
              data = water.cs,
              family = "binomial",
              na.action = na.fail)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m.time) 
plot(simulationOutput)
testDispersion(m.time) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.time)$Event) #  good

#------------------------------------------------------------------------------#

# null model

m.null <- glm(WaterNeonicDetection ~ 1,
              family = "binomial",
              data = water.cs)

#------------------------------------------------------------------------------#

# *drought candidate model ----
m.drought <- glm(WaterNeonicDetection ~ SPEI + Event,
              family = "binomial",
              data = water.cs)


#------------------------------------------------------------------------------#
#                 STAGE II: evaluate support for hypotheses                 ----                        
#------------------------------------------------------------------------------#

models <- list(
  Time        = m.time,
  Agriculture = m.ag,
  Precip      = m.precip,
  Hydrology   = m.hydrology,
  Null        = m.null,
  Drought     = m.drought
)

aictab(cand.set = models)

# Model selection based on AICc:
#   
#             K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Agriculture 5  85.56       0.00   0.74   0.74 -37.45
# Hydrology   6  89.43       3.87   0.11   0.84 -38.25
# Time        4  90.66       5.10   0.06   0.90 -41.11
# Precip      5  90.92       5.36   0.05   0.95 -40.13
# Drought     5  91.03       5.48   0.05   1.00 -40.19
# Null        1 110.58      25.03   0.00   1.00 -54.27


#------------------------------------------------------------------------------#
#                             PLOT TOP MODEL(s)                             ----                        
#------------------------------------------------------------------------------#


