#---------------------------------#
# Predictors of Neonics in Plasma #
#   Analysis by Shelby McCahon    #
#      Created: 2026-01-21        #
#     Modified: 2026-01-21        #
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

# ANALYSIS NOTES:
# Event is included in every model due to strong temporal influence
# Stage I: Identify the best candidate model for each hypothesis using AICc
# -> If multiple models from stage I are within 2 deltaAICc, select the top and
# most parsimonious model

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

# combine temporary and seasonal permanence classes to increase sample size
birds <- birds %>% 
  mutate(Permanence = case_when(
    Permanence %in% c("Temporary", "Seasonal") ~ "Temporary/Seasonal",
    Permanence == "Semipermanent" ~ "Semipermanent",
    Permanence == "Permanent" ~ "Permanent"))

# log transform nearest crop distance due to right skew
birds <- birds %>% 
  mutate(CropDistance = NearestCropDistance_m + 0.0001) %>% 
  mutate(LogCropDistance = log(CropDistance))

# log transform nearest wetland distance due to right skew
birds <- birds %>% 
  mutate(WetlandDistance = Dist_Closest_Wetland_m + 0.0001) %>% 
  mutate(LogWetlandDistance = log(WetlandDistance))

# log transform precipitation amount due to right skew
birds <- birds %>% 
  mutate(LogPrecipAmount = PrecipitationAmount_7days + 0.0001) %>% 
  mutate(LogPrecipAmount = log(LogPrecipAmount))

# log transform precipitation event due to right skew
birds <- birds %>% 
  mutate(LogPrecipDays = DaysSinceLastPrecipitation_5mm + 0.0001) %>% 
  mutate(LogPrecipDays = log(LogPrecipDays))

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
#                   identify best model for each hypothesis                 ----                        
#------------------------------------------------------------------------------#

# *agriculture candidate models ----

###... LANDSCAPE-SCALE EXPOSURE HYPOTHESIS ----
# (detections are influenced by overall agricultural intensity)
# Mechanism: Birds in landsacpes with more cropland cover are more likely to
# encounter neonic residues
m1 <- glm(PlasmaDetection ~ PercentAg + Event, 
          data = birds.cs,
          family = "binomial")

###... LOCAL CROP PROXIMITY HYPOTHESIS ----
# (detections are influenced by local proximity to crops)
# Mechanism: Birds close to crop fields are more likely to be exposued

# log transformation improves fit (wt. = 0.95)

m2 <- glm(PlasmaDetection ~ LogCropDistance + Event, 
          data = birds.cs,
          family = "binomial")

###... CROP-TYPE SPECIFICITY HYPOTHESIS ----
# (detections are influenced by different crop-type exposure)
# Mechanism: Certain crop types are more heavily treated with neonics which
# affects exposure and detection
m3 <- glm(PlasmaDetection ~ DominantCrop + Event, 
          data = birds.cs,
          family = "binomial")

###... LOCAL AND LANDSCAPE EXPOSURE HYPOTHESIS ----
# (detections are influenced by local and landscape-scale cropland cover)
# Mechanism: Exposure depends on both proximity and landscale-level
# cropland cover intensity
m4 <- glm(PlasmaDetection ~ PercentAg + LogCropDistance + Event, 
          data = birds.cs,
          family = "binomial")

car::vif(m4) # vif < 2

model_names <- paste0("m", 1:4)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
  
#    K   AICc Delta_AICc AICcWt Cum.Wt     LL
# m2 5 161.02       0.00   0.61   0.61 -75.33
# m4 6 162.27       1.25   0.33   0.94 -74.88
# m3 7 166.21       5.19   0.05   0.98 -75.76
# m1 5 168.27       7.25   0.02   1.00 -78.95

# final candidate model
m.ag <- glm(PlasmaDetection ~ LogCropDistance + Event,
            family = "binomial",
            data = birds.cs)

#------------------------------------------------------------------------------#

# *contaminant candidate model ----
###... PESTICIDE EXPOSURE HYPOTHESIS ----
# (detections are influenced by the detection of pesticides in the wetland)
# Mechanism: Birds exposed to wetlands that have pesticide detections are more 
# likely to accumulate neonics in their plasma
m.pesticide <- glm(PlasmaDetection ~ EnvDetection + Event,
               family = "binomial",
               data = birds.cs)


#------------------------------------------------------------------------------#

# NOTES
# precipitation amount and days since last precipitation event are correlated
# (-0.75) so neither were included in the same model

# *precipitation candidate models ----
###...RECENT PRECIPITATION HYPOTHESIS ----
# (detections are influenced by recent precipitation events)
# Mechanism: Higher chances for agricultural runoff
# log transformation does not improve fit (wt = 0.32)
m1 <- glm(PlasmaDetection ~ DaysSinceLastPrecipitation_5mm + Event,
                   family = "binomial",
                   data = birds.cs)

###...PRECIPITATION AMOUNT RUNOFF HYPOTHESIS ----
# (detections are influenced by recent precipitation amount)
# Mechanism: More precipitation runoff from the amount of precipitation
# log transformation does not improve fit (wt = 0.46)
m2 <- glm(PlasmaDetection ~ PrecipitationAmount_7days + Event,
          family = "binomial",
          data = birds.cs)

###...SNOWFALL & SNOWMELT TRANSPORT HYPOTHESIS ----
# (detections are influenced by annual snowfall levels)
# Mechanism: Snowfall influences seasonal water availability and snowmelt 
# is known to transport neonics into wetlands
m3 <- glm(PlasmaDetection ~ AnnualSnowfall_in + Event,
          family = "binomial",
          data = birds.cs)

### TOTAL PRECIPITATION RUNOFF HYPOTHESIS ----
# (detections are influenced by annual levels of precipitation)
# Mechanism: Annual levels of precipitation influence cumulative neonic
# exposure and transport into wetlands
m4 <- glm(PlasmaDetection ~ AnnualPrecipitation_in + Event,
          family = "binomial",
          data = birds.cs)

model_names <- paste0("m", 1:4)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K   AICc Delta_AICc AICcWt Cum.Wt     LL
# m4 5 166.50       0.00   0.36   0.36 -78.07
# m1 5 166.69       0.19   0.33   0.69 -78.16
# m2 5 168.08       1.58   0.16   0.86 -78.86
# m3 5 168.35       1.85   0.14   1.00 -78.99

m.precip <- glm(PlasmaDetection ~ AnnualPrecipitation_in + Event,
            family = "binomial",
            data = birds.cs)

#------------------------------------------------------------------------------#

# *drought candidate model ----
###...LONG_TERM MOISUTRE HYPOTHESIS ----
# (detections are influenced by changes in average wetness)
# Mechanism: During wetter-than-average conditions, exposure may be higher
# due to more runoff
m.drought <- glm(PlasmaDetection ~ SPEI + Event,
                 family = "binomial",
                 data = birds.cs)

#------------------------------------------------------------------------------#

# *life history candidate models ----

# GET RID OF DOWITCHERS FOR THIS SET
m1 <- glm(PlasmaDetection ~ Species + Event,
          family = "binomial",
          data = birds.cs)

m2 <- glm(PlasmaDetection ~ Sex + Event,
          family = "binomial",
          data = birds.cs)


m3 <- glm(PlasmaDetection ~ MigStatus + Event,
          family = "binomial",
          data = birds.cs)


m4 <- glm(PlasmaDetection ~ Species + Sex + Event,
          family = "binomial",
          data = birds.cs)

m5 <- glm(PlasmaDetection ~ MigStatus + Sex + Event,
          family = "binomial",
          data = birds.cs)


#------------------------------------------------------------------------------#

# *hydrology candidate models ----

m1 <- glm(PlasmaDetection ~ Permanence + Event,
          family = "binomial",
          data = birds.cs)

m2 <- glm(PlasmaDetection ~ Porosity + Event,
          family = "binomial",
          data = birds.cs)

m3 <- glm(PlasmaDetection ~ Permanence + Event,
          family = "binomial",
          data = birds.cs)

# CHECK FOR LOG TRANSFORMATION
m4 <- glm(PlasmaDetection ~ Dist_Closest_Wetland_m + Event,
          family = "binomial",
          data = birds.cs)

