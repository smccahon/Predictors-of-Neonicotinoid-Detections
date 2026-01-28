#---------------------------------#
# Predictors of Neonics in Plasma #
#   Analysis by Shelby McCahon    #
#      Created: 2026-01-21        #
#     Modified: 2026-01-28        #
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
# STAGE I: Identify the best candidate model for each hypothesis using AICc

###  If multiple models from stage I are supported (within 2 deltaAICc), we 
### selected the most supported and parsimonious model. 


# STAGE II:
# Model selection of the top hypotheses

# I identified correlations from an earlier analysis. Season & SPEI (r = 0.90)
# and Precipitation Amount & Days Since Last Precipitation Event (r = -0.75) 

# variables with correlation > 0.80 not considered in the same model

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
  mutate(LogPrecipAmount = PrecipitationAmount_7days_cm + 0.0001) %>% 
  mutate(LogPrecipAmount = log(LogPrecipAmount))

# log transform precipitation event due to right skew
birds <- birds %>% 
  mutate(LogPrecipDays = DaysSinceLastPrecipitation_5mm + 0.0001) %>% 
  mutate(LogPrecipDays = log(LogPrecipDays))

# filter out long-billed dowitchers to test for the effect of species
# birds <- birds %>% 
#   filter(Species != "Longbilled Dowitcher")

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
#               STAGE I: identify best model for each hypothesis            ----                        
#------------------------------------------------------------------------------#

# *agriculture candidate models ----

###... LANDSCAPE-SCALE EXPOSURE HYPOTHESIS ----
# (detections are influenced by overall agricultural intensity)
# Mechanism: Birds in landscapes with more cropland cover are more likely to
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
# Mechanism: Exposure depends on both proximity and landscape-level
# cropland cover intensity
m4 <- glm(PlasmaDetection ~ PercentAg + LogCropDistance + Event, 
          data = birds.cs,
          family = "binomial")

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

car::vif(m.ag) # vif < 2

# model validation --> good with more simulations
simulationOutput <- simulateResiduals(fittedModel = m.ag, n = 2000) 
plot(simulationOutput)
testDispersion(m.ag) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.ag)$LogCropDistance) # good
plotResiduals(simulationOutput, form = model.frame(m.ag)$Event) # good

#------------------------------------------------------------------------------#

# *contaminant candidate model ----
###... PESTICIDE EXPOSURE HYPOTHESIS ----
# (detections are influenced by the detection of pesticides in the wetland)
# Mechanism: Birds exposed to wetlands that have pesticide detections are more 
# likely to accumulate neonics in their plasma
m.pesticide <- glm(PlasmaDetection ~ EnvDetection + Event,
               family = "binomial",
               data = birds.cs)

car::vif(m.pesticide) # vif < 3

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m.pesticide, n = 2000) 
plot(simulationOutput)
testDispersion(m.pesticide) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.pesticide)$EnvDetection) # good
plotResiduals(simulationOutput, form = model.frame(m.pesticide)$Event) # good

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
m2 <- glm(PlasmaDetection ~ PrecipitationAmount_7days_cm + Event,
          family = "binomial",
          data = birds.cs)

###...SNOWFALL & SNOWMELT TRANSPORT HYPOTHESIS ----
# (detections are influenced by annual snowfall levels)
# Mechanism: Snowfall influences seasonal water availability and snowmelt 
# is known to transport neonics into wetlands
m3 <- glm(PlasmaDetection ~ AnnualSnowfall_cm + Event,
          family = "binomial",
          data = birds.cs)

###... TOTAL PRECIPITATION RUNOFF HYPOTHESIS ----
# (detections are influenced by annual levels of precipitation)
# Mechanism: Annual levels of precipitation influence cumulative neonic
# exposure and transport into wetlands
m4 <- glm(PlasmaDetection ~ AnnualPrecipitation_cm + Event,
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

# does including both m1 and m4 improve fit? no
# m1 <- glm(PlasmaDetection ~ AnnualPrecipitation_cm + 
#             DaysSinceLastPrecipitation_5mm + Event,
#           family = "binomial",
#           data = birds.cs)
# 
# m2 <- glm(PlasmaDetection ~ AnnualPrecipitation_cm + Event,
#           family = "binomial",
#           data = birds.cs)
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)

m.precip <- glm(PlasmaDetection ~ AnnualPrecipitation_cm + Event,
            family = "binomial",
            data = birds.cs)

car::vif(m.precip) # vif < 2

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m.precip, n = 1000) 
plot(simulationOutput)
testDispersion(m.precip) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.precip)$EnvDetection) # good
plotResiduals(simulationOutput, form = model.frame(m.precip)$AnnualPrecipitation_cm) # good

#------------------------------------------------------------------------------#

# *drought candidate model ----
###...LONG_TERM MOISUTRE HYPOTHESIS ----
# (detections are influenced by changes in average wetness)
# Mechanism: During wetter-than-average conditions, exposure may be higher
# due to more runoff
m.drought <- glm(PlasmaDetection ~ SPEI + Event,
                 family = "binomial",
                 data = birds.cs)

car::vif(m.drought) # vif < 3

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m.drought, n = 1000) 
plot(simulationOutput)
testDispersion(m.drought) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.drought)$EnvDetection) # good
plotResiduals(simulationOutput, form = model.frame(m.drought)$SPEI) # good

#------------------------------------------------------------------------------#

# *life history candidate models ----
# results do not include LBDO due to no variation in detection for that species

###...SPECIES HYPOTHESIS
# (detections will differ by species)
# Mechanism: Shorebirds have different diets and forage at different depths and
# times which could affect exposure
m1 <- glm(PlasmaDetection ~ Species + Event,
          family = "binomial",
          data = birds.cs)

###...SEX HYPOTHESIS
# (detections will differ by sex)
# Mechanism: In some migratory shorebird species, one sex migrates first (e.g.,
# Lesser Yellowlegs). Some sexes might be differentially exposed to neonics. 
m2 <- glm(PlasmaDetection ~ Sex + Event,
          family = "binomial",
          data = birds.cs)

###... MIGRANT STATUS HYPOTHESIS
# (detections will differ between migrants and residents)
# Mechanism: Migrant species may have a higher probability of neonic exposure
# due to the limited amount of time they have to forage
m3 <- glm(PlasmaDetection ~ MigStatus + Event,
          family = "binomial",
          data = birds.cs)

###... SPECIES AND SEX HYPOTHESIS
# (detections will vary between species AND sex)
m4 <- glm(PlasmaDetection ~ Species + Sex + Event,
          family = "binomial",
          data = birds.cs)


###... MIGRANT STATUS AND SEX
m5 <- glm(PlasmaDetection ~ MigStatus + Sex + Event,
          family = "binomial",
          data = birds.cs)

# Model selection based on AICc:
#   
#     K   AICc Delta_AICc AICcWt Cum.Wt     LL
# m2  5 161.22       0.00   0.58   0.58 -75.42
# m3  5 163.16       1.94   0.22   0.80 -76.39
# m5  6 163.36       2.13   0.20   1.00 -75.41
# m1 14 178.17      16.94   0.00   1.00 -73.68
# m4 15 179.67      18.44   0.00   1.00 -73.22

model_names <- paste0("m", 1:5)
models <- mget(model_names)
aictab(models, modnames = model_names)

m.lifehistory <- glm(PlasmaDetection ~ Sex + Event,
                     family = "binomial",
                     data = birds.cs)

car::vif(m.lifehistory) # vif < 2

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m.lifehistory, n = 1000) 
plot(simulationOutput)
testDispersion(m.lifehistory) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.lifehistory)$Event) # good
plotResiduals(simulationOutput, form = model.frame(m.lifehistory)$Sex) # good
#------------------------------------------------------------------------------#

# *hydrology candidate models ----

###... PERMANENCE HYPOTHESIS ----
# (detections will be influenced by how long wetlands retain water)
# Mechanism: Birds captured in permanence wetlands may have a higher
# probability of accumulating neonics because permanence wetlands may
# accumulate more pesticides over time
m1 <- glm(PlasmaDetection ~ Permanence + Event,
          family = "binomial",
          data = birds.cs)

###... SHORELINE POROSITY HYPTHESIS ----
# (detections will be higher in more porous wetlands)
# Mechanism: Porosity affects water infiltration and runoff
m2 <- glm(PlasmaDetection ~ Porosity + Event,
          family = "binomial",
          data = birds.cs)

###... PROXIMITY TO WETLAND HYPOTHESIS ----
# (detections will be influenced by wetland connectivity)
# Mechanism: Birds captured in isolated wetlands will have less detections
# compared to birds captured in connected wetlands because they have more
# exposure opportunities prior to capture
# log transformation does NOT improve model fit (wt = 0.41)
m3 <- glm(PlasmaDetection ~ Dist_Closest_Wetland_m + Event,
          family = "binomial",
          data = birds.cs)

###... WETLAND HYDROLOGY HYPOTHESIS ----
# (detections will be influenced both by wetland connectivity and permanence)
# Mechanism: Birds captured in larger wetlands that are more connected to
# others will have a higher probability of exposure
m4 <- glm(PlasmaDetection ~ Permanence + Dist_Closest_Wetland_m + Event,
          family = "binomial",
          data= birds.cs)

model_names <- paste0("m", 1:4)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K   AICc Delta_AICc AICcWt Cum.Wt     LL
# m3 5 161.02       0.00   0.57   0.57 -75.33
# m4 7 161.71       0.69   0.40   0.97 -73.51
# m2 5 168.06       7.04   0.02   0.99 -78.85
# m1 6 168.75       7.73   0.01   1.00 -78.12

m.hydrology <- glm(PlasmaDetection ~ Dist_Closest_Wetland_m + Event ,
                   family = "binomial",
                   data = birds.cs)

car::vif(m.hydrology) # vif < 2

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m.hydrology, n = 1000) 
plot(simulationOutput)
testDispersion(m.hydrology) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m.hydrology)$Event) # good
plotResiduals(simulationOutput, form = model.frame(m.hydrology)$Dist_Closest_Wetland_m) # good

#------------------------------------------------------------------------------#

# *temporal candidate models ----

###... FORAGING TIME HYPOTHESIS ----
# (detections increase during foraging periods)
# Mechanism: Birds are more likely to be exposed during periods of 
# foraging 
# No model violation issues with time (given cyclical nature)
m1 <- glm(PlasmaDetection ~ time_hours + Event,
          family = "binomial",
          data = birds.cs)

###... SAMPLING EVENT HYPOTHESIS/INFORMED NULL ----
# (detections are only influenced by the season and year we captured birds in)
m2 <- glm(PlasmaDetection ~ Event,
          family = "binomial",
          data = birds.cs)

model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)

m.time <- glm(PlasmaDetection ~ Event,
              family = "binomial",
              data = birds.cs)

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
# 
# Model selection based on AICc:
#   
#             K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Agriculture 5 161.02       0.00   0.44   0.44 -75.33
# Hydrology   5 161.02       0.00   0.44   0.89 -75.33
# Time        4 166.31       5.30   0.03   0.92 -79.04
# LifeHistory 5 166.32       5.30   0.03   0.95 -77.98
# Precip      5 166.50       5.48   0.03   0.98 -78.07
# Pesticide   5 168.42       7.40   0.01   0.99 -79.03
# Drought     5 168.43       7.41   0.01   1.00 -79.03

#------------------------------------------------------------------------------#
#                   STAGE III: assessing combined effects                   ----                        
#------------------------------------------------------------------------------#

m.combined <- glm(PlasmaDetection ~ LogCropDistance + Event + 
              Dist_Closest_Wetland_m,
              family = "binomial",
              data = birds.cs)

models <- list(
  Agriculture = m.ag,
  Hydrology   = m.hydrology,
  Combined    = m.combined)

aictab(cand.set = models)

# Model selection based on AICc:
#   
#             K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Combined    6 160.54       0.00   0.39   0.39 -74.02
# Agriculture 5 161.02       0.47   0.31   0.69 -75.33
# Hydrology   5 161.02       0.48   0.31   1.00 -75.33 

summary(m.combined)

car::vif(m.combined) # vif < 2

#------------------------------------------------------------------------------#
#                               model average                               ----                        
#------------------------------------------------------------------------------#

model_avg <- model.avg(m.ag, m.hydrology, m.time, m.lifehistory)

summary(model_avg)
confint(model_avg)

# Model-averaged coefficients:  
#   (full average) 
#                        Estimate Std. Error Adjusted SE z value Pr(>|z|)    
# (Intercept)            -0.62018    0.45122     0.45436   1.365 0.172270    
# LogCropDistance        -0.27650    0.33337     0.33389   0.828 0.407618    
# EventFall 2023         -0.96080    0.54156     0.54547   1.761 0.078167 .  
# EventSpring 2022        4.21093    1.12903     1.13707   3.703 0.000213 ***
# EventSpring 2023       -0.64941    0.69138     0.69583   0.933 0.350671    
# Dist_Closest_Wetland_m  0.29093    0.34835     0.34887   0.834 0.404324    
# SexM                   -0.01941    0.12819     0.12849   0.151 0.879908   

#                             2.5 %     97.5 %
# (Intercept)            -1.5107032  0.2703519
# LogCropDistance        -1.0388668 -0.1446794
# EventFall 2023         -2.0298980  0.1082976
# EventSpring 2022        1.9823231  6.4395378
# EventSpring 2023       -2.0132199  0.7143937
# Dist_Closest_Wetland_m  0.1697762  1.0769380
# SexM                   -1.3861720  0.2084972

#------------------------------------------------------------------------------#
#                               graph results                               ----                        
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
           "Sex (Female)"), 
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
                 "Dist_Closest_Wetland_m",
                 "LogCropDistance",
                 "SexM",
                 "Sex (Female)")),
  labels = rev(c("Sampling Event (Spring 2022)",
                 "Sampling Event (Spring 2023)",
                 "Sampling Event (Fall 2023)",
                 "Sampling Event (Fall 2021)",
                 "Distance to Nearest Wetland",
                 "Log(Distance to Nearest Crop)",
                 "Sex (Male)",
                 "Sex (Female)")))

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

# save figure
ggsave(filename = "Model-averaged Plasma Neonic Results_2026-01-22.tiff",
       path = "C:/Users/shelb/OneDrive - University of Idaho/Masters Project/Analysis/Predictors of Neonics/figures", 
       width = 13, height = 8, units = "in", device='tiff', dpi=300)

#------------------------------------------------------------------------------#
#                             PLOT TOP MODEL(s)                             ----                        
#------------------------------------------------------------------------------#

### ...agriculture model ----
# convert to factor to add points to figure
birds$PlasmaDetectionNum <- ifelse(birds$PlasmaDetection  == "Y", 
                                   1, 0)

# add colors
cols <- c(
  "Fall 2021"   = "#875907",  
  "Fall 2023"   = "#FADBA4", 
  "Spring 2023" = "#98CCF0", 
  "Spring 2022" = "#0D3957"  
)


m <- glm(PlasmaDetection ~ LogCropDistance + Event,
         family = "binomial",
         data = birds)

d <- expand.grid(
  LogCropDistance = seq(
    min(birds$LogCropDistance),
    max(birds$LogCropDistance),
    length = 1000),
  Event = unique(birds$Event))

# not necessary, but nice to make sure it's using the right type
d$yhat <- predict(m, newdata = d, type = "response")

d <- cbind(d, trtools::glmint(m, newdata = d))

head(d)

ggplot(d, aes(x = LogCropDistance, y = yhat, 
              color = Event, fill = Event)) +
  geom_line(linewidth = 1.2) + 
  geom_ribbon(aes(ymin = low, ymax = upp), 
              alpha = 0.4, color = NA) +
  facet_wrap(~Event) +
  theme_classic() +
  labs(x ="Log(Distance to Nearest Crop)", 
       y = "Probability of Neonicotinoid Detection\nin Shorebird Plasma (%)") +
  theme(strip.text = element_text(face = "bold", color = "black",
                                  hjust = 0.5, size = 15),
        strip.background = element_rect(fill = "lightgrey", 
                                        linetype = "solid",
                                        color = "black", linewidth = 0.5),
        panel.spacing = unit(20, 'points'),
        panel.border = element_rect(fill = "transparent",
                                    color = "black", linewidth = 0.5),
        axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  # geom_point(data = birds, aes(x = LogCropDistance,
  #                                y = PlasmaDetection), size = 2) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(expand = c(0,0), 
                     limits = c(-10,8.2)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)

# save high res figure
ggsave(filename = "ALL SPECIES AGRICULTURE Plasma Neonic Detection Results_2026-01-28.tiff",
       path = "C:/Users/shelb/OneDrive - University of Idaho/Masters Project/Analysis/Predictors of Neonics/figures", 
       width = 13, height = 8, units = "in", device='tiff', dpi=300)


#------------------------------------------------------------------------------#


### ...hydrology model ----
# convert to factor to add points to figure
birds$PlasmaDetectionNum <- ifelse(birds$PlasmaDetection  == "Y", 
                                   1, 0)

# add colors
cols <- c(
  "Fall 2021"   = "#875907",  
  "Fall 2023"   = "#FADBA4", 
  "Spring 2023" = "#98CCF0", 
  "Spring 2022" = "#0D3957"  
)


m <- glm(PlasmaDetection ~ Dist_Closest_Wetland_m + Event,
         family = "binomial",
         data = birds)

d <- expand.grid(
  Dist_Closest_Wetland_m = seq(
    min(birds$Dist_Closest_Wetland_m),
    max(birds$Dist_Closest_Wetland_m),
    length = 1000),
  Event = unique(birds$Event))

# not necessary, but nice to make sure it's using the right type
d$yhat <- predict(m, newdata = d, type = "response")

d <- cbind(d, trtools::glmint(m, newdata = d))

head(d)

ggplot(d, aes(x = Dist_Closest_Wetland_m, y = yhat, 
              color = Event, fill = Event)) +
  geom_line(linewidth = 1.2) + 
  geom_ribbon(aes(ymin = low, ymax = upp), 
              alpha = 0.4, color = NA) +
  facet_wrap(~Event) +
  theme_classic() +
  labs(x ="Distance to Nearest Wetland", 
       y = "Probability of Neonicotinoid Detection\nin Shorebird Plasma (%)") +
  theme(strip.text = element_text(face = "bold", color = "black",
                                  hjust = 0.5, size = 15),
        strip.background = element_rect(fill = "lightgrey", 
                                        linetype = "solid",
                                        color = "black", linewidth = 0.5),
        panel.spacing = unit(20, 'points'),
        panel.border = element_rect(fill = "transparent",
                                    color = "black", linewidth = 0.5),
        axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  # geom_point(data = birds, aes(x = LogCropDistance,
  #                                y = PlasmaDetection), size = 2) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(expand = c(0,0), 
                     limits = c(0,410)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)

# save high res figure
ggsave(filename = "ALL SPECIES HYDROLOGY Plasma Neonic Detection Results_2026-01-28.tiff",
       path = "C:/Users/shelb/OneDrive - University of Idaho/Masters Project/Analysis/Predictors of Neonics/figures", 
       width = 13, height = 8, units = "in", device='tiff', dpi=300)
