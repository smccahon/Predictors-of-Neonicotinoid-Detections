#--------------------------------------#
# Predictors of Neonics in LEYE Plasma #
#       Analysis by Shelby McCahon     #
#        Created: 2025-11-13           #
#        Modified: 2025-11-13          #
#--------------------------------------#

# load packages
library(tidyverse)
library(AICcmodavg)
library(MASS)
library(ggplot2)
library(dplyr)
library(nlme)
library(lme4)
library(glmmTMB)
library(DHARMa)

# read in data
leye <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")

# subset to lesser yellowlegs (n = 54)
leye <- leye %>% 
  filter(Species == "Lesser Yellowlegs")

# manipulate data
leye <- leye %>% 
  mutate(Permanence = case_when(
    Permanence %in% c("Temporary", "Seasonal") ~ "Temporary/Seasonal",
    Permanence == "Semipermanent" ~ "Semipermanent",
    Permanence == "Permanent" ~ "Permanent",
    TRUE ~ NA_character_), 
  PlasmaDetection = case_when(
    PlasmaDetection == "Y" ~ 1,
    PlasmaDetection == "N" ~ 0,
    TRUE ~ NA_real_),
  EnvDetection = factor(EnvDetection,
                        levels = c("Y", "N")),
  Sex = factor(Sex,
               levels = c("M", "F"),
               labels = c("Male", "Female")),
  Age = factor(Age,
               levels = c("Adult", "Juvenile")),
  DominantCrop = factor(DominantCrop,
                        levels = c("Grassland",
                                   "Soybean",
                                   "Wheat")),
  Event = factor(Event,
                 levels = c("Fall 2021",
                            "Spring 2022",
                            "Fall 2023"),
                 labels = c("2021 Fall",
                            "2022 Spring",
                            "2023 Fall")))
  
leye$Site <- as.factor(leye$Site)

# filter data
# leye <- leye %>%
#   filter(!is.na(PlasmaDetection)) %>% 
#   group_by(Site) %>% 
#   filter(n() >= 3) %>% # removes 10 leye...
#   ungroup()
  
# standardize data except for response
leye.cs <- leye %>%
  mutate(across(where(is.numeric) & !matches("PlasmaDetection"), scale))


#------------------------------------------------------------------------------#
#                          identify correlations                            ----
#------------------------------------------------------------------------------#

 cor_results <- leye %>%
   # select only numeric columns
   dplyr::select(where(is.numeric)) %>% 
   
   # compute correlation matrix
   cor(use = "pairwise.complete.obs") %>%
   
   # turn matrix into tidy tibble
   as.data.frame() %>%
   rownames_to_column(var = "Var1") %>%
   pivot_longer(-Var1, names_to = "Var2", values_to = "Correlation") %>%
   
   # remove self-correlations and duplicates
   filter(Var1 < Var2) %>%
   
   # keep only strong correlations
   filter(abs(Correlation) >= 0.6) %>%
   
   # arrange by strength of correlation
   arrange(desc(abs(Correlation))) %>%
   
   print()

### ...identify correlations among continuous variables ----
#------------------------------------------------------------------------------#
# Var1                           Var2                      Correlation
#   1 seconds_since_midnight         time_hours                      1    
# 2 Biomass                        NearestCropDistance_m           0.993
# 3 PecSizeBest                    Standardized.Pec.NoEvent        0.989
# 4 BCI.NoEvent                    Mass                            0.963
# 5 Beta                           FatteningIndex                 -0.932
# 6 DaysSinceLastPrecipitation_5mm PrecipitationAmount_7days      -0.849
# 7 FatteningIndex                 Tri                             0.777
# 8 BCI.NoEvent                    Standardized.Pec.NoEvent        0.716
# 9 Mass                           PecSizeBest                     0.712
# 10 BCI.NoEvent                    PecSizeBest                     0.708
# 11 Diversity                      Julian                         -0.706
# 12 AnnualSnowfall_in              SPEI                            0.701
# 13 Mass                           Standardized.Pec.NoEvent        0.697
# 14 Biomass                        PercentAg                      -0.692
# 15 DaysSinceLastPrecipitation_5mm Dist_Closest_Wetland_m          0.667
# 16 NearestCropDistance_m          PercentAg                      -0.662
# 17 Diversity                      SPEI                           -0.658
# 18 Diversity                      NearestCropDistance_m           0.620
# 19 Dist_Closest_Wetland_m         PrecipitationAmount_7days      -0.614
# 20 Julian                         PlasmaDetection                -0.605
#------------------------------------------------------------------------------#

### ...identify correlations among two categorical variables                ----

#------------------------------------------------------------------------------#
# season and permanence
# event and everything...
#------------------------------------------------------------------------------#
kruskal.test(AnnualSnowfall_in ~ Event, data = leye) # corr
kruskal.test(SPEI ~ Event, data = leye) # corr
kruskal.test(PrecipitationAmount_7days ~ Event, data = leye) # corr
kruskal.test(Julian ~ Event, data = leye) # corr
kruskal.test(Percent_Exposed_Shoreline ~ Event, data = leye) # NOT corr
kruskal.test(PercentAg ~ Event, data = leye) # corr
kruskal.test(NearestCropDistance_m ~ Event, data = leye) # corr

chisq.test(table(leye$Event, leye$Permanence)) # corr
chisq.test(table(leye$Event, leye$MigStatus))  # corr
chisq.test(table(leye$Event, leye$Sex)) # NOT corr
chisq.test(table(leye$Event, leye$EnvDetection)) # NOT corr
chisq.test(table(leye$Event, leye$Age)) # corr

chisq.test(table(leye$Season, leye$Permanence)) # corr
chisq.test(table(leye$Event, leye$Permanence))  # corr
chisq.test(table(leye$Event, leye$DominantCrop)) # NOT corr

m <- glm(PlasmaDetection ~ Event + DominantCrop, data = leye,
         family = binomial(link = "logit"))

### ...identify correlations among continuous and categorical variables     ----

#------------------------------------------------------------------------------#
#  season and julian
#  season and event
#  season and SPEI
#  julian and event
#  % cropland and dominant crop                                        
#  % cropland and nearest crop distance
#  drought and event
#------------------------------------------------------------------------------#

# Mann-Whitney U test (2 groups, t-test is parametric)
t.test(Julian ~ Season, data = leye)
wilcox.test(Julian ~ Season, data = leye)
wilcox.test(SPEI ~ Season, data = leye)

# Kruskal-Wallis test (3+ groups, ANOVA is parametric)
summary(aov(AnnualSnowfall_in ~ Event, data = leye))
kruskal.test(Julian ~ Event, data = leye)
kruskal.test(SPEI ~ Season, data = leye)
kruskal.test(AnnualSnowfall_in ~ Event, data = leye)

kruskal.test(PercentAg ~ DominantCrop, data = leye) 
kruskal.test(PercentAg ~ NearestCropDistance_m, data = leye)

#------------------------------------------------------------------------------#
#               identify ideal random effect structure                      ----
#------------------------------------------------------------------------------#

m1 <- glmmTMB(PlasmaDetection ~ (1|Event/Site), 
              data = leye.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ (1|Event), 
              data = leye.cs,
              family = binomial(link = "logit"))

m3 <- glmmTMB(PlasmaDetection ~ (1|Site),
              data = leye.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:3)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K    AICc Delta_AICc AICcWt Cum.Wt     LL
# m2 2 174.89       0.00   0.53   0.53 -85.41
# m1 3 175.17       0.28   0.46   1.00 -84.51
# m3 2 186.67      11.77   0.00   1.00 -91.30

#------------------------------------------------------------------------------#
#    model selection for corr. variables  (site & event as a random effect) ----
#------------------------------------------------------------------------------#

# agriculture variable --> Nearest Crop Distance is best (wt. = 0.45)
m1 <- glmmTMB(PlasmaDetection ~ PercentAg + (1|Event/Site), 
              data = leye.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ DominantCrop + (1|Event/Site), 
              data = leye.cs,
              family = binomial(link = "logit"))

m3 <- glmmTMB(PlasmaDetection ~ NearestCropDistance_m + (1|Event/Site),
              data = leye.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:3)
models <- mget(model_names)
aictab(models, modnames = model_names)

# agriculture variable --> Nearest Crop Distance is best (wt. = 0.68)
m1 <- glmmTMB(PlasmaDetection ~ PercentAg + (1|Site), 
              data = leye.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ DominantCrop + (1|Site), 
              data = leye.cs,
              family = binomial(link = "logit"))

m3 <- glmmTMB(PlasmaDetection ~ NearestCropDistance_m + (1|Site),
              data = leye.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:3)
models <- mget(model_names)
aictab(models, modnames = model_names)

# precipitation variable --> Precipitation Amount is best (wt = 0.60)
m1 <- glmmTMB(PlasmaDetection ~ DaysSinceLastPrecipitation_5mm + (1|Event/Site), 
              data = leye.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ PrecipitationAmount_7days + (1|Event/Site), 
              data = leye.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)

# precipitation variable --> Precipitation Amount is best (wt = 0.92)
m1 <- glmmTMB(PlasmaDetection ~ DaysSinceLastPrecipitation_5mm + (1|Site), 
              data = leye.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ PrecipitationAmount_7days + (1|Site), 
              data = leye.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)


# date variable --> Season is better but both are within 2 delta AICc
m1 <- glmmTMB(PlasmaDetection ~ Season + (1|Event/Site), 
              data = leye.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ Julian + (1|Event/Site), 
              data = leye.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)

# date variable --> Julian is better but both are within 2 delta AICc
m1 <- glmmTMB(PlasmaDetection ~ Season + (1|Site), 
              data = leye.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ Julian + (1|Site), 
              data = leye.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)


# species variable --> Migration Status is better (wt. = 0.88)
m1 <- glmmTMB(PlasmaDetection ~ MigStatus + (1|Event/Site), 
              data = leye.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ Species + (1|Event/Site), 
              data = leye.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)

# species variable --> Migration Status is better (wt. = 0.66)
m1 <- glmmTMB(PlasmaDetection ~ MigStatus + (1|Site), 
              data = leye.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ Species + (1|Site), 
              data = leye.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)

#------------------------------------------------------------------------------#
#       model selection exercise (event/site as a nested random effect)     ----
#------------------------------------------------------------------------------#
# event is extremely influential....but I can't include it in every model
# because it's correlated with pretty much everything! Decided to use
# a nested random effect of event and site

m1 <- glmmTMB(PlasmaDetection ~ NearestCropDistance_m + (1|Event/Site), 
              family = "binomial", 
              data = leye.cs)
m2 <- glmmTMB(PlasmaDetection ~ SPEI + (1|Event/Site), 
              family = "binomial", 
              data = leye.cs) 
m3 <- glmmTMB(PlasmaDetection ~ AnnualSnowfall_in + (1|Event/Site), 
              family = "binomial", 
              data = leye.cs)
m4 <- glmmTMB(PlasmaDetection ~ PrecipitationAmount_7days + (1|Event/Site), 
              family = "binomial", 
              data = leye.cs) 
m5 <- glmmTMB(PlasmaDetection ~ Percent_Exposed_Shoreline + (1|Event/Site), 
              family = "binomial", 
              data = leye.cs)
m6 <- glmmTMB(PlasmaDetection ~ Permanence + (1|Event/Site), 
              family = "binomial", 
              data = leye.cs)
m7 <- glmmTMB(PlasmaDetection ~ MigStatus + (1|Event/Site), 
              family = "binomial", 
              data = leye.cs)
# m8 IS THE BEST MODEL
m8 <- glmmTMB(PlasmaDetection ~ Sex + (1|Event/Site), 
              family = "binomial", 
              data = leye.cs)
m9 <- glmmTMB(PlasmaDetection ~ EnvDetection + (1|Event/Site), 
              family = "binomial", 
              data = leye.cs)
m10 <- glmmTMB(PlasmaDetection ~ Julian + (1|Event/Site), 
              family = "binomial", 
              data = leye.cs)
m11 <- glmmTMB(PlasmaDetection ~ (1|Event/Site), 
              family = "binomial", 
              data = leye.cs) 

model_names <- paste0("m", 1:11)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#     K   AICc Delta_AICc AICcWt Cum.Wt     LL
# m8  4 172.50       0.00   0.41   0.41 -82.12
# m11 3 175.17       2.67   0.11   0.52 -84.51
# m1  4 175.45       2.95   0.09   0.62 -83.60
# m10 4 176.12       3.63   0.07   0.68 -83.94
# m4  4 176.16       3.66   0.07   0.75 -83.95
# m2  4 176.46       3.97   0.06   0.81 -84.11
# m7  4 176.62       4.12   0.05   0.86 -84.19
# m3  4 176.97       4.48   0.04   0.90 -84.36
# m9  4 177.22       4.72   0.04   0.94 -84.48
# m5  4 177.27       4.78   0.04   0.98 -84.51
# m6  5 178.59       6.10   0.02   1.00 -84.11

# model diagnostics
# residuals do not like great...but that's expected with little variance in RE
simulationOutput <- simulateResiduals(fittedModel = m8) 
plot(simulationOutput)
testDispersion(m8) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 


#------------------------------------------------------------------------------#
#            model selection exercise (site as a random effect)             ----
#------------------------------------------------------------------------------#
m1 <- glmmTMB(PlasmaDetection ~ NearestCropDistance_m + (1|Site), 
              family = "binomial", 
              data = leye.cs)
m2 <- glmmTMB(PlasmaDetection ~ SPEI + (1|Site), 
              family = "binomial", 
              data = leye.cs) 
m3 <- glmmTMB(PlasmaDetection ~ AnnualSnowfall_in + (1|Site), 
              family = "binomial", 
              data = leye.cs)
m4 <- glmmTMB(PlasmaDetection ~ PrecipitationAmount_7days + (1|Site), 
              family = "binomial", 
              data = leye.cs) 
m5 <- glmmTMB(PlasmaDetection ~ Percent_Exposed_Shoreline + (1|Site), 
              family = "binomial", 
              data = leye.cs)
m6 <- glmmTMB(PlasmaDetection ~ Permanence + (1|Site), 
              family = "binomial", 
              data = leye.cs)
m7 <- glmmTMB(PlasmaDetection ~ MigStatus + (1|Site), 
              family = "binomial", 
              data = leye.cs)
m8 <- glmmTMB(PlasmaDetection ~ Sex + (1|Site), 
              family = "binomial", 
              data = leye.cs)
m9 <- glmmTMB(PlasmaDetection ~ EnvDetection + (1|Site), 
              family = "binomial", 
              data = leye.cs)
m10 <- glmmTMB(PlasmaDetection ~ Julian + (1|Site), 
               family = "binomial", 
               data = leye.cs)
# residuals look good
m11 <- glmmTMB(PlasmaDetection ~ Event + (1|Site), 
               family = "binomial", 
               data = leye.cs)
m12 <- glmmTMB(PlasmaDetection ~ 1 + (1|Site), 
               family = "binomial", 
               data = leye.cs) 

model_names <- paste0("m", 1:12)
models <- mget(model_names)
aictab(models, modnames = model_names)

# there are no other informative parameters. Just sampling event

# Model selection based on AICc:
#   
#     K   AICc Delta_AICc AICcWt Cum.Wt     LL
# m11 5 165.07       0.00      1      1 -77.35
# m10 3 178.05      12.98      0      1 -85.95
# m2  3 179.82      14.74      0      1 -86.83
# m4  3 180.74      15.66      0      1 -87.29
# m8  3 184.16      19.08      0      1 -89.00
# m1  3 186.16      21.08      0      1 -90.00
# m12 2 186.67      21.59      0      1 -91.30
# m6  4 186.98      21.91      0      1 -89.37
# m7  3 187.70      22.62      0      1 -90.77
# m3  3 188.46      23.38      0      1 -91.15
# m5  3 188.66      23.58      0      1 -91.26
# m9  3 188.74      23.66      0      1 -91.29

# model diagnostics
simulationOutput <- simulateResiduals(fittedModel = m11) 
plot(simulationOutput)
testDispersion(m11) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 


