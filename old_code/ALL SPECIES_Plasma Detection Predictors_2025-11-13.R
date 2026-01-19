#---------------------------------#
# Predictors of Neonics in Plasma #
#   Analysis by Shelby McCahon    #
#      Created: 2025-11-13        #
#     Modified: 2025-11-13        #
#---------------------------------#

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
library(performance)

# read in data
birds <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")

# manipulate data
birds <- birds %>% 
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
  DominantCrop = factor(DominantCrop,
                        levels = c("Grassland",
                                   "Canola",
                                   "Soybean",
                                   "Wheat")),
  Event = factor(Event,
                 levels = c("Fall 2021",
                            "Spring 2022",
                            "Spring 2023",
                            "Fall 2023"),
                 labels = c("2021 Fall",
                            "2022 Spring",
                            "2023 Spring",
                            "2023 Fall")),
  Species = factor(Species,
                   levels = c("American Avocet", 
                              "Killdeer", 
                              "Least Sandpiper",
                              "Lesser Yellowlegs",
                              "Longbilled Dowitcher",
                              "Marbled Godwit",
                              "Pectoral Sandpiper",
                              "Semipalmated Sandpiper",
                              "Short-billed",
                              "Willet",
                              "Wilsons Phalarope"),
                   labels = c("American Avocet", 
                              "Killdeer", 
                              "Least Sandpiper",
                              "Lesser Yellowlegs",
                              "Long-billed Dowitcher",
                              "Marbled Godwit",
                              "Pectoral Sandpiper",
                              "Semipalmated Sandpiper",
                              "Short-billed Dowitcher",
                              "Willet",
                              "Wilson's Phalarope")))
  
birds$Site <- as.factor(birds$Site)

# filter data
birds <- birds %>%
  filter(!is.na(PlasmaDetection)) %>% 
  group_by(Site) %>% 
  filter(n() >= 3) %>% # removed 7 birds
  ungroup()
  
# standardize data except for response
birds.cs <- birds %>%
  mutate(across(where(is.numeric) & !matches("PlasmaDetection"), scale))


#------------------------------------------------------------------------------#
#                          identify correlations                            ----
#------------------------------------------------------------------------------#

 cor_results <- birds %>%
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
#  Var1                           Var2                      Correlation
#  seconds_since_midnight         time_hours                      1.000    
#  Biomass                        NearestCropDistance_m           0.985
#  FatteningIndex                 Tri                             0.783
#  Beta                           FatteningIndex                 -0.783
#  DaysSinceLastPrecipitation_5mm PrecipitationAmount_7days      -0.765
#  PecSizeBest                    Standardized.Pec.NoEvent        0.752
#  Julian                         SPEI                           -0.672
#  DaysSinceLastPrecipitation_5mm Julian                          0.628
#------------------------------------------------------------------------------#

### ...identify correlations among two categorical variables                ----

#------------------------------------------------------------------------------#
# season and permanence
# event and everything...
#------------------------------------------------------------------------------#
kruskal.test(AnnualSnowfall_in ~ Event, data = birds) # corr
kruskal.test(SPEI ~ Event, data = birds) # corr
kruskal.test(PrecipitationAmount_7days ~ Event, data = birds) # corr
kruskal.test(Julian ~ Event, data = birds) # corr
kruskal.test(Percent_Exposed_Shoreline ~ Event, data = birds) # corr

chisq.test(table(birds$Event, birds$Permanence)) # corr
chisq.test(table(birds$Event, birds$MigStatus))  # corr
chisq.test(table(birds$Event, birds$Sex)) # corr
chisq.test(table(birds$Event, birds$EnvDetection)) # almost corr


chisq.test(table(birds$Season, birds$Permanence)) 
chisq.test(table(birds$Event, birds$Permanence)) 
chisq.test(table(birds$Event, birds$DominantCrop))

m <- glm(PlasmaDetection ~ Event + DominantCrop, data = birds,
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
t.test(Julian ~ Season, data = birds)
wilcox.test(Julian ~ Season, data = birds)
wilcox.test(SPEI ~ Season, data = birds)

# Kruskal-Wallis test (3+ groups, ANOVA is parametric)
summary(aov(AnnualSnowfall_in ~ Event, data = birds))
kruskal.test(Julian ~ Event, data = birds)
kruskal.test(SPEI ~ Season, data = birds)
kruskal.test(AnnualSnowfall_in ~ Event, data = birds)

kruskal.test(PercentAg ~ DominantCrop, data = birds) 
kruskal.test(PercentAg ~ NearestCropDistance_m, data = birds)

#------------------------------------------------------------------------------#
#               identify ideal random effect structure                      ----
#------------------------------------------------------------------------------#

m1 <- glmmTMB(PlasmaDetection ~ (1|Event/Site), 
              data = birds.cs, 
              family = binomial(link = "logit"),
              REML = TRUE)

m2 <- glmmTMB(PlasmaDetection ~ (1|Event), 
              data = birds.cs,
              family = binomial(link = "logit"),
              REML = TRUE)

m3 <- glmmTMB(PlasmaDetection ~ (1|Site),
              data = birds.cs,
              family = binomial(link = "logit"),
              REML = TRUE)

m4 <- glmmTMB(PlasmaDetection ~ 1,
              data = birds.cs,
              family = binomial(link = "logit"),
              REML = TRUE)

m5 <- glmmTMB(PlasmaDetection ~ (1|Event) + (1|Event:Site), 
              data = birds.cs,
              family = binomial(link = "logit"),
              REML = TRUE)

m6 <- glmmTMB(PlasmaDetection ~ (1|Event:Site), 
              data = birds.cs,
              family = binomial(link = "logit"),
              REML = TRUE)

model_names <- paste0("m", 1:6)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#   K   AICc Delta_AICc AICcWt Cum.Wt      LL
# m2 2 174.89       0.00   0.53   0.53  -85.41
# m1 3 175.17       0.28   0.46   1.00  -84.51
# m3 2 186.67      11.77   0.00   1.00  -91.30
# m4 1 215.99      41.09   0.00   1.00 -106.98

birds.cs$SiteWithinEvent <- interaction(birds.cs$Event, birds.cs$Site)

levels(birds.cs$SiteWithinEvent)       # shows unique site-event combinations
length(levels(birds.cs$SiteWithinEvent))  # counts them
sum(table(birds.cs$SiteWithinEvent) > 0)


#------------------------------------------------------------------------------#
#    model selection for corr. variables  (site & event as a random effect) ----
#------------------------------------------------------------------------------#

# agriculture variable --> Nearest Crop Distance is best (wt. = 0.45)
m1 <- glmmTMB(PlasmaDetection ~ PercentAg + (1|Event/Site), 
              data = birds.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ DominantCrop + (1|Event/Site), 
              data = birds.cs,
              family = binomial(link = "logit"))

m3 <- glmmTMB(PlasmaDetection ~ NearestCropDistance_m + (1|Event/Site),
              data = birds.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:3)
models <- mget(model_names)
aictab(models, modnames = model_names)

# agriculture variable --> Nearest Crop Distance is best (wt. = 0.68)
m1 <- glmmTMB(PlasmaDetection ~ PercentAg + (1|Site), 
              data = birds.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ DominantCrop + (1|Site), 
              data = birds.cs,
              family = binomial(link = "logit"))

m3 <- glmmTMB(PlasmaDetection ~ NearestCropDistance_m + (1|Site),
              data = birds.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:3)
models <- mget(model_names)
aictab(models, modnames = model_names)

# precipitation variable --> Precipitation Amount is best (wt = 0.60)
m1 <- glmmTMB(PlasmaDetection ~ DaysSinceLastPrecipitation_5mm + (1|Event/Site), 
              data = birds.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ PrecipitationAmount_7days + (1|Event/Site), 
              data = birds.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)

# precipitation variable --> Precipitation Amount is best (wt = 0.92)
m1 <- glmmTMB(PlasmaDetection ~ DaysSinceLastPrecipitation_5mm + (1|Site), 
              data = birds.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ PrecipitationAmount_7days + (1|Site), 
              data = birds.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)


# date variable --> Season is better but both are within 2 delta AICc
m1 <- glmmTMB(PlasmaDetection ~ Season + (1|Event/Site), 
              data = birds.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ Julian + (1|Event/Site), 
              data = birds.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)

# date variable --> Julian is better but both are within 2 delta AICc
m1 <- glmmTMB(PlasmaDetection ~ Season + (1|Site), 
              data = birds.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ Julian + (1|Site), 
              data = birds.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)


# species variable --> Migration Status is better (wt. = 0.88)
m1 <- glmmTMB(PlasmaDetection ~ MigStatus + (1|Event/Site), 
              data = birds.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ Species + (1|Event/Site), 
              data = birds.cs,
              family = binomial(link = "logit"))

model_names <- paste0("m", 1:2)
models <- mget(model_names)
aictab(models, modnames = model_names)

# species variable --> Migration Status is better (wt. = 0.66)
m1 <- glmmTMB(PlasmaDetection ~ MigStatus + (1|Site), 
              data = birds.cs, 
              family = binomial(link = "logit"))

m2 <- glmmTMB(PlasmaDetection ~ Species + (1|Site), 
              data = birds.cs,
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
              data = birds.cs)
m2 <- glmmTMB(PlasmaDetection ~ SPEI + (1|Event/Site), 
              family = "binomial", 
              data = birds.cs) 
m3 <- glmmTMB(PlasmaDetection ~ AnnualSnowfall_in + (1|Event/Site), 
              family = "binomial", 
              data = birds.cs)
m4 <- glmmTMB(PlasmaDetection ~ PrecipitationAmount_7days + (1|Event/Site), 
              family = "binomial", 
              data = birds.cs) 
m5 <- glmmTMB(PlasmaDetection ~ Percent_Exposed_Shoreline + (1|Event/Site), 
              family = "binomial", 
              data = birds.cs)
m6 <- glmmTMB(PlasmaDetection ~ Permanence + (1|Event/Site), 
              family = "binomial", 
              data = birds.cs)
m7 <- glmmTMB(PlasmaDetection ~ MigStatus + (1|Event/Site), 
              family = "binomial", 
              data = birds.cs)
# m8 IS THE BEST MODEL
m8 <- glmmTMB(PlasmaDetection ~ Sex + (1|Event/Site), 
              family = "binomial", 
              data = birds.cs)
m9 <- glmmTMB(PlasmaDetection ~ EnvDetection + (1|Event/Site), 
              family = "binomial", 
              data = birds.cs)
m10 <- glmmTMB(PlasmaDetection ~ Julian + (1|Event/Site), 
              family = "binomial", 
              data = birds.cs)
m11 <- glmmTMB(PlasmaDetection ~ (1|Event/Site), 
              family = "binomial", 
              data = birds.cs) 
m12 <- glmmTMB(PlasmaDetection ~ time_hours + (1|Event/Site), 
               family = "binomial", 
               data = birds.cs) 

model_names <- paste0("m", 1:12)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
  
#     K   AICc Delta_AICc AICcWt Cum.Wt     LL
# m8  4 172.50       0.00   0.39   0.39 -82.12
# m11 3 175.17       2.67   0.10   0.50 -84.51
# m1  4 175.45       2.95   0.09   0.59 -83.60
# m10 4 176.12       3.63   0.06   0.65 -83.94
# m4  4 176.16       3.66   0.06   0.71 -83.95
# m2  4 176.46       3.97   0.05   0.77 -84.11
# m7  4 176.62       4.12   0.05   0.82 -84.19
# m12 4 176.75       4.25   0.05   0.87 -84.25
# m3  4 176.97       4.48   0.04   0.91 -84.36
# m9  4 177.22       4.72   0.04   0.95 -84.48
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
              data = birds.cs)
m2 <- glmmTMB(PlasmaDetection ~ SPEI + (1|Site), 
              family = "binomial", 
              data = birds.cs) 
m3 <- glmmTMB(PlasmaDetection ~ AnnualSnowfall_in + (1|Site), 
              family = "binomial", 
              data = birds.cs)
m4 <- glmmTMB(PlasmaDetection ~ PrecipitationAmount_7days + (1|Site), 
              family = "binomial", 
              data = birds.cs) 
m5 <- glmmTMB(PlasmaDetection ~ Percent_Exposed_Shoreline + (1|Site), 
              family = "binomial", 
              data = birds.cs)
m6 <- glmmTMB(PlasmaDetection ~ Permanence + (1|Site), 
              family = "binomial", 
              data = birds.cs)
m7 <- glmmTMB(PlasmaDetection ~ MigStatus + (1|Site), 
              family = "binomial", 
              data = birds.cs)
m8 <- glmmTMB(PlasmaDetection ~ Sex + (1|Site), 
              family = "binomial", 
              data = birds.cs)
m9 <- glmmTMB(PlasmaDetection ~ EnvDetection + (1|Site), 
              family = "binomial", 
              data = birds.cs)
m10 <- glmmTMB(PlasmaDetection ~ Julian + (1|Site), 
               family = "binomial", 
               data = birds.cs)
# residuals look good
m11 <- glmmTMB(PlasmaDetection ~ Event + (1|Site), 
               family = "binomial", 
               data = birds.cs)
m12 <- glmmTMB(PlasmaDetection ~ time_hours + (1|Site), 
               family = "binomial", 
               data = birds.cs)
m13 <- glmmTMB(PlasmaDetection ~ 1 + (1|Site), 
               family = "binomial", 
               data = birds.cs) 

model_names <- paste0("m", 1:13)
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
# m13 2 186.67      21.59      0      1 -91.30
# m6  4 186.98      21.91      0      1 -89.37
# m7  3 187.70      22.62      0      1 -90.77
# m3  3 188.46      23.38      0      1 -91.15
# m5  3 188.66      23.58      0      1 -91.26
# m12 3 188.71      23.64      0      1 -91.28
# m9  3 188.74      23.66      0      1 -91.29

# model diagnostics
simulationOutput <- simulateResiduals(fittedModel = m11) 
plot(simulationOutput)
testDispersion(m11) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 


# test with multiple linear regression
m1 <- glmmTMB(PlasmaDetection ~ AnnualSnowfall_in + Event + (1|Site), 
              family = "binomial", 
              data = birds.cs)
summary(m1)
check_collinearity(m1)
