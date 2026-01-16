#---------------------------------#
#  Predictors of Neonics in Water #
#   Analysis by Shelby McCahon    #
#      Created: 2026-01-15        #
#     Modified: 2026-01-16        #
#---------------------------------#

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
# seasonal wetlands all had detections, so I combined it with temporary wetlands
# event highly influential so I included it in all models and used it as informed null

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
# combine temporary and seasonal permanence classes
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

# standardize data
water.cs <- water %>%
  mutate(across(where(is.numeric), scale))

# set reference levels
water.cs$DominantCrop <- relevel(water.cs$DominantCrop,
                                 ref = "Grassland")

water.cs$Permanence <- relevel(water.cs$Permanence,
                               ref = "Temporary/Seasonal")

#------------------------------------------------------------------------------#
#                   identify best model for each hypothesis                 ----                        
#------------------------------------------------------------------------------# 

# include all variables that show up in models with support (delta AICc < 2)

# best agricultural model
global.model <- glm(WaterNeonicDetection ~ PercentAg + DominantCrop +
                      NearestCropDistance_m + Buffered + Event,
                    family = "binomial",
                    data = water.cs,
                    na.action = na.fail)

# check for collinearity
car::vif(global.model) # vif < 4

dredge(global.model)

# Model selection table 
#       (Int) Bff DmC Evn    NCD_m     PrA df  logLik  AICc delta weight
# 30 -0.19340   +       +  2.12000 0.97510  7 -35.248  85.8  0.00  0.439
# 29 -0.80140           +  0.98180 1.21700  6 -37.177  87.3  1.53  0.204
# 14  0.01437   +       +  1.69800          6 -37.607  88.1  2.39  0.133
# 21 -0.78440           +          0.49630  5 -39.601  89.9  4.11  0.056
# 5  -0.69310           +                   4 -41.112  90.7  4.90  0.038
# 16 -0.06119   +   +   +  2.26500         10 -34.255  91.1  5.31  0.031
# 6  -0.61710   +       +                   5 -40.498  91.7  5.90  0.023
# 22 -0.80410   +       +          0.53500  6 -39.591  92.1  6.36  0.018
# 32 -0.16120   +   +   +  2.36800 0.61750 11 -33.680  92.5  6.71  0.015
# 13 -0.68740           +  0.06482          5 -41.085  92.8  7.07  0.013
# 31 -0.87890       +   +  1.05100 0.90680 10 -35.551  93.7  7.90  0.008
# 7  -0.85090       +   +                   8 -38.051  93.7  7.98  0.008
# 15 -0.89210       +   +  0.51900          9 -36.874  93.8  8.06  0.008
# 23 -0.84710       +   +          0.08525  9 -38.028  96.1 10.37  0.002
# 8  -0.86120   +   +   +                   9 -38.048  96.2 10.41  0.002
# 24 -0.88620   +   +   +          0.15790 10 -37.998  98.6 12.80  0.001
# 26 -0.59680   +          1.44400 0.73450  4 -48.093 104.6 18.86  0.000
# 10 -0.34060   +          1.18400          3 -49.871 106.0 20.25  0.000
# 25 -1.24200              0.70910 1.04900  3 -50.107 106.5 20.72  0.000
# 17 -1.18700                      0.54900  2 -51.908 107.9 22.19  0.000
# 2  -0.86750   +                           2 -52.412 109.0 23.20  0.000
# 18 -1.06600   +                  0.41010  3 -51.749 109.8 24.00  0.000
# 1  -1.11200                               1 -54.270 110.6 24.83  0.000
# 12 -0.79060   +   +      1.31400          7 -48.339 111.9 26.18  0.000
# 9  -1.11300             -0.05361          2 -54.246 112.6 26.86  0.000
# 28 -0.74630   +   +      1.44300 0.57010  8 -47.620 112.9 27.12  0.000
# 3  -1.60900       +                       5 -51.689 114.0 28.28  0.000
# 27 -1.45700       +      0.67940 0.85780  7 -49.599 114.5 28.70  0.000
# 4  -1.29900   +   +                       6 -51.206 115.3 29.59  0.000
# 19 -1.42900       +              0.33940  6 -51.241 115.4 29.66  0.000
# 11 -1.72800       +      0.23580          6 -51.337 115.6 29.85  0.000
# 20 -1.29200   +   +              0.20800  7 -51.091 117.4 31.68  0.000
# Models ranked by AICc(x) 

# can remove dominant crop
m1 <- glm(WaterNeonicDetection ~ PercentAg + NearestCropDistance_m + 
            Buffered + Event,
          family = "binomial",
          data = water.cs)

summary(m1)
confint(m1)

# model validation --> some pattern but n.s.
simulationOutput <- simulateResiduals(fittedModel = m1) 
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m1)$PercentAg) #  good
plotResiduals(simulationOutput, form = model.frame(m1)$NearestCropDistance_m) # some pattern
plotResiduals(simulationOutput, form = model.frame(m1)$Buffered)  # good
plotResiduals(simulationOutput, form = model.frame(m1)$Event) # good

#------------------------------------------------------------------------------#

# best wetland dynamics model
global.model <- glm(WaterNeonicDetection ~ Permanence + 
                      Dist_Closest_Wetland_m + Event,
                    family = "binomial",
                    na.action = na.fail,
                    data = water.cs)

# check for collinearity
car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#     (Int) Dst_Cls_Wtl_m Evn Prm df  logLik  AICc delta weight
# 7  1.1090                 +   +  6 -38.247  89.4  0.00  0.467
# 3 -0.6931                 +      4 -41.112  90.7  1.23  0.252
# 8  1.0010       -0.1591   +   +  7 -38.101  91.5  2.03  0.169
# 4 -0.7623       -0.2115   +      5 -40.835  92.3  2.90  0.109
# 5 -0.3514                     +  3 -47.226 100.7 11.28  0.002
# 6 -0.3663       -0.3415       +  4 -46.522 101.5 12.05  0.001
# 1 -1.1120                        1 -54.270 110.6 21.15  0.000
# 2 -1.1330       -0.2869          2 -53.724 111.6 22.15  0.000
# Models ranked by AICc(x) 

# distance to nearest wetland can be removed
m2 <- glm(WaterNeonicDetection ~ Permanence + Event,
          data = water.cs,
          family = "binomial",
          na.action = na.fail)

summary(m3)
confint(m3)

# model validation --> pattern but not significant
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$Permanence) #  good
plotResiduals(simulationOutput, form = model.frame(m2)$Event) # good


#------------------------------------------------------------------------------#

# best precipitation model
global.model <- glm(WaterNeonicDetection ~ AnnualSnowfall_in + 
                      AnnualPrecipitation_in + DaysSinceLastPrecipitation_5mm +
                      PrecipitationAmount_7days + Event,
                    family = "binomial",
                    na.action = na.fail,
                    data = water.cs)

# check for collinearity
car::vif(global.model) # vif < 2

dredge(global.model)

# Model selection table 
#      (Int)   AnP_in  AnS_in DSL_5mm Evn   PrA_7dy df  logLik  AICc delta weight
# 9  -0.6931                            +            4 -41.112  90.7  0.00  0.120
# 10 -1.0650 -0.50470                   +            5 -40.130  90.9  0.26  0.106
# 11 -1.6890          -0.5883           +            5 -40.134  90.9  0.27  0.105
# 26 -1.4470 -0.66470                   + -0.506000  6 -39.319  91.6  0.91  0.076
# 27 -2.2100          -0.7427           + -0.470700  6 -39.389  91.7  1.05  0.071
# 30 -1.4950 -0.65150         -0.5013   + -0.813800  7 -38.278  91.8  1.16  0.068
# 31 -2.2450          -0.7332 -0.4669   + -0.761600  7 -38.326  91.9  1.25  0.064
# 13 -0.6369                  -0.2568   +            5 -40.655  92.0  1.31  0.063
# 25 -0.8481                            + -0.287800  5 -40.792  92.2  1.59  0.055
# 29 -0.8984                  -0.5059   + -0.599300  6 -39.663  92.3  1.60  0.054
# 15 -1.5730          -0.5460 -0.2050   +            6 -39.829  92.6  1.93  0.046
# 14 -0.9936 -0.46090         -0.2008   +            6 -39.857  92.6  1.99  0.045
# 12 -1.5680 -0.32280 -0.3769           +            6 -39.866  92.7  2.01  0.044
# 28 -2.1210 -0.44910 -0.4741           + -0.555000  7 -38.924  93.1  2.45  0.035
# 32 -2.1650 -0.43950 -0.4706 -0.4781   + -0.852400  8 -37.890  93.4  2.76  0.030
# 16 -1.4780 -0.28960 -0.3608 -0.1892   +            7 -39.614  94.5  3.83  0.018
# 1  -1.1120                                         1 -54.270 110.6 19.92  0.000
# 3  -1.1340          -0.2925                        2 -53.505 111.1 20.48  0.000
# 5  -1.1290                  -0.2598                2 -53.835 111.8 21.14  0.000
# 7  -1.1530          -0.2988 -0.2773                3 -53.045 112.3 21.69  0.000
# 17 -1.1130                               0.061350  2 -54.235 112.6 21.94  0.000
# 2  -1.1120 -0.01712                                2 -54.268 112.7 22.00  0.000
# 4  -1.1360  0.11490 -0.3374                        3 -53.405 113.1 22.41  0.000
# 19 -1.1340          -0.2896              0.039110  3 -53.491 113.2 22.58  0.000
# 21 -1.1320                  -0.3125     -0.075440  3 -53.799 113.9 23.20  0.000
# 6  -1.1290 -0.03454         -0.2640                3 -53.824 113.9 23.25  0.000
# 23 -1.1590          -0.3093 -0.3640     -0.112500  4 -52.963 114.4 23.70  0.000
# 8  -1.1540  0.09180 -0.3335 -0.2652                4 -52.980 114.4 23.74  0.000
# 18 -1.1140 -0.04647                      0.078260  3 -54.219 114.7 24.04  0.000
# 20 -1.1360  0.11970 -0.3400             -0.009243  4 -53.404 115.2 24.58  0.000
# 22 -1.1320 -0.01186         -0.3102     -0.070060  4 -53.798 116.0 25.37  0.000
# 24 -1.1660  0.18150 -0.3875 -0.4126     -0.205000  5 -52.768 116.2 25.54  0.000
# Models ranked by AICc(x) 

# retain all predictors
m3 <- glm(WaterNeonicDetection ~ AnnualSnowfall_in + 
                      DaysSinceLastPrecipitation_5mm + AnnualPrecipitation_in +
                      PrecipitationAmount_7days + Event,
                    family = "binomial",
                    na.action = na.fail,
                    data = water.cs)

summary(m3)
confint(m3)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$AnnualSnowfall_in) #  good
plotResiduals(simulationOutput, form = model.frame(m3)$DaysSinceLastPrecipitation_5mm) # pattern but n.s.
plotResiduals(simulationOutput, form = model.frame(m3)$Event) # good
plotResiduals(simulationOutput, form = model.frame(m3)$PrecipitationAmount_7days) # good
plotResiduals(simulationOutput, form = model.frame(m3)$AnnualPrecipitation_in) # good

#------------------------------------------------------------------------------#

# best drought model
m4 <- glm(WaterNeonicDetection ~ SPEI + Event,
          data = water.cs,
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

plotResiduals(simulationOutput, form = model.frame(m4)$Event) # good
plotResiduals(simulationOutput, form = model.frame(m4)$SPEI) # good

#------------------------------------------------------------------------------#

# informed null
m5 <- glm(WaterNeonicDetection ~ Event,
          data = water.cs,
          family = "binomial",
          na.action = na.fail)

summary(m5)
confint(m5)

# model validation --> good
simulationOutput <- simulateResiduals(fittedModel = m5) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m5)$Event) #  good

#------------------------------------------------------------------------------#

# STAGE II: Cross-hypothesis

# AIC model selection
model_names <- paste0("m", 1:5)
models <- mget(model_names)
aictab(models, modnames = model_names)

# Model selection based on AICc:
#   
#    K  AICc Delta_AICc AICcWt Cum.Wt     LL
# m1 7 85.76       0.00   0.75   0.75 -35.25 ** agriculture
# m2 6 89.43       3.67   0.12   0.87 -38.25 ** wetland dynamics
# m5 4 90.66       4.90   0.06   0.93 -41.11 ** informed null
# m4 5 91.03       5.28   0.05   0.98 -40.19
# m3 8 93.42       7.66   0.02   1.00 -37.89

#------------------------------------------------------------------------------#

# STAGE III: model average
model_avg <- model.avg(m1, m2, m5)

summary(model_avg)
confint(model_avg)

# Model-averaged coefficients:  
#   (full average) 
#                         Estimate Std. Error Adjusted SE z value Pr(>|z|)  
# (Intercept)             -0.06133    0.97855     0.98884   0.062   0.9505  
# PercentAg                0.78288    0.57832     0.58267   1.344   0.1791  
# NearestCropDistance_m    1.70172    1.17990     1.18779   1.433   0.1520  
# BufferedY               -3.34274    3.01470     3.04344   1.098   0.2721  
# EventFall 2023          -3.30963    1.37320     1.39139   2.379   0.0174 *
# EventSpring 2022         1.18100    1.26840     1.28094   0.922   0.3565  
# EventSpring 2023         0.21751    1.00498     1.01501   0.214   0.8303  
# PermanencePermanent     -0.27734    0.82045     0.82291   0.337   0.7361  
# PermanenceSemipermanent -0.17899    0.57655     0.57923   0.309   0.7573  

#                               2.5 %      97.5 %
# (Intercept)             -1.99941429  1.87676140
# PercentAg                0.02418989  1.92607568
# NearestCropDistance_m    0.28982713  3.94937678
# BufferedY               -9.74840583  1.42121911
# EventFall 2023          -6.03670711 -0.58256270
# EventSpring 2022        -1.32959000  3.69159964
# EventSpring 2023        -1.77187785  2.20690714
# PermanencePermanent     -4.31061211 -0.02372223
# PermanenceSemipermanent -3.27432296  0.47704909



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
  term = c("Sampling Event (Fall 2021)", 
           "Wetland Permanence (Temporary/Seasonal)",
           "No Buffer Presence"), 
  estimate = 0,
  conf_low = 0,
  conf_high = 0
)

df <- rbind(df, ref_levels)

# relabel terms
df$term <- factor(
  df$term,
  levels = rev(c("NearestCropDistance_m",
                 "PercentAg",
                 "BufferedY",
                 "No Buffer Presence",
                 "EventSpring 2022",
                 "EventSpring 2023",
                 "EventFall 2023",
                 "Sampling Event (Fall 2021)",
                 "PermanenceSemipermanent",
                 "PermanencePermanent",
                 "Wetland Permanence (Temporary/Seasonal)")),
  labels = rev(c("Nearest Crop Distance",
                 "% Surrounding Cropland Cover",
                 "Buffer Presence",
                 "No Buffer Presence",
                 "Sampling Event (Spring 2022)",
                 "Sampling Event (Spring 2023)",
                 "Sampling Event (Fall 2023)",
                 "Sampling Event (Fall 2021)",
                 "Wetland Permanence (Semi-permanent)",
                 "Wetland Permanence (Permanent)",
                 "Wetland Permanence (Temporary/Seasonal)")))

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
  xlim(-10,5) +
  labs(
    x = "Standardized Parameter Estimate",
    y = "")


# plot relationships
ggplot(water, aes(y = NearestCropDistance_m,
                  x = WaterNeonicDetection)) +
  geom_boxplot() + geom_point()

cor(water$NearestCropDistance_m, water$PercentAg) # -0.72
