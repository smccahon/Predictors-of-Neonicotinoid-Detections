#-----------------------------------#
#  Predictors of Neonics (plotting) #
#     Analysis by Shelby McCahon    #
#        Created: 2026-01-30        #
#       Modified: 2026-01-30        #
#-----------------------------------#

# load packages
library(ggplot2)
library(tidyverse)
library(trtools)
library(patchwork)
library(emmeans)

#------------------------------------------------------------------------------#
#                               WATER ANALYSIS                              ----
#------------------------------------------------------------------------------#

### prep water data ----
water <- read.csv("cleaned_data/wetland_data_cleaned_2026-01-28.csv")

# reorder event
water$Event <- factor(
  water$Event,
  levels = c("Fall 2021",
             "Spring 2022",
             "Spring 2023",
             "Fall 2023"),
  labels = c("2021 Fall",
             "2022 Spring",
             "2023 Spring",
             "2023 Fall"))

# combine temporary and seasonal wetland permanence classes
water <- water %>% 
  mutate(Permanence = case_when(
    Permanence %in% c("Temporary", "Seasonal") ~ "Temporary/Seasonal",
    Permanence == "Semipermanent" ~ "Semi-permanent",
    Permanence == "Permanent" ~ "Permanent"))

# reorder permanence
water$Permanence <- factor(
  water$Permanence,
  levels = c("Temporary/Seasonal",
             "Semi-permanent",
             "Permanent"),
  labels = c("Temp./Seas.",
             "Semi-perm.",
             "Perm."))

# convert character to factor
water <- water %>% 
  mutate_if(is.character, as.factor)

# filter out wetlands without neonic information
water <- water %>%
  filter(!is.na(WaterNeonicDetection))

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
  mutate(LogPrecipAmount = PrecipitationAmount_7days_cm + 0.0001) %>% 
  mutate(LogPrecipAmount = log(LogPrecipAmount))

# log transform precipitation events due to skew
water <- water %>% 
  mutate(LogPrecipDays = DaysSinceLastPrecipitation_5mm + 0.0001) %>% 
  mutate(LogPrecipDays = log(LogPrecipDays))

# set reference levels
water$DominantCrop <- relevel(water$DominantCrop,
                                 ref = "Grassland")

water$Permanence <- relevel(water$Permanence,
                               ref = "Temp./Seas.")

### top water models ----

m.ag <- glm(WaterNeonicDetection ~ LogCropDistance + Event, 
            data = water,
            family = "binomial")

m.precip <- glm(WaterNeonicDetection ~ AnnualPrecipitation_cm + Event,
                family = "binomial",
                data = water)

m.drought <- glm(WaterNeonicDetection ~ SPEI + Event,
                 family = "binomial",
                 data = water)

m.time <- glm(WaterNeonicDetection ~ Event,
              data = water,
              family = "binomial",
              na.action = na.fail)

m.hydrology <- glm(WaterNeonicDetection ~ Permanence + Event,
                   family = "binomial",
                   data = water)

#------------------------------------------------------------------------------#
#                                 PLOTTING                                  ----
#------------------------------------------------------------------------------#

# color theme for event
cols <- c(
  "2021 Fall"   = "#875907",  
  "2023 Fall"   = "#FADBA4", 
  "2023 Spring" = "#98CCF0", 
  "2022 Spring" = "#0D3957"  
)

###... ag model ----
d <- expand.grid(
  LogCropDistance = seq(
    min(water$LogCropDistance),
    max(water$LogCropDistance),
    length = 1000),
  Event = "2022 Spring")

d <- cbind(d, trtools::glmint(m.ag, newdata = d))

head(d)

# colored by xaxis
w.ag <- ggplot(d, aes(x = LogCropDistance, y = fit)) +
  geom_line(linewidth = 1.2, col = "goldenrod3") + 
  geom_ribbon(aes(ymin = low, ymax = upp), 
              alpha = 0.4, color = NA, fill = "goldenrod3", show.legend = FALSE) +
  theme_classic() +
  labs(x ="Log(Dist. to Nearest Crop)", 
       y = "Detection Probability in Water (%)") +
  theme(axis.title.x = element_text(size = 12,
                                    margin = margin(t = 5)),
        axis.title.y = element_text(size = 12,
                                    margin = margin(r = 5)),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(expand = c(0,0), 
                     limits = c(-10,7.8))

###... hydrology model ----
d <- expand.grid(
  Permanence = unique(water$Permanence),
  Event = "2022 Spring")

d <- cbind(d, trtools::glmint(m.hydrology, newdata = d))

head(d)

# colored by xaxis
w.hydrology <- ggplot(d, aes(x = Permanence, y = fit)) +
  geom_point(size = 4, col = "blue4") + 
  geom_errorbar(aes(ymin = low, ymax = upp), width = 0,
                col = "blue4",
                size = 1) +
  theme_classic() +
  labs(x ="Wetland Permanence", 
       y = NULL) +
  theme(axis.title.x = element_text(size = 12,
                                    margin = margin(t = 5)),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(linetype = "dashed"),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(
    labels = c("Temp./Seas." = "Temp./\nSeas.",
               "Semi-perm." = "Semi-\nperm.",
               "Perm." = "Perm.")
  )


###... time model ----
d <- expand.grid(
  Event = unique(water$Event))

d <- cbind(d, trtools::glmint(m.time, newdata = d))

head(d)

# colored by xaxis
w.time <- ggplot(d, aes(x = Event, y = fit)) +
  geom_point(size = 4, col = "seagreen3") + 
  geom_errorbar(aes(ymin = low, ymax = upp), width = 0,
                col = "seagreen3",
                size = 1) +
  theme_classic() +
  labs(x ="Sampling Event", 
       y = NULL) +
  theme(axis.title.x = element_text(size = 12,
                                    margin = margin(t = 5)),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(linetype = "dashed"),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(
    labels = c("2021 Fall" = "2021\nFall",
               "2022 Spring" = "2022\nSpring",
               "2023 Spring" = "2023\nSpring",
               "2023 Fall" = "2023\nFall"))

###... precip model ----
d <- expand.grid(
  AnnualPrecipitation_cm = seq(
    min(water$AnnualPrecipitation_cm),
    max(water$AnnualPrecipitation_cm),
    length = 1000),
  Event = "2022 Spring")

d <- cbind(d, trtools::glmint(m.precip, newdata = d))

head(d)

# colored by xaxis
w.precip <- ggplot(d, aes(x = AnnualPrecipitation_cm, y = fit)) +
  geom_line(linewidth = 1.2, col = "skyblue3") + 
  geom_ribbon(aes(ymin = low, ymax = upp), 
              alpha = 0.4, color = NA, fill = "skyblue3", show.legend = FALSE) +
  theme_classic() +
  labs(x ="Annual Precipitation (cm)", 
       y = NULL) +
  theme(axis.title.x = element_text(size = 12,
                                    margin = margin(t = 5)),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(linetype = "dashed"),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(expand = c(0,0), 
                     limits = c(34,63),
                     breaks = c(35,40,45,50,55,60))


###... drought model ----
d <- expand.grid(
  SPEI = seq(
    min(water$SPEI),
    max(water$SPEI),
    length = 1000),
  Event = "2022 Spring")

d <- cbind(d, trtools::glmint(m.drought, newdata = d))

head(d)

# colored by xaxis
w.drought <- ggplot(d, aes(x = SPEI, y = fit)) +
  geom_line(linewidth = 1.2, col = "lightsalmon4") + 
  geom_ribbon(aes(ymin = low, ymax = upp), 
              alpha = 0.4, color = NA, fill = "lightsalmon4", show.legend = FALSE) +
  theme_classic() +
  labs(x ="Drought (SPEI)", 
       y = NULL) +
  theme(axis.title.x = element_text(size = 12,
                                    margin = margin(t = 5)),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(linetype = "dashed"),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(expand = c(0,0), 
                     limits = c(-2,2.5),
                     breaks = c(-2, -1, 0, 1, 2, 3))


###... combine plots ----

w.all <-  w.ag + w.hydrology + w.time + w.precip + w.drought + plot_layout(ncol = 5)





#------------------------------------------------------------------------------#
#                               INVERT ANALYSIS                             ----
#------------------------------------------------------------------------------#

### prep invert data ----
invert <- read.csv("cleaned_data/invert_data_cleaned_2025-08-11.csv")

# reorder event
invert$Season <- factor(
  invert$Season,
  levels = c("Spring",
             "Fall"))

# reorder permanence
invert$Permanence <- factor(
  invert$Permanence,
  levels = c("Temporary",
             "Seasonal",
             "Semipermanent",
             "Permanent"),
  labels = c("Temp.",
             "Seas.",
             "Semi-perm.",
             "Perm."))

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
  mutate(LogPrecipAmount = PrecipitationAmount_7days_cm) %>% 
  mutate(LogPrecipAmount = log(LogPrecipAmount))

# log transform precipitation events due to skew
invert <- invert %>% 
  mutate(LogPrecipDays = DaysSinceLastPrecipitation_5mm + 0.0001) %>% 
  mutate(LogPrecipDays = log(LogPrecipDays))

# set reference levels
invert$DominantCrop <- relevel(invert$DominantCrop,
                                  ref = "Grassland")

invert$Permanence <- relevel(invert$Permanence,
                                ref = "Temp.")

### top invert models ----

m.ag <- glm(InvertPesticideDetection ~ LogCropDistance, 
            data = invert,
            family = "binomial")

m.hydrology <- glm(InvertPesticideDetection ~ Permanence,
                   family = "binomial",
                   data = invert)

m.waterquality <- glm(InvertPesticideDetection ~ pH_probe,
                      family = "binomial",
                      data= invert)

m.precip <- glm(InvertPesticideDetection ~ PrecipitationAmount_7days_cm,
                family = "binomial",
                data = invert)

m.time <- glm(InvertPesticideDetection ~ Season,
              data = invert,
              family = "binomial")


#------------------------------------------------------------------------------#
#                                 PLOTTING                                  ----
#------------------------------------------------------------------------------#

###... precip model ----
d <- expand.grid(
  PrecipitationAmount_7days_cm = seq(
    min(invert$PrecipitationAmount_7days_cm),
    max(invert$PrecipitationAmount_7days_cm),
    length = 1000))

d <- cbind(d, trtools::glmint(m.precip, newdata = d))

head(d)

# colored by xaxis
i.precip <- ggplot(d, aes(x = PrecipitationAmount_7days_cm, y = fit)) +
  geom_line(linewidth = 1.2, col = "skyblue3") + 
  geom_ribbon(aes(ymin = low, ymax = upp), 
              alpha = 0.4, color = NA, fill = "skyblue3", show.legend = FALSE) +
  theme_classic() +
  labs(x ="Precipitation in Last 7 Days (cm)", 
       y = "Detection Probability in Invertebrates (%)") +
  theme(axis.title.x = element_text(size = 12,
                                    margin = margin(t = 5)),
        axis.title.y = element_text(size = 12,
                                    margin = margin(r = 5)),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(expand = c(0,0), 
                     limits = c(-0.15,5.5),
                     breaks = c(0,1,2,3,4,5))


###... ag model ----
d <- expand.grid(
  LogCropDistance = seq(
    min(invert$LogCropDistance),
    max(invert$LogCropDistance),
    length = 1000))

d <- cbind(d, trtools::glmint(m.ag, newdata = d))

head(d)

# colored by xaxis
i.ag <- ggplot(d, aes(x = LogCropDistance, y = fit)) +
  geom_line(linewidth = 1.2, col = "goldenrod3") + 
  geom_ribbon(aes(ymin = low, ymax = upp), 
              alpha = 0.4, color = NA, fill = "goldenrod3", show.legend = FALSE) +
  theme_classic() +
  labs(x ="Log(Dist. to Nearest Crop)", 
       y = NULL) +
  theme(axis.title.x = element_text(size = 12,
                                    margin = margin(t = 5)),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(linetype = "dashed"),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(expand = c(0,0), 
                     limits = c(-9.8,7.2))

###... hydrology model ----
d <- expand.grid(
  Permanence = unique(invert$Permanence))

d <- cbind(d, trtools::glmint(m.hydrology, newdata = d))

head(d)

# colored by xaxis
i.hydrology <- ggplot(d, aes(x = Permanence, y = fit)) +
  geom_point(size = 4, col = "blue4") + 
  geom_errorbar(aes(ymin = low, ymax = upp), width = 0,
                col = "blue4",
                size = 1) +
  theme_classic() +
  labs(x ="Wetland Permanence", 
       y = NULL) +
  theme(axis.title.x = element_text(size = 12,
                                    margin = margin(t = 5)),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(linetype = "dashed"),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(
    labels = c("Temp." = "Temp.",
               "Seas." = "Seas.",
               "Semi-perm." = "Semi-\nperm.",
               "Perm."= "Perm.")
  )

###... time model ----
d <- expand.grid(
  Season = unique(invert$Season))

d <- cbind(d, trtools::glmint(m.time, newdata = d))

head(d)

# colored by xaxis
i.time <- ggplot(d, aes(x = Season, y = fit)) +
  geom_point(size = 4, col = "seagreen3") + 
  geom_errorbar(aes(ymin = low, ymax = upp), width = 0,
                col = "seagreen3",
                size = 1) +
  theme_classic() +
  labs(x ="Season", 
       y = NULL) +
  theme(axis.title.x = element_text(size = 12,
                                    margin = margin(t = 5)),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(linetype = "dashed"),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = scales::percent_format(accuracy = 1))



###... pH model ----
d <- expand.grid(
  pH_probe = seq(
    min(invert$pH_probe),
    max(invert$pH_probe),
    length = 1000))

d <- cbind(d, trtools::glmint(m.waterquality, newdata = d))

head(d)

# colored by xaxis
i.waterquality <- ggplot(d, aes(x = pH_probe, y = fit)) +
  geom_line(linewidth = 1.2, col = "thistle4") + 
  geom_ribbon(aes(ymin = low, ymax = upp), 
              alpha = 0.4, color = NA, fill = "thistle4", show.legend = FALSE) +
  theme_classic() +
  labs(x ="pH", 
       y = NULL) +
  theme(axis.title.x = element_text(size = 12,
                                    margin = margin(t = 5)),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(linetype = "dashed"),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = scales::percent_format(accuracy = 1))


###... combine plots ----

i.all <-  i.precip + i.hydrology + i.waterquality +
  i.ag + i.time + plot_layout(ncol = 5)


w.all/plot_spacer()/i.all +
  plot_layout(heights = c(1, 0.2, 1))

ggsave(filename = "Water & Invert Detection Figure_2026-01-30.tiff",
       path = "C:/Users/shelb/OneDrive - University of Idaho/Masters Project/Analysis/Predictors of Neonics/figures", 
       width = 13, height = 8, units = "in", device='tiff', dpi=300)
