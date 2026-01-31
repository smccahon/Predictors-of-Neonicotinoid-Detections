#--------------------------------------------------#
#  Predictors of Neonics (plotting for shorebirds) #
#              Analysis by Shelby McCahon          #
#                 Created: 2026-01-30              #
#                Modified: 2026-01-30              #
#--------------------------------------------------#

# load packages
library(ggplot2)
library(tidyverse)
library(trtools)
library(patchwork)
library(emmeans)

#------------------------------------------------------------------------------#
#                               LEYE ANALYSIS                               ----
#------------------------------------------------------------------------------#

### prep leye data ----
leye <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")

# combine temporary and seasonal permanence classes
leye <- leye %>% 
  mutate(Permanence = case_when(
    Permanence %in% c("Temporary", "Seasonal") ~ "Temporary/Seasonal",
    Permanence == "Semipermanent" ~ "Semipermanent",
    Permanence == "Permanent" ~ "Permanent"))

# reorder detection
leye$EnvDetection <- factor(
  leye$EnvDetection,
  levels = c("Y","N"))

# Convert Julian to months
month_starts <- as.Date(paste0("2023-", 1:12, "-01"))  # Non-leap year example
julian_breaks <- as.numeric(format(month_starts, "%j"))  
month_labels <- format(month_starts, "%B") # use %b for first three letters
month_labels[month_labels == "September"] <- "Sept."

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
  mutate(LogPrecipAmount = PrecipitationAmount_7days_cm + 0.0001) %>% 
  mutate(LogPrecipAmount = log(LogPrecipAmount))

# log transform precipitation event due to right skew
leye <- leye %>% 
  mutate(LogPrecipDays = DaysSinceLastPrecipitation_5mm) %>% 
  mutate(LogPrecipDays = log(LogPrecipDays))

# set reference levels
leye$DominantCrop <- relevel(leye$DominantCrop,
                                ref = "Grassland")

leye$Permanence <- relevel(leye$Permanence,
                              ref = "Temporary/Seasonal")

### top leye models ----

m.ag <- glm(PlasmaDetection ~ PercentAg + Julian,
            family = "binomial",
            data = leye)

m.pesticide <- glm(PlasmaDetection ~ EnvDetection + Julian,
                   family = "binomial",
                   data = leye)

m.precip <- glm(PlasmaDetection ~ AnnualPrecipitation_cm + Julian,
                family = "binomial",
                data = leye)

m.drought <- glm(PlasmaDetection ~ SPEI + Julian,
                 family = "binomial",
                 data = leye)

m.lifehistory <- glm(PlasmaDetection ~ Age + Julian,
                     family = "binomial",
                     data = leye)

m.hydrology <- glm(PlasmaDetection ~ LogWetlandDistance + Julian,
                   family = "binomial",
                   data=leye)

m.time <- glm(PlasmaDetection ~ Julian,
              family = "binomial",
              data = leye)


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
  PercentAg = seq(
    min(leye$PercentAg),
    max(leye$PercentAg),
    length = 1000),
  Julian = mean(leye$Julian))

d <- cbind(d, trtools::glmint(m.ag, newdata = d))

head(d)

# colored by xaxis
l.ag <- ggplot(d, aes(x = PercentAg, y = fit)) +
  geom_line(linewidth = 1.2, col = "goldenrod3") + 
  geom_ribbon(aes(ymin = low, ymax = upp), 
              alpha = 0.4, color = NA, fill = "goldenrod3", show.legend = FALSE) +
  theme_classic() +
  labs(x ="Surrounding Cropland Cover (%)", 
       y = "Detection Probability\nin Lesser Yellowlegs Plasma (%)") +
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
                     labels = scales::percent_format(accuracy = 1))

###... pesticide exposure model ----
d <- expand.grid(
  EnvDetection = unique(leye$EnvDetection),
  Julian = mean(leye$Julian))

d <- cbind(d, trtools::glmint(m.pesticide, newdata = d))

head(d)

# colored by xaxis
l.pesticide <- ggplot(d, aes(x = EnvDetection, y = fit)) +
  geom_point(size = 4, col = "orangered4") + 
  geom_errorbar(aes(ymin = low, ymax = upp), width = 0,
                col = "orangered4",
                size = 1) +
  theme_classic() +
  labs(x ="Wetland Pesticide Detection", 
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
    labels = c("Y" = "Yes",
               "N" = "No")
  )

###... time model ----
d <- expand.grid(
  Julian = seq(
    min(leye$Julian),
    max(leye$Julian),
    length = 1000))

d <- cbind(d, trtools::glmint(m.time, newdata = d))

head(d)

# colored by xaxis
l.time <- ggplot(d, aes(x = Julian, y = fit)) +
  geom_line(linewidth = 1.2, col = "seagreen3") + 
  geom_ribbon(aes(ymin = low, ymax = upp), 
              alpha = 0.4, color = NA, fill = "seagreen3", show.legend = FALSE) +
  theme_classic() +
  labs(x ="Month", 
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
  scale_x_continuous(breaks = julian_breaks, labels = month_labels)

###... precipitation model ----
d <- expand.grid(
  AnnualPrecipitation_cm = seq(
    min(leye$AnnualPrecipitation_cm),
    max(leye$AnnualPrecipitation_cm),
    length = 1000),
  Julian = mean(leye$Julian)
)

d <- cbind(d, trtools::glmint(m.precip, newdata = d))

head(d)

# colored by xaxis
l.precip <- ggplot(d, aes(x = AnnualPrecipitation_cm, y = fit)) +
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
                     breaks = c(40,45,50,55,60))
  

###... hydrology model ----
d <- expand.grid(
  LogWetlandDistance = seq(
    min(leye$LogWetlandDistance),
    max(leye$LogWetlandDistance),
    length = 1000),
  Julian = mean(leye$Julian)
)

d <- cbind(d, trtools::glmint(m.hydrology, newdata = d))

head(d)

# colored by xaxis
l.hydrology <- ggplot(d, aes(x = LogWetlandDistance, y = fit)) +
  geom_line(linewidth = 1.2, col = "blue4") + 
  geom_ribbon(aes(ymin = low, ymax = upp), 
              alpha = 0.4, color = NA, fill = "blue4", show.legend = FALSE) +
  theme_classic() +
  labs(x ="Log(Dist. to Nearest Wetland)", 
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

###... drought model ----
d <- expand.grid(
  SPEI = seq(
    min(leye$SPEI),
    max(leye$SPEI),
    length = 1000),
  Julian = mean(leye$Julian)
)

d <- cbind(d, trtools::glmint(m.drought, newdata = d))

head(d)

# colored by xaxis
l.drought <- ggplot(d, aes(x = SPEI, y = fit)) +
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
                     labels = scales::percent_format(accuracy = 1))

###... life history model ----
d <- expand.grid(
  Age = unique(leye$Age),
  Julian = mean(leye$Julian))

d <- cbind(d, trtools::glmint(m.lifehistory, newdata = d))

head(d)

# colored by xaxis
l.lifehistory <- ggplot(d, aes(x = Age, y = fit)) +
  geom_point(size = 4, col = "orchid3") + 
  geom_errorbar(aes(ymin = low, ymax = upp), width = 0,
                col = "orchid3",
                size = 1) +
  theme_classic() +
  labs(x ="Age", 
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


#------------------------------------------------------------------------------#
#                             ALL SPECIES ANALYSIS                          ----
#------------------------------------------------------------------------------#

### prep bird data ----
birds <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")

# combine temporary and seasonal permanence classes to increase sample size
birds <- birds %>% 
  mutate(Permanence = case_when(
    Permanence %in% c("Temporary", "Seasonal") ~ "Temporary/Seasonal",
    Permanence == "Semipermanent" ~ "Semipermanent",
    Permanence == "Permanent" ~ "Permanent"))

# reorder detection
birds$EnvDetection <- factor(
  birds$EnvDetection,
  levels = c("Y","N"))


# reorder event
birds$Event <- factor(
  birds$Event,
  levels = c("Fall 2021",
             "Spring 2022",
             "Spring 2023",
             "Fall 2023"),
  labels = c("2021 Fall",
             "2022 Spring",
             "2023 Spring",
             "2023 Fall"))

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

# set reference levels
birds$DominantCrop <- relevel(birds$DominantCrop,
                                 ref = "Grassland")

birds$Permanence <- relevel(birds$Permanence,
                               ref = "Temporary/Seasonal")

### top leye models ----

m.ag <- glm(PlasmaDetection ~ LogCropDistance + Event,
            family = "binomial",
            data = birds)

m.pesticide <- glm(PlasmaDetection ~ EnvDetection + Event,
                   family = "binomial",
                   data = birds)

m.precip <- glm(PlasmaDetection ~ AnnualPrecipitation_cm + Event,
                family = "binomial",
                data = birds)

m.drought <- glm(PlasmaDetection ~ SPEI + Event,
                 family = "binomial",
                 data = birds)

m.lifehistory <- glm(PlasmaDetection ~ Sex + Event,
                     family = "binomial",
                     data = birds)

m.hydrology <- glm(PlasmaDetection ~ Dist_Closest_Wetland_m + Event,
                   family = "binomial",
                   data = birds)

m.time <- glm(PlasmaDetection ~ Event,
              family = "binomial",
              data = birds)


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
    min(birds$LogCropDistance),
    max(birds$LogCropDistance),
    length = 1000),
  Event = "2022 Spring")

d <- cbind(d, trtools::glmint(m.ag, newdata = d))

head(d)

# colored by xaxis
b.ag <- ggplot(d, aes(x = LogCropDistance, y = fit)) +
  geom_line(linewidth = 1.2, col = "goldenrod3") + 
  geom_ribbon(aes(ymin = low, ymax = upp), 
              alpha = 0.4, color = NA, fill = "goldenrod3", show.legend = FALSE) +
  theme_classic() +
  labs(x ="Log(Dist. to Nearest Crop)", 
       y = "Detection Probability\nin All Species Plasma (%)") +
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
                     labels = scales::percent_format(accuracy = 1))

###... pesticide exposure model ----
d <- expand.grid(
  EnvDetection = unique(birds$EnvDetection),
  Event = "2022 Spring")

d <- cbind(d, trtools::glmint(m.pesticide, newdata = d))

head(d)

# colored by xaxis
b.pesticide <- ggplot(d, aes(x = EnvDetection, y = fit)) +
  geom_point(size = 4, col = "orangered4") + 
  geom_errorbar(aes(ymin = low, ymax = upp), width = 0,
                col = "orangered4",
                size = 1) +
  theme_classic() +
  labs(x ="Wetland Pesticide Detection", 
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
    labels = c("Y" = "Yes",
               "N" = "No")
  )

###... time model ----
d <- expand.grid(
  Event = unique(birds$Event))

d <- cbind(d, trtools::glmint(m.time, newdata = d))

head(d)

# colored by xaxis
b.time <- ggplot(d, aes(x = Event, y = fit)) +
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

###... precipitation model ----
d <- expand.grid(
  AnnualPrecipitation_cm = seq(
    min(birds$AnnualPrecipitation_cm),
    max(birds$AnnualPrecipitation_cm),
    length = 1000),
  Event = "2022 Spring"
)

d <- cbind(d, trtools::glmint(m.precip, newdata = d))

head(d)

# colored by xaxis
b.precip <- ggplot(d, aes(x = AnnualPrecipitation_cm, y = fit)) +
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

###... hydrology model ----
d <- expand.grid(
  Dist_Closest_Wetland_m = seq(
    min(birds$Dist_Closest_Wetland_m),
    max(birds$Dist_Closest_Wetland_m),
    length = 1000),
  Event = "2022 Spring"
)

d <- cbind(d, trtools::glmint(m.hydrology, newdata = d))

head(d)

# colored by xaxis
b.hydrology <- ggplot(d, aes(x = Dist_Closest_Wetland_m, y = fit)) +
  geom_line(linewidth = 1.2, col = "blue4") + 
  geom_ribbon(aes(ymin = low, ymax = upp), 
              alpha = 0.4, color = NA, fill = "blue4", show.legend = FALSE) +
  theme_classic() +
  labs(x ="Dist. to Nearest Wetland (m)", 
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

###... drought model ----
d <- expand.grid(
  SPEI = seq(
    min(birds$SPEI),
    max(birds$SPEI),
    length = 1000),
  Event = "2022 Spring"
)

d <- cbind(d, trtools::glmint(m.drought, newdata = d))

head(d)

# colored by xaxis
b.drought <- ggplot(d, aes(x = SPEI, y = fit)) +
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
                     labels = scales::percent_format(accuracy = 1))

###... life history model ----
d <- expand.grid(
  Sex = unique(birds$Sex),
  Event = "2022 Spring")

d <- cbind(d, trtools::glmint(m.lifehistory, newdata = d))

head(d)

# colored by xaxis
b.lifehistory <- ggplot(d, aes(x = Sex, y = fit)) +
  geom_point(size = 4, col = "orchid3") + 
  geom_errorbar(aes(ymin = low, ymax = upp), width = 0,
                col = "orchid3",
                size = 1) +
  theme_classic() +
  labs(x ="Sex", 
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
    labels = c("F" = "Female",
               "M" = "Male"))


###... combine plots ----

l.all <-  l.ag + l.pesticide + l.time + l.precip + l.hydrology +
  l.drought + l.lifehistory + plot_layout(ncol = 7)

b.all <-  b.ag + b.hydrology + b.time + b.lifehistory +
  b.precip + b.pesticide + b.drought +  plot_layout(ncol = 7)

b.all/plot_spacer()/l.all +
  plot_layout(heights = c(1, 0.2, 1))

ggsave(filename = "PLASMA Detection Figure_2026-01-30.tiff",
       path = "C:/Users/shelb/OneDrive - University of Idaho/Masters Project/Analysis/Predictors of Neonics/figures", 
       width = 20, height = 6, units = "in", device='tiff', dpi=300)

