

library(readxl)
library(dplyr)
library(car)
library(pwr)


# ANOVA_Analysis.R

# Load necessary libraries
library(dplyr)
library(afex)
library(apa)
library(papaja)

# Load the data
dose_response_data <- read.csv('Datasets/DoseResponseStatistics.csv')

# Convert Condition to factor & Need a subject number for ANOVA
dose_response_data <- dose_response_data %>%
  mutate(Condition = as.factor(Condition)) %>%
  mutate(Subject = seq(1, nrow(dose_response_data)))


# Run ANOVA
dose_response_anova_afex <- aov_car(Distance ~ Condition + Error(Subject), data = dose_response_data)

# Create apa object
apa_anova_results <- apa_print(dose_response_anova_afex)




