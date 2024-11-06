

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



# I need to check assumption of equal variances

levene_dose_response <- leveneTest(Distance ~ Condition, dose_response_data)


# I need to check assumption of normality

shapiro_test_dose_response <- shapiro.test(dose_response_data$Distance)


# because assumption of normality is violated, I need to perform a non-parametric test

kruskal_test <- kruskal.test(Distance ~ Condition, data = dose_response_data)





######## anova format for future use ##########

  #dose_response_anova <- aov_car(Distance ~ Condition + Error(Subject), data = dose_response_data)



