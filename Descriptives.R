

library(readxl)
library(dplyr)
library(car)
library(pwr)


# ANOVA_Analysis.R

# Load necessary libraries
library(dplyr)
library(afex)
library(apa)

# Load the data
dose_response_data <- read.csv('Datasets/DoseResponseStatistics.csv')

# Convert Condition to factor
dose_response_data <- dose_response_data %>%
  mutate(Condition = as.factor(Condition))

# Run ANOVA
anova_model <- lm(Distance ~ Condition, data = dose_response_data)

dose_response_anova <- anova(anova_model)
dose_response_anova
# Save the summary of the ANOVA
anova_summary <- summary(dose_response_anova)

# Save the ANOVA summary as an R object
save(anova_summary, file = "anova_summary.RData")

anova_results <- anova_summary[[1]] 

# Extract F-value and p-value
f_value <- signif(anova_results["Condition", "F value"], 3)
p_value <- signif(anova_results["Condition", "Pr(>F)"], 3)

df1 <- anova_results$Df[1]  # Degrees of freedom
df2 <- anova_results$Df[2]  # Degrees of freedom

# Print extracted values for APA-style reporting
f_value
p_value

# I need to check assumption of equal variances

levene_dose_response <- leveneTest(Distance ~ Condition, dose_response_data)

levene_f_value <- signif(levene_dose_response$`F value`[1], 3)
levene_p_value <- signif(levene_dose_response$`Pr(>F)`[1], 3)
df1_levene <- levene_dose_response$Df[1]  # Numerator df
df2_levene <- levene_dose_response$Df[2]  # Denominator df

# I need to check assumption of normality

shapiro_dose_response <-  shapiro.test(residuals(dose_response_anova))

shapiro_dose_response

shapiro_w_value <- signif(shapiro_dose_response$statistic, 3)
shapiro_p_value <- signif(shapiro_dose_response$p.value, 3)



