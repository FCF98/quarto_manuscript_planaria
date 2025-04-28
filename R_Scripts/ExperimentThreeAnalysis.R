library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest) # for p-values in mixed models



# I need to read in the data

Exp3_data <- read.csv("Datasets/ExperimentThreeData.csv")


#remove low moving planaria

# Calculate 5th percentile thresholds
baseline_threshold <- quantile(Exp3_data$baseline_distance, 0.05)
test_threshold <- quantile(Exp3_data$Test_distance, 0.05)

# Filter the data
filtered_data <- Exp3_data[
  Exp3_data$baseline_distance > baseline_threshold & 
    Exp3_data$Test_distance > test_threshold,
]

# First, let's reshape the data to long format

# Create long format data
long_data <- tidyr::pivot_longer(
  filtered_data,
  cols = c("Baseline_active_proportion", "Test_active_proportion"),
  names_to = "Time",
  values_to = "Active_Preference"
)

# Clean up the Time variable
long_data$Time <- factor(ifelse(long_data$Time == "Baseline_active_proportion", 
                                "Baseline", "Test"))

# Fit mixed effects model
model <- lmer(Active_Preference ~ Condition * Time + (1|Subject), 
              data = long_data)

# Create visualization
ggplot(long_data, aes(x = Time, y = Active_Preference, 
                      color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_point(alpha = 0.3, position = position_jitter(width = 0.1)) +
  theme_classic() +
  scale_color_manual(values = c("blue", "red")) +
  labs(
    title = "Change in Active Preference Over Time",
    y = "Proportion of Time in Active Surface",
    x = "Time Point"
  )

# Print model summary
print("Mixed Effects Model Results:")
print(summary(model))

# Calculate and print means for each condition and time point
print("\nMean Active Preference by Condition and Time:")
aggregate(Active_Preference ~ Condition + Time, data = long_data, 
          FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))