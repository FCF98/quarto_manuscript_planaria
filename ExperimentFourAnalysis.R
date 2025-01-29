#### Experiment 4 - conditioning with improved Y mazes without divot


library(readxl)
library(tidyr)
library(dplyr)
library(lme4)
library(ggplot2)
library(hrbrthemes)
library(effectsize)
library(emmeans)
library(DHARMa)
library(ggsignif)
library(patchwork)
library(ez)
library(rstatix)

Exp4_data <- read_excel("Datasets/ExperimentFourSummaryData.xlsx")

# Restructure data
Exp4_data_long <- Exp4_data %>%
  mutate(across(ends_with("_arm_%"), 
                ~as.numeric(as.character(.)))) %>%
  pivot_longer(
    cols = c(`baseline_active_arm_%`, `endpoint_active_arm_%`),
    names_to = "Time",
    values_to = "ActiveCount"
  ) %>%
  mutate(
    Time = factor(Time, 
                  levels = c("baseline_active_arm_%", "endpoint_active_arm_%"),
                  labels = c("Baseline", "Endpoint")),
    Subject = factor(Subject),
    ActiveCount = if_else(!is.na(ActiveCount), 
                          as.integer(ActiveCount), 
                          NA_integer_)
  ) %>%
  mutate(
    TotalTrials = 6,
    InactiveCount = if_else(!is.na(ActiveCount),
                            TotalTrials - ActiveCount,
                            NA_integer_),
    ActiveArmProportion = if_else(!is.na(ActiveCount),
                                  ActiveCount / TotalTrials,
                                  NA_real_)
  ) %>%
  dplyr::filter(is.na(ActiveCount) | ActiveCount <= TotalTrials)

# Set contrasts
contrasts(Exp4_data_long$Time) <- contr.sum

# Fit simplified model for within-subject comparison
Exp4_m1 <- glmer(cbind(ActiveCount, InactiveCount) ~ Time + (1|Subject), 
                 data = Exp4_data_long, 
                 family = "binomial",
                 na.action = na.omit)

# Print model summary
summary(Exp4_m1)

# Get overall effect of Time
exp4_decisions_model_output <- car::Anova(Exp4_m1, type = "III")
print(exp4_decisions_model_output)

# Get estimated means and comparison
exp4_time_comparison <- emmeans(Exp4_m1, 
                                specs = pairwise ~ Time,
                                adjust = "bonferroni",
                                type = "response")

# Print results
print(exp4_time_comparison)


######### Plotting learning across conditioning days ##########



# Reshape the data for plotting
Exp4_Data_long_days <- Exp4_data %>%  # Changed from 'data' to 'Exp4_data'
  select(Subject,
         `Baseline_day1`, 
         `Baseline_day2`, 
         `Conditioning_day1`, 
         `Conditioning_day2`,
         `Conditioning_day3`, 
         `Conditioning_day4`, 
         `Conditioning_day5`) %>%
  pivot_longer(
    cols = c(`Baseline_day1`:`Conditioning_day5`),
    names_to = "TimePoint",
    values_to = "ActiveArmChoices"
  ) %>%
  mutate(
    TimePoint = factor(TimePoint,
                       levels = c("Baseline_day1", "Baseline_day2", 
                                  "Conditioning_day1", "Conditioning_day2", 
                                  "Conditioning_day3", "Conditioning_day4", 
                                  "Conditioning_day5")),
    Proportion = ActiveArmChoices / 3
  ) %>%
  group_by(TimePoint) %>%
  summarise(
    mean_prop = mean(Proportion, na.rm = TRUE),
    se = sd(Proportion, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Create the visualization
ggplot(Exp4_Data_long_days, aes(x = TimePoint, y = mean_prop)) +
  geom_line(linewidth = 1.5, color = "#D55E00", group = 1) +
  geom_point(size = 3, color = "#D55E00") +
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
                width = 0.2, color = "#D55E00") +
  theme_minimal() +
  labs(
    title = "Average Proportion of Active Arm Choices Across Days",
    y = "Proportion of Active Arm Choices",
    x = "Day"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylim(0, 1)



######## assessing significance of active arm entries across days #######


# First reshape the data to wide format for repeated measures ANOVA
Exp4_repeated_data <- Exp4_data %>%
  select(Subject,
         `Baseline_day1`, 
         `Baseline_day2`, 
         `Conditioning_day1`, 
         `Conditioning_day2`,
         `Conditioning_day3`, 
         `Conditioning_day4`, 
         `Conditioning_day5`) %>%
  mutate(across(-Subject, ~./3)) # Convert to proportions

# Conduct repeated measures ANOVA
# First need to get data in correct format for ezANOVA
long_data <- Exp4_repeated_data %>%
  pivot_longer(-Subject, 
               names_to = "TimePoint", 
               values_to = "Proportion") %>%
  mutate(TimePoint = factor(TimePoint, 
                            levels = c("Baseline_day1", "Baseline_day2",
                                       "Conditioning_day1", "Conditioning_day2",
                                       "Conditioning_day3", "Conditioning_day4",
                                       "Conditioning_day5")))

# Run repeated measures ANOVA using ez package
library(ez)
rm_anova <- ezANOVA(
  data = long_data,
  dv = Proportion,
  wid = Subject,
  within = TimePoint,
  detailed = TRUE
)

# Print ANOVA results
print(rm_anova)

# If ANOVA is significant, conduct planned comparisons
# Let's compare baseline vs conditioning days using paired t-tests

# Average baseline
baseline_avg <- rowMeans(Exp4_repeated_data[,c("Baseline_day1", "Baseline_day2")], na.rm = TRUE)

# Average last two conditioning days
final_cond_avg <- rowMeans(Exp4_repeated_data[,c("Conditioning_day4", "Conditioning_day5")], na.rm = TRUE)

# Paired t-test
t_test_result <- t.test(baseline_avg, final_cond_avg, paired = TRUE)

print(t_test_result)




# Run repeated measures ANOVA
rm_anova <- ezANOVA(
  data = long_data,
  dv = Proportion,
  wid = Subject,
  within = TimePoint,
  detailed = TRUE
)

# Print ANOVA results
print(rm_anova)

# Conduct pairwise comparisons with p-value adjustment
pwc <- long_data %>%
  pairwise_t_test(
    Proportion ~ TimePoint,
    paired = TRUE,
    p.adjust.method = "bonferroni"  # You can change this to "fdr" or other methods
  )

# Print pairwise comparisons
print(pwc, n = 22)

# You can also do specific planned comparisons
# Example: Compare average baseline to each conditioning day
baseline_comparisons <- long_data %>%
  filter(TimePoint %in% c("Baseline_day1", "Baseline_day2")) %>%
  group_by(Subject) %>%
  summarise(baseline_mean = mean(Proportion, na.rm = TRUE)) %>%
  full_join(
    long_data %>%
      filter(grepl("Conditioning", TimePoint)),
    by = "Subject"
  ) %>%
  group_by(TimePoint) %>%
  summarise(
    t_stat = t.test(Proportion, baseline_mean, paired = TRUE)$statistic,
    p_value = t.test(Proportion, baseline_mean, paired = TRUE)$p.value
  )

print(baseline_comparisons)