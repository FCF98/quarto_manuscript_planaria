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
  mutate(across(ends_with("_arm_entries"), 
                ~as.numeric(as.character(.)))) %>%
  pivot_longer(
    cols = c(`baseline_active_arm_entries`, `endpoint_active_arm_entries`),
    names_to = "Time",
    values_to = "ActiveCount"
  ) %>%
  mutate(
    Time = factor(Time, 
                  levels = c("baseline_active_arm_entries", "endpoint_active_arm_entries"),
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






######## Visualisation of 10 subjects selected to be bisected #######


# Reshape the data for plotting with selected subjects
Exp4_Data_long_days <- Exp4_data %>%
  # Filter for specific subjects
  filter(Subject %in% c(1, 2, 3, 6, 17, 22, 28, 38, 39, 42)) %>%
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
    title = "Average Proportion of Active Arm Choices Across Days\n(Selected Subjects)",
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



# Read the Excel file
Ten_subjects_data <- read_excel("for claude v1.xlsx")

# Define the selected subjects
selected_subjects <- c(1, 2, 36, 17, 22, 28, 38, 39, 42)

# Filter for selected subjects
selected_data <- Ten_subjects_data %>%
  filter(Subject %in% selected_subjects)

# Function to calculate daily proportions (grouping trials by day)
calculate_daily_prop <- function(data, day_columns) {
  data %>%
    rowwise() %>%
    summarize(
      active_entries = sum(across(all_of(day_columns)) == `Active arm`, na.rm = TRUE),
      total_trials = sum(!is.na(across(all_of(day_columns)))),
      proportion = active_entries / total_trials
    ) %>%
    pull(proportion)
}

# Group baseline trials by day (3 trials per day)
baseline_day1 <- c("B1", "B2", "B3")
baseline_day2 <- c("B4", "B5", "B6")

# Group conditioning trials by day (3 trials per day)
conditioning_day1 <- c("C1", "C2", "C3")
conditioning_day2 <- c("C4", "C5", "C6")
conditioning_day3 <- c("C7", "C8", "C9")
conditioning_day4 <- c("C10", "C11", "C12")
conditioning_day5 <- c("C13", "C14", "C15")

# Calculate proportions for each day
baseline_props <- data.frame(
  Day = c("Baseline_day1", "Baseline_day2"),
  Proportion = c(
    mean(calculate_daily_prop(selected_data, baseline_day1)),
    mean(calculate_daily_prop(selected_data, baseline_day2))
  ),
  Phase = "Baseline",
  TimePoint = 1:2
)

conditioning_props <- data.frame(
  Day = paste0("Conditioning_day", 1:5),
  Proportion = c(
    mean(calculate_daily_prop(selected_data, conditioning_day1)),
    mean(calculate_daily_prop(selected_data, conditioning_day2)),
    mean(calculate_daily_prop(selected_data, conditioning_day3)),
    mean(calculate_daily_prop(selected_data, conditioning_day4)),
    mean(calculate_daily_prop(selected_data, conditioning_day5))
  ),
  Phase = "Conditioning",
  TimePoint = 3:7
)

# Combine all training data
all_days <- rbind(baseline_props, conditioning_props)

# Process regenerate data
regenerate_data <- Ten_subjects_data %>%
  filter(grepl("[HT]$", Subject)) %>%
  mutate(Type = ifelse(grepl("H$", Subject), "Head", "Tail"))

# Function to calculate SEM for regenerates
calculate_regenerate_sem <- function(data, phase_cols) {
  # Calculate proportion for each subject
  subject_props <- data %>%
    rowwise() %>%
    summarize(
      active_entries = sum(across(all_of(phase_cols)) == `Active arm`, na.rm = TRUE),
      total_trials = sum(!is.na(across(all_of(phase_cols)))),
      proportion = active_entries / total_trials
    ) %>%
    pull(proportion)
  
  # Calculate SEM
  sd(subject_props, na.rm = TRUE) / sqrt(length(subject_props))
}

# Calculate test and reinstatement data with SEM
test_data <- regenerate_data %>%
  group_by(Type) %>%
  summarize(
    Proportion = mean(calculate_phase_prop(regenerate_data[regenerate_data$Type == first(Type), ], c("T1", "T2", "T3"))),
    SEM = calculate_regenerate_sem(regenerate_data[regenerate_data$Type == first(Type), ], c("T1", "T2", "T3"))
  ) %>%
  mutate(TimePoint = 8,
         Phase = "Test")

reinst_data <- regenerate_data %>%
  group_by(Type) %>%
  summarize(
    Proportion = mean(calculate_phase_prop(regenerate_data[regenerate_data$Type == first(Type), ], c("R1", "R2", "R3"))),
    SEM = calculate_regenerate_sem(regenerate_data[regenerate_data$Type == first(Type), ], c("R1", "R2", "R3"))
  ) %>%
  mutate(TimePoint = 9,
         Phase = "Reinstatement")

create_combined_plot <- function(type_label) {
  # Ensure the last conditioning day has the same columns as test and reinstatement data
  last_conditioning <- transform(tail(all_days, 1), 
                                 TimePoint = 7,
                                 Phase = "Conditioning")
  
  combined_data <- bind_rows(
    last_conditioning,
    test_data %>% filter(Type == type_label),
    reinst_data %>% filter(Type == type_label)
  )
  
  ggplot() +
    # Training data line and points
    geom_line(data = all_days, 
              aes(x = TimePoint, y = Proportion),
              color = "orange", size = 1) +
    geom_point(data = all_days,
               aes(x = TimePoint, y = Proportion),
               color = "orange", size = 3) +
    geom_errorbar(data = all_days,
                  aes(x = TimePoint, 
                      y = Proportion,
                      ymin = Proportion - sd(Proportion)/sqrt(length(selected_subjects)),
                      ymax = Proportion + sd(Proportion)/sqrt(length(selected_subjects))),
                  width = 0.2, color = "orange") +
    # Test and reinstatement points with error bars
    geom_line(data = combined_data,
              aes(x = TimePoint, y = Proportion),
              color = "orange") +
    geom_point(data = test_data %>% filter(Type == type_label),
               aes(x = TimePoint, y = Proportion),
               color = "orange", size = 3) +
    geom_errorbar(data = test_data %>% filter(Type == type_label),
                  aes(x = TimePoint,
                      ymin = Proportion - SEM,
                      ymax = Proportion + SEM),
                  width = 0.2, color = "orange") +
    geom_point(data = reinst_data %>% filter(Type == type_label),
               aes(x = TimePoint, y = Proportion),
               color = "orange", size = 3) +
    geom_errorbar(data = reinst_data %>% filter(Type == type_label),
                  aes(x = TimePoint,
                      ymin = Proportion - SEM,
                      ymax = Proportion + SEM),
                  width = 0.2, color = "orange") +
    scale_x_continuous(breaks = 1:9,
                       labels = c("Baseline_day1", "Baseline_day2",
                                  "Conditioning_day1", "Conditioning_day2", 
                                  "Conditioning_day3", "Conditioning_day4",
                                  "Conditioning_day5",
                                  "Test", "Reinstatement")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = paste(type_label, "Regenerate Performance"),
         x = "Day",
         y = "Proportion of Active Arm Choices") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}
# Create and save the plots
p1 <- create_combined_plot("Head")
p2 <- create_combined_plot("Tail")
