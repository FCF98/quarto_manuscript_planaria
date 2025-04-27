# Planaria Y-Maze Analysis for Experiment 4
# Load required libraries
library(readr)
library(dplyr)
library(tidyverse)
library(ggsignif)
library(lme4)
library(emmeans)
library(grid)
library(gridExtra)
library(car)
library(ggplot2)
library(readxl)
library(patchwork) # For combining plots

# Read in the data
Exp4_full_data <- read_csv("Datasets/ExperimentFourFullData.csv")


#=================================================================
# Creating theme for all plots
#=================================================================

consistent_theme <- function() {
  theme_classic() +
    theme(
      text = element_text(family = "Times New Roman", size = 18),
      axis.title = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(margin = margin(r = 20)),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.text = element_text(size = 16, color = "black"),
      legend.position = "top",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      axis.line = element_line(color = "black", linewidth = 0.8),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
}

#=================================================================
# PART 1: ANALYSIS OF INTACT SUBJECTS (PRE-REGENERATION)
#=================================================================

# Filter for original intact subjects (subjects with numeric IDs)
# Use grepl to identify purely numeric subject IDs
Orig_subjects <- Exp4_full_data %>%
  filter(grepl("^[0-9]+$", as.character(Subject)))

# Restructure data for intact analysis (comparing baseline vs endpoint)
Exp4_full_data_long <- Orig_subjects %>%
  mutate(across(ends_with("_arm_entries"), 
                ~as.numeric(as.character(.)))) %>%
  pivot_longer(
    cols = c(baseline_active_arm_entries, endpoint_active_arm_entries),
    names_to = "Time",
    values_to = "ActiveCount"
  ) %>%
  mutate(
    Time = factor(Time, 
                  levels = c("baseline_active_arm_entries", "endpoint_active_arm_entries"),
                  labels = c("Baseline", "Endpoint")),
    Subject = factor(Subject),
    # Since we only have one condition, create a dummy Condition variable
    Condition = "Drug",
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
  filter(is.na(ActiveCount) | ActiveCount <= TotalTrials)

# Convert to wide format for paired t-test
Exp4_wide_format <- Exp4_full_data_long %>%
  select(Subject, Time, ActiveArmProportion) %>%
  pivot_wider(
    names_from = Time,
    values_from = ActiveArmProportion
  )

# Display wide format data structure to verify
print("Wide format data structure:")
print(str(Exp4_wide_format))


#=================================
#apa formatted results
#================================

# Perform paired t-test
paired_ttest_model <- t.test(
  Exp4_wide_format$Baseline,
  Exp4_wide_format$Endpoint,
  paired = TRUE
)

# Calculate Cohen's d for paired t-test
t_value <- paired_ttest_model$statistic
df <- paired_ttest_model$parameter  # Extract degrees of freedom
n <- length(Exp4_wide_format$Baseline)
cohens_d <- t_value / sqrt(n)

# Format p-value according to APA style
p_formatted <- ifelse(paired_ttest_model$p.value < .001, "< .001", 
                      paste0("= ", gsub("0\\.", ".", round(paired_ttest_model$p.value, 3))))

# Create ttest_results dataframe with test statistic and degrees of freedom
ttest_results <- data.frame(
  contrast = "Baseline / Endpoint",
  t_statistic = t_value,
  df = df,
  cohens_d = abs(cohens_d),
  p.value = paired_ttest_model$p.value,
  apa_result = paste0("*t*(", round(df, 2), ") = ", round(t_value, 2), 
                      ", *d* = ", round(abs(cohens_d), 2), 
                      ", *p* ", p_formatted)
)

print(ttest_results)


# Convert wide format back to long for plotting
Exp4_plot_data <- Exp4_wide_format %>%
  pivot_longer(
    cols = c(Baseline, Endpoint),
    names_to = "Time",
    values_to = "ActiveArmProportion"
  ) %>%
  mutate(Time = factor(Time, levels = c("Baseline", "Endpoint")))

# Visualization 1: Bar plot with error bars
intact_plot <- Exp4_plot_data %>%
  group_by(Time) %>%
  summarise(
    mean_prop = mean(ActiveArmProportion, na.rm = TRUE),
    se = sd(ActiveArmProportion, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = Time, y = mean_prop)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6, 
           fill = "#FF8C00") + # Using the treatment color from original script
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.8) +
  # Add significance annotation based on t-test result
  geom_signif(
    annotations = ifelse(paired_ttest_model$p.value < 0.05, "*", "ns"),
    xmin = 1, xmax = 2,
    y_position = 0.7,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  labs(
    title = "Active Arm Entries Before and After Conditioning",
    y = "Proportion of Active Arm Entries",
    x = "Time Point"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  consistent_theme()

print(intact_plot)


#=================================================================
# PART 2: LEARNING ACROSS CONDITIONING DAYS
#=================================================================

# Reshape the data for plotting conditioning days
Exp4_full_data_long_days <- Orig_subjects %>%
  select(Subject,
         Baseline_day1, 
         Baseline_day2, 
         Conditioning_day1, 
         Conditioning_day2,
         Conditioning_day3, 
         Conditioning_day4, 
         Conditioning_day5) %>%
  pivot_longer(
    cols = c(Baseline_day1:Conditioning_day5),
    names_to = "TimePoint",
    values_to = "ActiveArmChoices"
  ) %>%
  mutate(
    TimePoint = factor(TimePoint,
                       levels = c("Baseline_day1", "Baseline_day2", 
                                  "Conditioning_day1", "Conditioning_day2", 
                                  "Conditioning_day3", "Conditioning_day4", 
                                  "Conditioning_day5"),
                       labels = c("BL1", "BL2", "CD1", "CD2", "CD3", "CD4", "CD5")),
    Subject = factor(Subject),
    Proportion = ActiveArmChoices / 3 # 3 trials per day
  )

# Create the learning curve visualization
learning_plot <- Exp4_full_data_long_days %>%
  group_by(TimePoint) %>%
  summarise(
    mean_prop = mean(Proportion, na.rm = TRUE),
    se = sd(Proportion, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = TimePoint, y = mean_prop, group = 1)) +
  geom_line(linewidth = 1.5, color = "#FF8C00") +
  geom_point(size = 4, color = "#FF8C00") +
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
                width = 0.2, linewidth = 0.8, color = "#FF8C00") +
  labs(
    title = "Active Arm Preference Throughout Conditioning",
    y = "Proportion of Active Arm Entries",
    x = "Day"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  consistent_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(learning_plot)

#=================================================================
# PART 3: REGENERATION ANALYSIS - UPDATED WITH NEW COLUMN NAMES
#=================================================================

# Filter for regenerated subjects (subjects with H or T suffix)
Regen_subjects <- Exp4_full_data %>%
  filter(grepl("[HT]$", Subject)) %>%
  mutate(
    # Extract original subject number and body part
    OrigSubject = as.numeric(str_remove(Subject, "[HT]$")),
    BodyPart = str_extract(Subject, "[HT]$"),
    BodyPart = factor(BodyPart, levels = c("H", "T"), labels = c("Head", "Tail")),
    Subject = factor(Subject)
  )

# Create a dataset for the regeneration memory test phase - UPDATED
Memory_test_data <- Regen_subjects %>%
  select(Subject, OrigSubject, BodyPart, T1, T2, T3, Regeneration_test_total) %>%
  filter(!is.na(Regeneration_test_total)) %>%
  mutate(
    TotalTrials = 3,
    InactiveCount = TotalTrials - Regeneration_test_total,
    ActiveProportion = Regeneration_test_total / TotalTrials
  )

# Since we only have one condition, we'll compare Head vs Tail directly
memory_test_ttest <- t.test(
  ActiveProportion ~ BodyPart, 
  data = Memory_test_data, 
  var.equal = TRUE
)

# Print test result
print("Memory Test: Head vs Tail comparison")
print(memory_test_ttest)

# Calculate effect size (Cohen's d)
head_mean <- mean(Memory_test_data$ActiveProportion[Memory_test_data$BodyPart == "Head"])
tail_mean <- mean(Memory_test_data$ActiveProportion[Memory_test_data$BodyPart == "Tail"])
pooled_sd <- sqrt(((sum((Memory_test_data$ActiveProportion[Memory_test_data$BodyPart == "Head"] - head_mean)^2) + 
                      sum((Memory_test_data$ActiveProportion[Memory_test_data$BodyPart == "Tail"] - tail_mean)^2)) / 
                     (nrow(Memory_test_data) - 2)))
memory_cohens_d <- (head_mean - tail_mean) / pooled_sd
print(paste("Memory Test Cohen's d:", memory_cohens_d))


#=================================================================
# PART 6: COMPARING REGENERATED PARTS TO ORIGINAL BASELINE SCORES
#=================================================================

# First, identify only the subjects that were cut
cut_subjects <- Exp4_full_data %>%
  filter(grepl("^[0-9]+$", as.character(Subject)) & `Selected to be cut` == "Y") %>%
  mutate(Subject = as.character(Subject))  # Ensure Subject is character type

# Create a baseline scores dataset for ONLY the cut subjects
cut_subjects_baseline <- cut_subjects %>%
  mutate(
    BaselineProportion = baseline_active_arm_entries / 6,
    Subject = as.character(Subject)  # Ensure Subject is character type for joining
  )

# Print the subjects that were cut
print("Subjects that were cut:")
print(cut_subjects$Subject)
print(paste("Number of cut subjects:", nrow(cut_subjects)))

# Create a dataset linking original subjects to their regenerated parts
Original_to_regen <- Exp4_full_data %>%
  filter(grepl("[HT]$", Subject)) %>%
  mutate(
    OrigSubject = str_remove(Subject, "[HT]$"),  # Keep as character for joining
    BodyPart = str_extract(Subject, "[HT]$"),
    BodyPart = factor(BodyPart, levels = c("H", "T"), labels = c("Head", "Tail"))
  ) %>%
  # Only join with the cut subjects' baseline data
  inner_join(
    cut_subjects_baseline %>% 
      select(Subject, BaselineProportion),
    by = c("OrigSubject" = "Subject")
  )

# Print the regenerated parts and their original subjects
print("Regenerated parts and their original subjects:")
print(Original_to_regen %>% select(Subject, OrigSubject, BodyPart))
print(paste("Number of regenerated parts:", nrow(Original_to_regen)))

# Create datasets for the regeneration and reinstatement phases - UPDATED
Memory_test_comparison <- Original_to_regen %>%
  mutate(
    RegeneratedProportion = Regeneration_test_total / 3,
    Change = RegeneratedProportion - BaselineProportion
  )

# Create datasets for the reinstatement phases - UPDATED
Reinstatement_comparison <- Original_to_regen %>%
  mutate(
    ReinstatedProportion = Reinstatment_test_total / 3,
    Change = ReinstatedProportion - BaselineProportion
  )

# Verify we're using the same baseline values for both comparisons
print("Regeneration analysis - Mean baseline value:")
print(mean(Memory_test_comparison$BaselineProportion, na.rm = TRUE))

print("Reinstatement analysis - Mean baseline value:")
print(mean(Reinstatement_comparison$BaselineProportion, na.rm = TRUE))

# Prepare data for plotting with consistent format
# REGENERATION PLOT DATA
regen_plot_data <- rbind(
  # Original baseline data
  cut_subjects_baseline %>%
    mutate(BodyStatus = "Original") %>%
    select(BodyStatus, Proportion = BaselineProportion),
  
  # Head data
  Memory_test_comparison %>%
    filter(BodyPart == "Head") %>%
    mutate(BodyStatus = "Head") %>%
    select(BodyStatus, Proportion = RegeneratedProportion),
  
  # Tail data
  Memory_test_comparison %>%
    filter(BodyPart == "Tail") %>%
    mutate(BodyStatus = "Tail") %>%
    select(BodyStatus, Proportion = RegeneratedProportion)
)

# REINSTATEMENT PLOT DATA - using SAME original baseline data
reinstate_plot_data <- rbind(
  # Original baseline data - SAME AS ABOVE
  cut_subjects_baseline %>%
    mutate(BodyStatus = "Original") %>%
    select(BodyStatus, Proportion = BaselineProportion),
  
  # Head data
  Reinstatement_comparison %>%
    filter(BodyPart == "Head") %>%
    mutate(BodyStatus = "Head") %>%
    select(BodyStatus, Proportion = ReinstatedProportion),
  
  # Tail data
  Reinstatement_comparison %>%
    filter(BodyPart == "Tail") %>%
    mutate(BodyStatus = "Tail") %>%
    select(BodyStatus, Proportion = ReinstatedProportion)
)

# Set factor levels for proper ordering
regen_plot_data$BodyStatus <- factor(regen_plot_data$BodyStatus, 
                                     levels = c("Original", "Head", "Tail"))
reinstate_plot_data$BodyStatus <- factor(reinstate_plot_data$BodyStatus, 
                                         levels = c("Original", "Head", "Tail"))

# Calculate summary statistics for plotting
regen_summary <- regen_plot_data %>%
  group_by(BodyStatus) %>%
  summarise(
    mean_prop = mean(Proportion, na.rm = TRUE),
    se = sd(Proportion, na.rm = TRUE) / sqrt(sum(!is.na(Proportion))),
    .groups = 'drop'
  )

reinstate_summary <- reinstate_plot_data %>%
  group_by(BodyStatus) %>%
  summarise(
    mean_prop = mean(Proportion, na.rm = TRUE),
    se = sd(Proportion, na.rm = TRUE) / sqrt(sum(!is.na(Proportion))),
    .groups = 'drop'
  )




# Create baseline data frame
Baseline_scores <- cut_subjects_baseline %>%
  select(Subject, BaselineProportion)

# Create datasets for paired analysis
# For regeneration test
Regeneration_vs_baseline <- Memory_test_comparison %>%
  select(Subject, OrigSubject, BodyPart, RegeneratedProportion) %>%
  rename(ActiveProportion = RegeneratedProportion)

# For reinstatement test
Reinstatement_vs_baseline <- Reinstatement_comparison %>%
  select(Subject, OrigSubject, BodyPart, ReinstatedProportion) %>%
  rename(ActiveProportion = ReinstatedProportion)

# Reshape data for paired analysis - separate for Head and Tail
# For Regeneration
head_regeneration_paired <- Regeneration_vs_baseline %>%
  filter(BodyPart == "Head") %>%
  inner_join(Baseline_scores, by = c("OrigSubject" = "Subject")) %>%
  select(OrigSubject, ActiveProportion, BaselineProportion)

tail_regeneration_paired <- Regeneration_vs_baseline %>%
  filter(BodyPart == "Tail") %>%
  inner_join(Baseline_scores, by = c("OrigSubject" = "Subject")) %>%
  select(OrigSubject, ActiveProportion, BaselineProportion)

# For Reinstatement
head_reinstatement_paired <- Reinstatement_vs_baseline %>%
  filter(BodyPart == "Head") %>%
  inner_join(Baseline_scores, by = c("OrigSubject" = "Subject")) %>%
  select(OrigSubject, ActiveProportion, BaselineProportion)

tail_reinstatement_paired <- Reinstatement_vs_baseline %>%
  filter(BodyPart == "Tail") %>%
  inner_join(Baseline_scores, by = c("OrigSubject" = "Subject")) %>%
  select(OrigSubject, ActiveProportion, BaselineProportion)

# Perform paired t-tests
# For Regeneration
head_regeneration_test <- head_regeneration_paired %>%
  with(t.test(ActiveProportion, BaselineProportion, paired = TRUE))

tail_regeneration_test <- tail_regeneration_paired %>%
  with(t.test(ActiveProportion, BaselineProportion, paired = TRUE))

# For Reinstatement
head_reinstatement_test <- head_reinstatement_paired %>%
  with(t.test(ActiveProportion, BaselineProportion, paired = TRUE))

tail_reinstatement_test <- tail_reinstatement_paired %>%
  with(t.test(ActiveProportion, BaselineProportion, paired = TRUE))

# Print results
print("PAIRED T-TESTS RESULTS:")

print("Head Regeneration vs Original Baseline:")
print(head_regeneration_test)

print("Tail Regeneration vs Original Baseline:")
print(tail_regeneration_test)

print("Head Reinstatement vs Original Baseline:")
print(head_reinstatement_test)

print("Tail Reinstatement vs Original Baseline:")
print(tail_reinstatement_test)

# Calculate effect sizes (Cohen's d)
# For Head Regeneration
head_regeneration_d <- head_regeneration_paired %>%
  summarise(
    mean_diff = mean(ActiveProportion - BaselineProportion, na.rm = TRUE),
    sd_diff = sd(ActiveProportion - BaselineProportion, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

# For Tail Regeneration
tail_regeneration_d <- tail_regeneration_paired %>%
  summarise(
    mean_diff = mean(ActiveProportion - BaselineProportion, na.rm = TRUE),
    sd_diff = sd(ActiveProportion - BaselineProportion, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

# For Head Reinstatement
head_reinstatement_d <- head_reinstatement_paired %>%
  summarise(
    mean_diff = mean(ActiveProportion - BaselineProportion, na.rm = TRUE),
    sd_diff = sd(ActiveProportion - BaselineProportion, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

# For Tail Reinstatement
tail_reinstatement_d <- tail_reinstatement_paired %>%
  summarise(
    mean_diff = mean(ActiveProportion - BaselineProportion, na.rm = TRUE),
    sd_diff = sd(ActiveProportion - BaselineProportion, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

print("Effect sizes:")
print("Head Regeneration (Cohen's d):")
print(head_regeneration_d)

print("Tail Regeneration (Cohen's d):")
print(tail_regeneration_d)

print("Head Reinstatement (Cohen's d):")
print(head_reinstatement_d)

print("Tail Reinstatement (Cohen's d):")
print(tail_reinstatement_d)


# Create a function to format t-test results in APA style
format_apa_results <- function(test_name, t_test_result, paired_data) {
  # Extract t-value and degrees of freedom
  t_value <- t_test_result$statistic
  df <- t_test_result$parameter
  
  # Calculate Cohen's d for paired data
  n <- nrow(paired_data)
  mean_diff <- mean(paired_data$ActiveProportion - paired_data$BaselineProportion, na.rm = TRUE)
  sd_diff <- sd(paired_data$ActiveProportion - paired_data$BaselineProportion, na.rm = TRUE)
  cohens_d <- mean_diff / sd_diff
  
  # Format p-value according to APA style
  p_value <- t_test_result$p.value
  p_formatted <- ifelse(p_value < .001, "< .001", 
                        paste0("= ", gsub("0\\.", ".", round(p_value, 3))))
  
  # Create results dataframe
  result_df <- data.frame(
    contrast = test_name,
    t_statistic = t_value,
    df = df,
    cohens_d = abs(cohens_d),
    p.value = p_value,
    apa_result = paste0("*t*(", round(df, 2), ") = ", round(t_value, 2), 
                        ", *d* = ", round(abs(cohens_d), 2), 
                        ", *p* ", p_formatted)
  )
  
  return(result_df)
}

# Create APA-formatted results for each test
head_regen_apa <- format_apa_results(
  "Head Regeneration vs Baseline", 
  head_regeneration_test, 
  head_regeneration_paired
)

tail_regen_apa <- format_apa_results(
  "Tail Regeneration vs Baseline", 
  tail_regeneration_test, 
  tail_regeneration_paired
)

head_reinstate_apa <- format_apa_results(
  "Head Reinstatement vs Baseline", 
  head_reinstatement_test, 
  head_reinstatement_paired
)

tail_reinstate_apa <- format_apa_results(
  "Tail Reinstatement vs Baseline", 
  tail_reinstatement_test, 
  tail_reinstatement_paired
)

# Combine all results into a single dataframe
regen_reinstate_results <- rbind(
  head_regen_apa,
  tail_regen_apa,
  head_reinstate_apa,
  tail_reinstate_apa
)

# Print the combined results
print(regen_reinstate_results)

# For easy referencing in inline code within Quarto
# These variables can be used with `r head_regen_result` etc. in your Quarto document
head_regen_result <- head_regen_apa$apa_result
tail_regen_result <- tail_regen_apa$apa_result
head_reinstate_result <- head_reinstate_apa$apa_result
tail_reinstate_result <- tail_reinstate_apa$apa_result

# Also save detailed result objects that you can reference for specific values
head_regen_details <- head_regen_apa
tail_regen_details <- tail_regen_apa
head_reinstate_details <- head_reinstate_apa
tail_reinstate_details <- tail_reinstate_apa


#=================================================================
# PART 8: VISUALIZATIONS FOR PAIRED COMPARISONS
#=================================================================

# Create visualization for head regeneration
head_regeneration_viz <- head_regeneration_paired %>%
  pivot_longer(
    cols = c(ActiveProportion, BaselineProportion),
    names_to = "Measurement",
    values_to = "Proportion"
  ) %>%
  mutate(
    Measurement = factor(Measurement, 
                         levels = c("BaselineProportion", "ActiveProportion"),
                         labels = c("Original", "Regenerated Head"))
  ) %>%
  ggplot(aes(x = Measurement, y = Proportion)) +
  # Add individual subject lines
  geom_line(aes(group = OrigSubject), alpha = 0.3, color = "#FF8C00") +
  # Add points
  geom_point(size = 3, color = "#FF8C00") +
  # Add means and error bars
  stat_summary(
    fun = mean, geom = "point", size = 5, shape = 18, color = "#FF8C00"
  ) +
  stat_summary(
    fun.data = function(x) {
      return(c(y = mean(x), ymin = mean(x) - sd(x)/sqrt(length(x)), 
               ymax = mean(x) + sd(x)/sqrt(length(x))))
    },
    geom = "errorbar", width = 0.2, linewidth = 1, color = "#FF8C00"
  ) +
  # Add p-value annotations
  annotate("text", x = 1.5, y = 0.95, 
           label = paste0("p = ", format(head_regeneration_test$p.value, digits = 3)), 
           color = "#FF8C00", fontface = "bold", size = 3.5) +
  # Customize aesthetics
  labs(
    title = "Head Regeneration vs Original Baseline",
    subtitle = "Lines connect individual subjects across measurements",
    y = "Proportion of Active Arm Entries",
    x = "Measurement"
  ) +
  ylim(0, 1) +
  consistent_theme()

# Create visualization for tail regeneration
tail_regeneration_viz <- tail_regeneration_paired %>%
  pivot_longer(
    cols = c(ActiveProportion, BaselineProportion),
    names_to = "Measurement",
    values_to = "Proportion"
  ) %>%
  mutate(
    Measurement = factor(Measurement, 
                         levels = c("BaselineProportion", "ActiveProportion"),
                         labels = c("Original", "Regenerated Tail"))
  ) %>%
  ggplot(aes(x = Measurement, y = Proportion)) +
  # Add individual subject lines
  geom_line(aes(group = OrigSubject), alpha = 0.3, color = "#159090") +
  # Add points
  geom_point(size = 3, color = "#159090") +
  # Add means and error bars
  stat_summary(
    fun = mean, geom = "point", size = 5, shape = 18, color = "#159090"
  ) +
  stat_summary(
    fun.data = function(x) {
      return(c(y = mean(x), ymin = mean(x) - sd(x)/sqrt(length(x)), 
               ymax = mean(x) + sd(x)/sqrt(length(x))))
    },
    geom = "errorbar", width = 0.2, linewidth = 1, color = "#159090"
  ) +
  # Add p-value annotations
  annotate("text", x = 1.5, y = 0.95, 
           label = paste0("p = ", format(tail_regeneration_test$p.value, digits = 3)), 
           color = "#159090", fontface = "bold", size = 3.5) +
  # Customize aesthetics
  labs(
    title = "Tail Regeneration vs Original Baseline",
    subtitle = "Lines connect individual subjects across measurements",
    y = "Proportion of Active Arm Entries",
    x = "Measurement"
  ) +
  ylim(0, 1) +
  consistent_theme()

# Create visualization for head reinstatement
head_reinstatement_viz <- head_reinstatement_paired %>%
  pivot_longer(
    cols = c(ActiveProportion, BaselineProportion),
    names_to = "Measurement",
    values_to = "Proportion"
  ) %>%
  mutate(
    Measurement = factor(Measurement, 
                         levels = c("BaselineProportion", "ActiveProportion"),
                         labels = c("Original", "Reinstated Head"))
  ) %>%
  ggplot(aes(x = Measurement, y = Proportion)) +
  # Add individual subject lines
  geom_line(aes(group = OrigSubject), alpha = 0.3, color = "#FF8C00") +
  # Add points
  geom_point(size = 3, color = "#FF8C00") +
  # Add means and error bars
  stat_summary(
    fun = mean, geom = "point", size = 5, shape = 18, color = "#FF8C00"
  ) +
  stat_summary(
    fun.data = function(x) {
      return(c(y = mean(x), ymin = mean(x) - sd(x)/sqrt(length(x)), 
               ymax = mean(x) + sd(x)/sqrt(length(x))))
    },
    geom = "errorbar", width = 0.2, linewidth = 1, color = "#FF8C00"
  ) +
  # Add p-value annotations
  annotate("text", x = 1.5, y = 0.95, 
           label = paste0("p = ", format(head_reinstatement_test$p.value, digits = 3)), 
           color = "#FF8C00", fontface = "bold", size = 3.5) +
  # Customize aesthetics
  labs(
    title = "Head Reinstatement vs Original Baseline",
    subtitle = "Lines connect individual subjects across measurements",
    y = "Proportion of Active Arm Entries",
    x = "Measurement"
  ) +
  ylim(0, 1) +
  consistent_theme()

# Create visualization for tail reinstatement
tail_reinstatement_viz <- tail_reinstatement_paired %>%
  pivot_longer(
    cols = c(ActiveProportion, BaselineProportion),
    names_to = "Measurement",
    values_to = "Proportion"
  ) %>%
  mutate(
    Measurement = factor(Measurement, 
                         levels = c("BaselineProportion", "ActiveProportion"),
                         labels = c("Original", "Reinstated Tail"))
  ) %>%
  ggplot(aes(x = Measurement, y = Proportion)) +
  # Add individual subject lines
  geom_line(aes(group = OrigSubject), alpha = 0.3, color = "#159090") +
  # Add points
  geom_point(size = 3, color = "#159090") +
  # Add means and error bars
  stat_summary(
    fun = mean, geom = "point", size = 5, shape = 18, color = "#159090"
  ) +
  stat_summary(
    fun.data = function(x) {
      return(c(y = mean(x), ymin = mean(x) - sd(x)/sqrt(length(x)), 
               ymax = mean(x) + sd(x)/sqrt(length(x))))
    },
    geom = "errorbar", width = 0.2, linewidth = 1, color = "#159090"
  ) +
  # Add p-value annotations
  annotate("text", x = 1.5, y = 0.95, 
           label = paste0("p = ", format(tail_reinstatement_test$p.value, digits = 3)), 
           color = "#159090", fontface = "bold", size = 3.5) +
  # Customize aesthetics
  labs(
    title = "Tail Reinstatement vs Original Baseline",
    subtitle = "Lines connect individual subjects across measurements",
    y = "Proportion of Active Arm Entries",
    x = "Measurement"
  ) +
  ylim(0, 1) +
  consistent_theme()



#=================================================================
# PART 9: FINAL COMBINED VISUALIZATION
#=================================================================

# Create a dataset for all four conditions
combined_data <- bind_rows(
  # Head Regeneration
  head_regeneration_paired %>%
    mutate(
      BodyPart = "Head",
      Phase = "Regeneration",
      Original = BaselineProportion,
      Regenerated = ActiveProportion
    ) %>%
    select(OrigSubject, BodyPart, Phase, Original, Regenerated),
  
  # Tail Regeneration
  tail_regeneration_paired %>%
    mutate(
      BodyPart = "Tail",
      Phase = "Regeneration",
      Original = BaselineProportion,
      Regenerated = ActiveProportion
    ) %>%
    select(OrigSubject, BodyPart, Phase, Original, Regenerated),
  
  # Head Reinstatement
  head_reinstatement_paired %>%
    mutate(
      BodyPart = "Head",
      Phase = "Reinstatement",
      Original = BaselineProportion,
      Regenerated = ActiveProportion
    ) %>%
    select(OrigSubject, BodyPart, Phase, Original, Regenerated),
  
  # Tail Reinstatement
  tail_reinstatement_paired %>%
    mutate(
      BodyPart = "Tail",
      Phase = "Reinstatement",
      Original = BaselineProportion,
      Regenerated = ActiveProportion
    ) %>%
    select(OrigSubject, BodyPart, Phase, Original, Regenerated)
)

# Reshape to long format for plotting
combined_long <- combined_data %>%
  pivot_longer(
    cols = c(Original, Regenerated),
    names_to = "Measurement",
    values_to = "Proportion"
  ) %>%
  mutate(
    BodyPart = factor(BodyPart, levels = c("Head", "Tail")),
    Phase = factor(Phase, levels = c("Regeneration", "Reinstatement")),
    Measurement = factor(Measurement, levels = c("Original", "Regenerated"))
  )


# Create a more visually appealing theme that matches the desired output
publication_theme <- function() {
  theme_classic() +
    theme(
      text = element_text(family = "Times New Roman", size = 14),
      axis.title = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(margin = margin(r = 20)),
      axis.text = element_text(size = 116, color = "black"),
      legend.position = c(0.85, 0.85),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      axis.line = element_line(color = "black", linewidth = 0.8),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "lightgray", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
}


# Create an improved version of the regeneration plot
# Create an improved version of the regeneration plot
improved_regeneration_plot <- ggplot(regen_summary, aes(x = BodyStatus, y = mean_prop)) +
  # Add error bars
  geom_errorbar(
    aes(ymin = mean_prop - se, ymax = mean_prop + se, color = BodyStatus),
    width = 0.2, 
    linewidth = 0.7,
    position = position_dodge(0.9)
  ) +
  # Add points with custom shapes - ADD FILL AESTHETIC
  geom_point(
    aes(shape = BodyStatus, color = BodyStatus, fill = BodyStatus),  # Added fill aesthetic
    size = 5,
    stroke = 0.7
  ) +
  # Set shapes manually
  scale_shape_manual(
    name = "Body Status",
    values = c(
      "Original" = 21,  # Circle
      "Head" = 24,      # Triangle
      "Tail" = 25       # CHANGED to 24 (filled triangle down) instead of 25
    ),
    labels = c(
      "Original" = "Intact Baseline",
      "Head" = "Head Regenerates",
      "Tail" = "Tail Regenerates"
    )
  ) +
  # Set colors manually
  scale_color_manual(
    name = "Body Status",
    values = c(
      "Original" = "#FF8C00", 
      "Head" = "#FF8C00",
      "Tail" = "#FF8C00"
    ),
    labels = c(
      "Original" = "Intact Baseline",
      "Head" = "Head Regenerates",
      "Tail" = "Tail Regenerates"
    )
  ) +
  # Add fill scale with the same colors
  scale_fill_manual(
    name = "Body Status",
    values = c(
      "Original" = "#FF8C00", 
      "Head" = "#FF8C00",
      "Tail" = "#FF8C00"
    ),
    labels = c(
      "Original" = "Intact Baseline",
      "Head" = "Head Regenerates",
      "Tail" = "Tail Regenerates"
    )
  ) +
  # Labels
  labs(
    title = "Regeneration Compared to Baseine",
    y = "Proportion of Active Arm Entries",
    x = "Body Section"
  ) +
  # Y-axis limits
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  # Apply publication theme
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 20),
    axis.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  # Combine duplicate legends (color, fill, and shape) into one
  guides(
    color = guide_legend(title = "Body Section", override.aes = list(size = 5)),
    fill = "none",  # Hide duplicate legend
    shape = "none"  # Hide duplicate legend
  )

# Create an improved version of the reinstatement plot
improved_reinstatement_plot <- ggplot(reinstate_summary, aes(x = BodyStatus, y = mean_prop)) +
  # Add error bars
  geom_errorbar(
    aes(ymin = mean_prop - se, ymax = mean_prop + se, color = BodyStatus),
    width = 0.2, 
    linewidth = 0.7,
    position = position_dodge(0.9)
  ) +
  # Add points with custom shapes - ADD FILL AESTHETIC
  geom_point(
    aes(shape = BodyStatus, color = BodyStatus, fill = BodyStatus),  # Added fill aesthetic
    size = 5,
    stroke = 0.7
  ) +
  # Set shapes manually
  scale_shape_manual(
    name = "Body Section",
    values = c(
      "Original" = 21,  # Circle
      "Head" = 24,      # Triangle
      "Tail" = 25       # CHANGED to 24 (filled triangle down) instead of 25
    ),
    labels = c(
      "Original" = "Intact Baseline",
      "Head" = "Head Regenerates",
      "Tail" = "Tail Regenerates"
    )
  ) +
  # Set colors manually
  scale_color_manual(
    name = "Body Section",
    values = c(
      "Original" = "#FF8C00", 
      "Head" = "#FF8C00",
      "Tail" = "#FF8C00"
    ),
    labels = c(
      "Original" = "Intact Baseline",
      "Head" = "Head Regenerates",
      "Tail" = "Tail Regenerates"
    )
  ) +
  # Add fill scale with the same colors
  scale_fill_manual(
    name = "Body Section",
    values = c(
      "Original" = "#FF8C00", 
      "Head" = "#FF8C00",
      "Tail" = "#FF8C00"
    ),
    labels = c(
      "Original" = "Intact Baseline",
      "Head" = "Head Regenerates",
      "Tail" = "Tail Regenerates"
    )
  ) +
  # Labels
  labs(
    title = "Reinstatement Active Arm Preference",
    y = "Proportion of Active Arm Entries",
    x = "Body Section"
  ) +
  # Y-axis limits
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  # Apply publication theme
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 20),
    axis.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  # Combine duplicate legends (color, fill, and shape) into one
  guides(
    color = guide_legend(title = "Body Status", override.aes = list(size = 5)),
    fill = "none",  # Hide duplicate legend
    shape = "none"  # Hide duplicate legend
  )

# Combine the plots using patchwork
combined_plots <- improved_regeneration_plot + improved_reinstatement_plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Memory Retention in Regenerated Planarians",
    theme = theme(
      plot.title = element_text(family = "Times New Roman", size = 20, face = "bold", hjust = 0.5),
      text = element_text(family = "Times New Roman", size = 18),
      axis.text = element_text(size = 15, color = "black")
    )
  )

# Print the plots
print(improved_regeneration_plot)
print(improved_reinstatement_plot)







###### the retention and reinstatement plots for regenerates

# Combine the regeneration and reinstatement plots

# Remvoe legend from one of the grapghs
reinstatement_plot_no_legend <- improved_reinstatement_plot +
  theme(legend.position = "none")

# Combine the plots, keeping only the regeneration_plot legend
Regen_combined_figure <- (improved_regeneration_plot + reinstatement_plot_no_legend) +
  plot_layout(
    guides = "collect",
    widths = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(face = "bold", size = 16, family = "Times New Roman"),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
  )

print(Regen_combined_figure)

#=======================================================================
# PART 10: Combined visualisation of both conditioning and regeneration
#=======================================================================

#### removing the uneeded legends
learning_plot_no_legend <- learning_plot + 
  theme(legend.position = "none")

intact_plot_no_legend <- intact_plot + 
  theme(legend.position = "none")

regeneration_plot_no_legend <- improved_regeneration_plot + 
  theme(legend.position = "none")

# Apply theme to each individual plot
learning_plot_no_legend <- learning_plot_no_legend +
  theme(
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold", family = "Times New Roman") 
  )

intact_plot_no_legend <- intact_plot_no_legend +
  theme(
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold", family = "Times New Roman") 
  )

regeneration_plot_no_legend <- regeneration_plot_no_legend +
  theme(
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold", family = "Times New Roman") 
  )

improved_reinstatement_plot <- improved_reinstatement_plot +
  theme(
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold", family = "Times New Roman") 
  )

##### Ordering the legend

reinstatement_plot_ordered <- improved_reinstatement_plot +
  theme(
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold", family = "Times New Roman")) +
  guides(
    # Set up the shape legend
    shape = guide_legend(
      order = 2, 
      title = "Legend",
      override.aes = list(stroke = 0.7, size = 5, color = "black", fill = "#FF8C00")
    ),
    color = "none"  # Hide separate color legend since it's redundant
    )

combined_figure <- (learning_plot_no_legend + intact_plot_no_legend) / 
  (regeneration_plot_no_legend + reinstatement_plot_ordered) +
  plot_layout(
    guides = "collect",
    widths = c(1, 1),
    heights = c(0.7, 0.7)
  ) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(face = "bold", size = 20, family = "Times New Roman"),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
  )


print(combined_figure)


# Save the panel
ggsave("Exp4_combined_figure.png", combined_figure, 
       width = 18, height = 16, dpi = 300) 


#=======================================================================
# PART 10: Extra descriptive information
#=======================================================================


# Looking for the spread of initial baseline preferences. 

Prefered_arm_count <- Exp4_full_data %>% filter(Subject %in% (1:42)) %>% count(Prefered_arm)

Active_arm_count <- Exp4_full_data %>% filter(Subject %in% (1:42)) %>% count(`Active arm`)

Left_active_arm_count <- Active_arm_count %>% filter(`Active arm` == "L") %>% pull(n)
Right_active_arm_count <- Active_arm_count %>% filter(`Active arm` == "R") %>% pull(n)


# Finding endpoint group mean

exp4_endpoint_mean <- Exp4_full_data_long %>% filter(Time %in% "Endpoint") %>% summarize(mean_entries = mean(ActiveArmProportion))

