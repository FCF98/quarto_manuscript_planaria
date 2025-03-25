# Planaria Y-Maze Analysis for Experiment 4
# Load required libraries
library(readr)
library(tidyverse)
library(lme4)
library(emmeans)
library(grid)
library(gridExtra)
library(car)
library(ggplot2)
library(readxl)
library(patchwork) # For combining plots

# Read in the data
Exp4_full_data <- read_csv("ExperimentFourFullData.csv")

#=================================================================
# Creating theme for all plots
#=================================================================

consistent_theme <- function() {
  theme_classic() +
    theme(
      text = element_text(family = "Times New Roman", size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(margin = margin(r = 20)),
      axis.text = element_text(size = 12, color = "black"),
      legend.position = "top",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
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

# Filter for original intact subjects (numeric subject IDs)
Orig_subjects <- Exp4_data %>%
  filter(is.numeric(Subject))

# Restructure data for intact analysis (comparing baseline vs endpoint)
Exp4_data_long <- Orig_subjects %>%
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

# For Drug condition only analysis, we'll use a paired t-test instead of GLMM
# This is simpler for a single treatment group
paired_ttest_model <- t.test(
  Exp4_data_long %>% filter(Time == "Baseline") %>% pull(ActiveArmProportion),
  Exp4_data_long %>% filter(Time == "Endpoint") %>% pull(ActiveArmProportion),
  paired = TRUE
)

# Print t-test results
print(paired_ttest_model)

# Calculate effect size (Cohen's d for paired data)
t_value <- paired_ttest_model$statistic
n <- length(Exp4_data_long %>% filter(Time == "Baseline") %>% pull(ActiveArmProportion))
cohen_d <- t_value / sqrt(n)
print(paste("Cohen's d:", cohen_d))

# Visualization for intact subjects
intact_plot <- Exp4_data_long %>%
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
    y_position = max(Exp4_data_long$ActiveArmProportion, na.rm = TRUE) + 0.1,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  labs(
    title = "Active Arm Choices Before and After Conditioning",
    y = "Proportion of Active Arm Choices",
    x = "Time"
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
Exp4_data_long_days <- Orig_subjects %>%
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
                                  "Conditioning_day5")),
    Subject = factor(Subject),
    Proportion = ActiveArmChoices / 3 # 3 trials per day
  )

# Create the learning curve visualization
learning_plot <- Exp4_data_long_days %>%
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
    y = "Proportion of Active Arm Choices",
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
# PART 3: REGENERATION ANALYSIS
#=================================================================

# Filter for regenerated subjects (subjects with H or T suffix)
Regen_subjects <- Exp4_data %>%
  filter(grepl("[HT]$", Subject)) %>%
  mutate(
    # Extract original subject number and body part
    OrigSubject = as.numeric(str_remove(Subject, "[HT]$")),
    BodyPart = str_extract(Subject, "[HT]$"),
    BodyPart = factor(BodyPart, levels = c("H", "T"), labels = c("Head", "Tail")),
    Subject = factor(Subject)
  )

# Create a dataset for the memory test phase
Memory_test_data <- Regen_subjects %>%
  select(Subject, OrigSubject, BodyPart, T1, T2, T3, Test_total) %>%
  filter(!is.na(Test_total)) %>%
  mutate(
    TotalTrials = 3,
    InactiveCount = TotalTrials - Test_total,
    ActiveProportion = Test_total / TotalTrials
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

# Visualization of memory test results
memory_plot <- Memory_test_data %>%
  group_by(BodyPart) %>%
  summarise(
    mean_prop = mean(ActiveProportion, na.rm = TRUE),
    se = sd(ActiveProportion, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = BodyPart, y = mean_prop)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6, fill = "#FF8C00") +
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.8) +
  # Add significance annotation based on t-test result
  geom_signif(
    annotations = ifelse(memory_test_ttest$p.value < 0.05, "*", "ns"),
    xmin = 1, xmax = 2,
    y_position = max(Memory_test_data$ActiveProportion, na.rm = TRUE) + 0.1,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  labs(
    title = "Memory Retention After Regeneration",
    y = "Proportion of Active Arm Choices",
    x = "Body Part"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  consistent_theme()

print(memory_plot)

#=================================================================
# PART 4: REINSTATEMENT ANALYSIS
#=================================================================

# Create a dataset for the reinstatement phase
Reinstatement_data <- Regen_subjects %>%
  select(Subject, OrigSubject, BodyPart, R1, R2, R3) %>%
  mutate(
    # Counting number of active arm choices in R1, R2, R3
    # Assuming R1, R2, R3 contain "L" or "R" and the active arm is defined
    ActiveArm = "R",  # Update this if it varies by subject
    Reinstatement = rowSums(across(c(R1, R2, R3), ~ . == ActiveArm), na.rm = TRUE),
    TotalTrials = 3,
    InactiveCount = TotalTrials - Reinstatement,
    ActiveProportion = Reinstatement / TotalTrials
  )

# Compare Head vs Tail for reinstatement
reinstatement_ttest <- t.test(
  ActiveProportion ~ BodyPart, 
  data = Reinstatement_data, 
  var.equal = TRUE
)

# Print test result
print("Reinstatement: Head vs Tail comparison")
print(reinstatement_ttest)

# Calculate effect size (Cohen's d)
head_mean_r <- mean(Reinstatement_data$ActiveProportion[Reinstatement_data$BodyPart == "Head"])
tail_mean_r <- mean(Reinstatement_data$ActiveProportion[Reinstatement_data$BodyPart == "Tail"])
pooled_sd_r <- sqrt(((sum((Reinstatement_data$ActiveProportion[Reinstatement_data$BodyPart == "Head"] - head_mean_r)^2) + 
                        sum((Reinstatement_data$ActiveProportion[Reinstatement_data$BodyPart == "Tail"] - tail_mean_r)^2)) / 
                       (nrow(Reinstatement_data) - 2)))
reinstatement_cohens_d <- (head_mean_r - tail_mean_r) / pooled_sd_r
print(paste("Reinstatement Cohen's d:", reinstatement_cohens_d))

# Visualization of reinstatement results
reinstatement_plot <- Reinstatement_data %>%
  group_by(BodyPart) %>%
  summarise(
    mean_prop = mean(ActiveProportion, na.rm = TRUE),
    se = sd(ActiveProportion, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = BodyPart, y = mean_prop)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6, fill = "#FF8C00") +
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.8) +
  # Add significance annotation based on t-test result
  geom_signif(
    annotations = ifelse(reinstatement_ttest$p.value < 0.05, "*", "ns"),
    xmin = 1, xmax = 2,
    y_position = max(Reinstatement_data$ActiveProportion, na.rm = TRUE) + 0.1,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  labs(
    title = "Reinstatement After Regeneration",
    y = "Proportion of Active Arm Choices",
    x = "Body Part"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  consistent_theme()

print(reinstatement_plot)

#=================================================================
# PART 5: COMPARISON OF TEST AND REINSTATEMENT
#=================================================================

# Create a combined dataset for test and reinstatement phases
Combined_regen_data <- bind_rows(
  Memory_test_data %>% 
    select(Subject, OrigSubject, BodyPart, ActiveProportion) %>%
    mutate(Phase = "Test"),
  Reinstatement_data %>% 
    select(Subject, OrigSubject, BodyPart, ActiveProportion) %>%
    mutate(Phase = "Reinstatement")
) %>%
  mutate(
    Phase = factor(Phase, levels = c("Test", "Reinstatement"))
  )

# Perform 2x2 mixed ANOVA (BodyPart as between-subjects factor, Phase as within-subjects)
# For this, we'll reshape the data to wide format for each subject
wide_data <- Combined_regen_data %>%
  pivot_wider(
    id_cols = c(Subject, OrigSubject, BodyPart),
    names_from = Phase,
    values_from = ActiveProportion
  )

# Run ANOVAs separately
aov_result <- aov(Test ~ BodyPart, data = wide_data)
print("ANOVA for Test phase:")
print(summary(aov_result))

aov_result2 <- aov(Reinstatement ~ BodyPart, data = wide_data)
print("ANOVA for Reinstatement phase:")
print(summary(aov_result2))

# For paired comparison between Test and Reinstatement
head_paired <- wide_data %>% 
  filter(BodyPart == "Head") %>%
  with(t.test(Test, Reinstatement, paired = TRUE))

tail_paired <- wide_data %>% 
  filter(BodyPart == "Tail") %>%
  with(t.test(Test, Reinstatement, paired = TRUE))

print("Head: Test vs Reinstatement comparison")
print(head_paired)

print("Tail: Test vs Reinstatement comparison")
print(tail_paired)

# Visualization for combined test and reinstatement
combined_plot <- Combined_regen_data %>%
  group_by(BodyPart, Phase) %>%
  summarise(
    mean_prop = mean(ActiveProportion, na.rm = TRUE),
    se = sd(ActiveProportion, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = Phase, y = mean_prop, fill = BodyPart)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6) +
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.8) +
  scale_fill_manual(values = c("Head" = "#FF8C00", "Tail" = "#159090")) +
  labs(
    title = "Memory Retention and Reinstatement After Regeneration",
    y = "Proportion of Active Arm Choices",
    x = "Phase"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  consistent_theme()

print(combined_plot)

#=================================================================
# PART 6: COMPARING REGENERATED PARTS TO ORIGINAL BASELINE SCORES
#=================================================================

# First, create a baseline scores dataset for original subjects
Baseline_scores <- Exp4_data %>%
  filter(is.numeric(Subject)) %>%
  mutate(
    BaselineProportion = baseline_active_arm_entries / 6
  )

# Find the original subjects that were cut and had regenerated parts
cut_subjects <- Exp4_data %>%
  filter(is.numeric(Subject) & `Selected to be cut` == "Y")

# Create a dataset linking original subjects to their regenerated parts
Original_to_regen <- Exp4_data %>%
  filter(grepl("[HT]$", Subject)) %>%
  mutate(
    OrigSubject = as.numeric(str_remove(Subject, "[HT]$")),
    BodyPart = str_extract(Subject, "[HT]$"),
    BodyPart = factor(BodyPart, levels = c("H", "T"), labels = c("Head", "Tail"))
  ) %>%
  inner_join(
    Baseline_scores %>% 
      select(Subject, BaselineProportion),
    by = c("OrigSubject" = "Subject")
  )

# Create a parallel dataset for memory test
Memory_vs_baseline <- Original_to_regen %>%
  select(Subject, OrigSubject, BodyPart, BaselineProportion, Test_total) %>%
  mutate(
    ActiveProportion = Test_total / 3
  )

# Create a parallel dataset for reinstatement
# Get the active arm for each original subject
active_arms <- Exp4_data %>%
  filter(is.numeric(Subject)) %>%
  select(Subject, `Active arm`) %>%
  rename(OrigSubject = Subject, ActiveArm = `Active arm`)

Reinstatement_vs_baseline <- Original_to_regen %>%
  select(Subject, OrigSubject, BodyPart, BaselineProportion, R1, R2, R3) %>%
  inner_join(active_arms, by = "OrigSubject") %>%
  mutate(
    # Count matches between R1, R2, R3 and ActiveArm
    R1_match = as.numeric(R1 == ActiveArm),
    R2_match = as.numeric(R2 == ActiveArm),
    R3_match = as.numeric(R3 == ActiveArm),
    Reinstatement_total = R1_match + R2_match + R3_match,
    ActiveProportion = Reinstatement_total / 3
  )

# Individual comparisons: Baseline vs Memory Test
# Create a dataset for paired comparisons with regenerated parts
Individual_comparison <- bind_rows(
  # Original baseline data
  Baseline_scores %>%
    filter(Subject %in% cut_subjects$Subject) %>%
    select(Subject, BaselineProportion) %>%
    rename(
      OrigSubject = Subject,
      ActiveProportion = BaselineProportion
    ) %>%
    mutate(
      BodyPart = "Original",
      TestPhase = "Baseline"
    ),
  
  # Memory test data for regenerated parts
  Memory_vs_baseline %>%
    select(OrigSubject, BodyPart, ActiveProportion) %>%
    mutate(TestPhase = "Memory")
) %>%
  mutate(
    OrigSubject = as.character(OrigSubject), # Ensure consistent type
    BodyPart = factor(BodyPart, levels = c("Original", "Head", "Tail")),
    TestPhase = factor(TestPhase, levels = c("Baseline", "Memory"))
  )

# Reshape data for paired analysis of Heads
paired_data_head <- Individual_comparison %>%
  filter(BodyPart %in% c("Original", "Head")) %>%
  select(OrigSubject, BodyPart, ActiveProportion) %>%
  pivot_wider(
    names_from = BodyPart,
    values_from = ActiveProportion
  )

# Perform paired t-test for heads
head_paired_test <- paired_data_head %>%
  with(t.test(Head, Original, paired = TRUE))

# Print results
print("Head vs Original Baseline (Paired t-test)")
print(head_paired_test)

# Calculate effect size (Cohen's d for paired data)
head_d <- paired_data_head %>%
  summarise(
    mean_diff = mean(Head - Original, na.rm = TRUE),
    sd_diff = sd(Head - Original, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

print("Effect size for Head vs Original (Cohen's d):")
print(head_d)

# Reshape data for paired analysis of Tails
paired_data_tail <- Individual_comparison %>%
  filter(BodyPart %in% c("Original", "Tail")) %>%
  select(OrigSubject, BodyPart, ActiveProportion) %>%
  pivot_wider(
    names_from = BodyPart,
    values_from = ActiveProportion
  )

# Perform paired t-test for tails
tail_paired_test <- paired_data_tail %>%
  with(t.test(Tail, Original, paired = TRUE))

# Print results
print("Tail vs Original Baseline (Paired t-test)")
print(tail_paired_test)

# Calculate effect size (Cohen's d for paired data)
tail_d <- paired_data_tail %>%
  summarise(
    mean_diff = mean(Tail - Original, na.rm = TRUE),
    sd_diff = sd(Tail - Original, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

print("Effect size for Tail vs Original (Cohen's d):")
print(tail_d)

# Create visualization for original vs regenerated head
head_comparison_viz <- Individual_comparison %>%
  filter(BodyPart %in% c("Original", "Head")) %>%
  ggplot(aes(x = BodyPart, y = ActiveProportion)) +
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
           label = paste0("p = ", format(head_paired_test$p.value, digits = 3)), 
           color = "#FF8C00", fontface = "bold", size = 3.5) +
  # Customize aesthetics
  labs(
    title = "Head Regeneration vs Original Baseline",
    subtitle = "Lines connect individual subjects across measurements",
    y = "Proportion of Active Arm Choices",
    x = "Body Status"
  ) +
  ylim(0, 1) +
  consistent_theme()

# Create visualization for original vs regenerated tail
tail_comparison_viz <- Individual_comparison %>%
  filter(BodyPart %in% c("Original", "Tail")) %>%
  ggplot(aes(x = BodyPart, y = ActiveProportion)) +
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
           label = paste0("p = ", format(tail_paired_test$p.value, digits = 3)), 
           color = "#159090", fontface = "bold", size = 3.5) +
  # Customize aesthetics
  labs(
    title = "Tail Regeneration vs Original Baseline",
    subtitle = "Lines connect individual subjects across measurements",
    y = "Proportion of Active Arm Choices",
    x = "Body Status"
  ) +
  ylim(0, 1) +
  consistent_theme()

# Print the visualizations
print(head_comparison_viz)
print(tail_comparison_viz)

#=================================================================
# PART 7: REINSTATEMENT VS BASELINE ANALYSIS
#=================================================================

# Now compare baseline to reinstatement
# Create a parallel dataset structure for Reinstatement
Reinstatement_comparison <- bind_rows(
  # Original baseline data
  Baseline_scores %>%
    filter(Subject %in% cut_subjects$Subject) %>%
    select(Subject, BaselineProportion) %>%
    rename(
      OrigSubject = Subject,
      ActiveProportion = BaselineProportion
    ) %>%
    mutate(
      BodyPart = "Original",
      TestPhase = "Baseline"
    ),
  
  # Reinstatement data for regenerated parts
  Reinstatement_vs_baseline %>%
    select(OrigSubject, BodyPart, ActiveProportion) %>%
    mutate(TestPhase = "Reinstatement")
) %>%
  mutate(
    OrigSubject = as.character(OrigSubject), # Ensure consistent type
    BodyPart = factor(BodyPart, levels = c("Original", "Head", "Tail")),
    TestPhase = factor(TestPhase, levels = c("Baseline", "Reinstatement"))
  )

# Reshape data for paired analysis - separate for Head and Tail
head_reinstatement_paired <- Reinstatement_comparison %>%
  filter(BodyPart %in% c("Original", "Head")) %>%
  select(OrigSubject, BodyPart, ActiveProportion) %>%
  pivot_wider(
    names_from = BodyPart,
    values_from = ActiveProportion
  )

tail_reinstatement_paired <- Reinstatement_comparison %>%
  filter(BodyPart %in% c("Original", "Tail")) %>%
  select(OrigSubject, BodyPart, ActiveProportion) %>%
  pivot_wider(
    names_from = BodyPart,
    values_from = ActiveProportion
  )

# Perform paired t-tests for Head vs Original
head_reinstatement_test <- head_reinstatement_paired %>%
  with(t.test(Head, Original, paired = TRUE))

# Perform paired t-tests for Tail vs Original
tail_reinstatement_test <- tail_reinstatement_paired %>%
  with(t.test(Tail, Original, paired = TRUE))

# Print results
print("Head Reinstatement vs Original Baseline:")
print(head_reinstatement_test)

print("Tail Reinstatement vs Original Baseline:")
print(tail_reinstatement_test)

# Calculate effect sizes (Cohen's d)
# For Head Reinstatement
head_reinstatement_d <- head_reinstatement_paired %>%
  summarise(
    mean_diff = mean(Head - Original, na.rm = TRUE),
    sd_diff = sd(Head - Original, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

# For Tail Reinstatement
tail_reinstatement_d <- tail_reinstatement_paired %>%
  summarise(
    mean_diff = mean(Tail - Original, na.rm = TRUE),
    sd_diff = sd(Tail - Original, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

print("Effect sizes for Head Reinstatement (Cohen's d):")
print(head_reinstatement_d)

print("Effect sizes for Tail Reinstatement (Cohen's d):")
print(tail_reinstatement_d)

# Create visualization for head reinstatement
head_reinstatement_viz <- Reinstatement_comparison %>%
  filter(BodyPart %in% c("Original", "Head")) %>%
  ggplot(aes(x = BodyPart, y = ActiveProportion)) +
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
    y = "Proportion of Active Arm Choices",
    x = "Body Status"
  ) +
  ylim(0, 1) +
  consistent_theme()

# Create visualization for tail reinstatement
tail_reinstatement_viz <- Reinstatement_comparison %>%
  filter(BodyPart %in% c("Original", "Tail")) %>%
  ggplot(aes(x = BodyPart, y = ActiveProportion)) +
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
    y = "Proportion of Active Arm Choices",
    x = "Body Status"
  ) +
  ylim(0, 1) +
  consistent_theme()

# Print the visualizations
print(head_reinstatement_viz)
print(tail_reinstatement_viz)

#=================================================================
# PART 8: COMBINED VISUALIZATION
#=================================================================

# Create a cleaner version of regeneration plot with only means and SE bars
regeneration_plot <- Individual_comparison %>%
  group_by(BodyPart, TestPhase) %>%
  summarise(
    Mean = mean(ActiveProportion, na.rm = TRUE),
    SE = sd(ActiveProportion, na.rm = TRUE)/sqrt(n()),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = BodyPart, y = Mean)) +
  # Add mean points with different shapes by TestPhase and BodyPart
  geom_point(
    aes(shape = interaction(TestPhase, BodyPart), fill = interaction(TestPhase, BodyPart)),
    size = 5, stroke = 0.7
  ) +
  # Add error bars
  geom_errorbar(
    aes(ymin = Mean - SE, ymax = Mean + SE, color = interaction(TestPhase, BodyPart)),
    width = 0.2, linewidth = 0.7
  ) +
  # Customize aesthetics
  scale_fill_manual(values = c(
    "Baseline.Original" = "#7f7f7f",
    "Memory.Head" = "#FF8C00", 
    "Memory.Tail" = "#159090"
  )) +
  scale_color_manual(values = c(
    "Baseline.Original" = "#7f7f7f",
    "Memory.Head" = "#FF8C00", 
    "Memory.Tail" = "#159090"
  )) +
  # Custom shape scale
  scale_shape_manual(
    name = "Test Phase",
    values = c(
      "Baseline.Original" = 21,  # Circle for baseline
      "Memory.Original" = 21,    # Circle for baseline (should not appear)
      "Memory.Head" = 24,        # Triangle pointing up for heads
      "Memory.Tail" = 25,        # Triangle pointing down for tails
      "Baseline.Head" = 24,      # Should not appear
      "Baseline.Tail" = 25       # Should not appear
    ),
    labels = c(
      "Baseline.Original" = "Intact Baseline",
      "Memory.Head" = "Head Regenerates",
      "Memory.Tail" = "Tail Regenerates"
    )
  ) +
  # Labels and titles
  labs(
    title = "Regeneration Active Arm Preference Compared to Baseline",
    y = "Proportion of Active Arm Choices",
    x = "Body Status"
  ) +
  ylim(0, 1) +
  # Use consistent theme
  consistent_theme() +
  theme(
    legend.position = "right",
  ) +
  guides(
    shape = guide_legend(title = "Legend", override.aes = list(stroke = 0.7)),
    fill = "none",  # Hide fill legend since it's redundant
    color = "none"  # Hide color legend since it's redundant
  )

# Create a cleaner version of reinstatement plot
reinstatement_plot <- Reinstatement_comparison %>%
  group_by(BodyPart, TestPhase) %>%
  summarise(
    Mean = mean(ActiveProportion, na.rm = TRUE),
    SE = sd(ActiveProportion, na.rm = TRUE)/sqrt(n()),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = BodyPart, y = Mean)) +
  # Add mean points with different shapes by TestPhase and BodyPart
  geom_point(
    aes(shape = interaction(TestPhase, BodyPart), fill = interaction(TestPhase, BodyPart)),
    size = 5, stroke = 0.7
  ) +
  # Add error bars
  geom_errorbar(
    aes(ymin = Mean - SE, ymax = Mean + SE, color = interaction(TestPhase, BodyPart)),
    width = 0.2, linewidth = 0.7
  ) +
  # Customize aesthetics
  scale_fill_manual(values = c(
    "Baseline.Original" = "#7f7f7f",
    "Reinstatement.Head" = "#FF8C00", 
    "Reinstatement.Tail" = "#159090"
  )) +
  scale_color_manual(values = c(
    "Baseline.Original" = "#7f7f7f",
    "Reinstatement.Head" = "#FF8C00", 
    "Reinstatement.Tail" = "#159090"
  )) +
  # Custom shape scale
  scale_shape_manual(
    name = "Test Phase",
    values = c(
      "Baseline.Original" = 21,  # Circle for baseline
      "Reinstatement.Original" = 21,  # Circle for baseline (should not appear)
      "Reinstatement.Head" = 24,  # Triangle pointing up for heads
      "Reinstatement.Tail" = 25,  # Triangle pointing down for tails
      "Baseline.Head" = 24,      # Should not appear
      "Baseline.Tail" = 25       # Should not appear
    ),
    labels = c(
      "Baseline.Original" = "Intact Baseline",
      "Reinstatement.Head" = "Head Regenerates",
      "Reinstatement.Tail" = "Tail Regenerates"
    )
  ) +
  # Labels and titles
  labs(
    title = "Reinstatement Active Arm Preference Compared to Baseline",
    y = "Proportion of Active Arm Choices",
    x = "Body Status"
  ) +
  ylim(0, 1) +
  # Use consistent theme
  consistent_theme() +
  theme(
    legend.position = "right",
  ) +
  guides(
    shape = guide_legend(title = "Legend", override.aes = list(stroke = 0.7)),
    fill = "none",  # Hide fill legend since it's redundant
    color = "none"  # Hide color legend since it's redundant
  )

# Remove legend from learning_plot and intact_plot
learning_plot_no_legend <- learning_plot + 
  theme(legend.position = "none")

intact_plot_no_legend <- intact_plot + 
  theme(legend.position = "none")

# Remove legend from regeneration_plot
regeneration_plot_no_legend <- regeneration_plot + 
  theme(legend.position = "none")

# Combined figure with ordered legends
combined_figure <- (learning_plot_no_legend + intact_plot_no_legend) / 
  (regeneration_plot_no_legend + reinstatement_plot) +
  plot_layout(
    guides = "collect",
    widths = c(1, 1),
    heights = c(0.7, 0.7)
  ) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(face = "bold", size = 16, family = "Times New Roman"),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
  )

print(combined_figure)

# Save with larger dimensions
ggsave(
  "Exp4_combined_figure.pdf", 
  combined_figure, 
  width = 16,
  height = 14,
  dpi = 300
)