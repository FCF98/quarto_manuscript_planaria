# Complete Planaria Y-Maze Analysis Including Regeneration
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
library(cowplot) # For get_legend and plot_grid functions

# Read in the data
Exp7_data <- read_excel("Datasets/Experiment_7_Full_Data.xlsx")


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

# Updated theme for regeneration and reinstatement plots
updated_theme <- function() {
  theme_classic() +
    theme(
      text = element_text(family = "Times New Roman", size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(margin = margin(r = 20)),
      axis.text = element_text(size = 12, color = "black"),
      legend.position = "right", # Updated to place legend on right side
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      axis.line = element_line(color = "black", linewidth = 0.8),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      # Add facet theming to match the figure
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 14, face = "bold", color = "black")
    )
}

#=================================================================
# PART 1: ANALYSIS OF INTACT SUBJECTS (PRE-REGENERATION)
#=================================================================

# Filter for original intact subjects (numeric subject IDs)
Orig_subjects <- Exp7_data %>%
  filter(!str_detect(Subject, "[HT]$"))

# Restructure data for intact analysis (comparing baseline vs endpoint)
Exp7_data_long <- Orig_subjects %>%
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
    Condition = factor(Condition),
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

# Set contrasts for intact analysis
contrasts(Exp7_data_long$Time) <- contr.sum
contrasts(Exp7_data_long$Condition) <- contr.sum

# Intact model - binomial GLMM
Intact_model <- glmer(cbind(ActiveCount, InactiveCount) ~ Condition * Time + (1|Subject), 
                      data = Exp7_data_long, 
                      family = "binomial",
                      na.action = na.omit)

# Print intact model summary
summary(Intact_model)

# Get overall effects for intact model
intact_model_output <- car::Anova(Intact_model, type = "III")
print(intact_model_output)

# Get estimated means and comparisons for both within and between groups
# Within-group comparisons (Time effects within each Condition)
intact_within_group_comparisons <- emmeans(Intact_model, 
                                           pairwise ~ Time | Condition,
                                           adjust = "bonferroni",
                                           type = "response")

# Between-group comparisons (Condition effects at each Time point)
intact_between_group_comparisons <- emmeans(Intact_model, 
                                            pairwise ~ Condition | Time,
                                            adjust = "bonferroni",
                                            type = "response")

# Print comparisons
print(summary(intact_within_group_comparisons$emmeans))
print(summary(intact_within_group_comparisons$contrasts))
print(summary(intact_between_group_comparisons$emmeans))
print(summary(intact_between_group_comparisons$contrasts))

# Visualization for intact subjects
intact_plot <- Exp7_data_long %>%
  group_by(Condition, Time) %>%
  summarise(
    mean_prop = mean(ActiveArmProportion, na.rm = TRUE),
    se = sd(ActiveArmProportion, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = Time, y = mean_prop, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6) + # Made bars narrower
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.8) +
  geom_signif(
    annotations = "*",
    xmin = 1.18, xmax = 2.18,
    y_position = 0.7,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  scale_fill_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  labs(
    title = "Active Arm Choices Before and After Conditioning",
    y = "Proportion of Active Arm Choices",
    x = "Time",
    fill = "Condition"
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
Exp7_data_long_days <- Orig_subjects %>%
  select(Subject, Condition,
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
    Condition = factor(Condition),
    Proportion = ActiveArmChoices / 3 # 3 trials per day
  )

# Create the learning curve visualization
learning_plot <- Exp7_data_long_days %>%
  group_by(Condition, TimePoint) %>%
  summarise(
    mean_prop = mean(Proportion, na.rm = TRUE),
    se = sd(Proportion, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = TimePoint, y = mean_prop, color = Condition, group = Condition)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
                width = 0.2, linewidth = 0.8) +
  scale_color_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  labs(
    title = "Average Proportion of Active Arm Choices",
    y = "Proportion of Active Arm Choices",
    x = "Day"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  consistent_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(learning_plot)

#=================================================================
# PART 3: REGENERATION ANALYSIS
#=================================================================

# Filter for regenerated subjects (subjects with H or T suffix)
Regen_subjects <- Exp7_data %>%
  filter(str_detect(Subject, "[HT]$")) %>%
  mutate(
    # Extract original subject number and body part
    OrigSubject = as.numeric(str_remove(Subject, "[HT]$")),
    BodyPart = str_extract(Subject, "[HT]$"),
    BodyPart = factor(BodyPart, levels = c("H", "T"), labels = c("Head", "Tail")),
    Condition = factor(Condition),
    Subject = factor(Subject)
  )

# Create a dataset for the memory test phase
Memory_test_data <- Regen_subjects %>%
  select(Subject, OrigSubject, BodyPart, Condition, Test) %>%
  filter(!is.na(Test)) %>%
  mutate(
    TotalTrials = 3,
    InactiveCount = TotalTrials - Test,
    ActiveProportion = Test / TotalTrials
  )

# Test if body part affects memory retention
memory_model <- glmer(cbind(Test, InactiveCount) ~ Condition * BodyPart + (1|OrigSubject), 
                      data = Memory_test_data, 
                      family = "binomial",
                      na.action = na.omit)

# Print memory model summary
summary(memory_model)

# Get overall effects for memory model
memory_model_output <- car::Anova(memory_model, type = "III")
print(memory_model_output)

# Get estimated means and comparisons
# Body part comparisons within each condition
memory_bodypart_comparisons <- emmeans(memory_model, 
                                       pairwise ~ BodyPart | Condition,
                                       adjust = "bonferroni",
                                       type = "response")

# Condition comparisons within each body part
memory_condition_comparisons <- emmeans(memory_model, 
                                        pairwise ~ Condition | BodyPart,
                                        adjust = "bonferroni",
                                        type = "response")

# Print memory comparison results
print(summary(memory_bodypart_comparisons$emmeans))
print(summary(memory_bodypart_comparisons$contrasts))
print(summary(memory_condition_comparisons$emmeans))
print(summary(memory_condition_comparisons$contrasts))

# Extract baseline performance for original intact subjects
Baseline_scores <- Orig_subjects %>%
  select(Subject, Condition, baseline_active_arm_entries) %>%
  mutate(
    Subject = as.character(Subject), # Convert to character to match with OrigSubject
    TotalBaselineTrials = 6,
    BaselineProportion = baseline_active_arm_entries / TotalBaselineTrials
  )

# Modify the Individual_comparison dataset to include all three body types with clear labels
Individual_comparison <- bind_rows(
  # Original baseline data
  Baseline_scores %>%
    select(Subject, Condition, BaselineProportion) %>%
    rename(
      OrigSubject = Subject,
      ActiveProportion = BaselineProportion
    ) %>%
    mutate(
      BodyPart = "Original",
      TestPhase = "Baseline"
    ),
  
  # Memory test data for regenerated parts
  Memory_test_data %>%
    select(OrigSubject, BodyPart, Condition, ActiveProportion) %>%
    mutate(
      OrigSubject = as.character(OrigSubject), # Convert to character to match first dataset
      TestPhase = "Memory"
    )
) %>%
  mutate(
    OrigSubject = as.character(OrigSubject), # Ensure consistent type
    BodyPart = factor(BodyPart, levels = c("Original", "Head", "Tail")),
    BodyStatus = factor(BodyPart, levels = c("Original", "Head", "Tail"), 
                        labels = c("Original", "Head", "Tail")), # For x-axis label
    TestPhase = factor(TestPhase, levels = c("Baseline", "Memory")),
    Condition = factor(Condition)
  )

# Create updated memory test visualization using points rather than bars
memory_plot <- Individual_comparison %>%
  group_by(Condition, BodyPart) %>%
  summarise(
    mean_prop = mean(ActiveProportion, na.rm = TRUE),
    se = sd(ActiveProportion, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = BodyPart, y = mean_prop, color = Condition, shape = BodyPart)) +
  # Use points instead of bars
  geom_point(size = 5) +
  # Error bars
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
                width = 0.2, linewidth = 1) +
  # Color and shape scales to match the figure
  scale_color_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  scale_shape_manual(values = c("Original" = 16, "Head" = 17, "Tail" = 25)) +
  # Update labels
  labs(
    title = "Regeneration Active Arm Preference Compared to Baseline",
    y = "Proportion of Active Arm Choices",
    x = "Body Status"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  updated_theme() +
  # Add facet by condition
  facet_wrap(~ Condition, scales = "free_x", nrow = 1) +
  # Customize the legend
  guides(
    color = guide_legend(title = "Condition"),
    shape = guide_legend(title = "Legend", 
                         override.aes = list(
                           shape = c(16, 17, 25),
                           color = rep("black", 3)
                         ))
  )

print(memory_plot)

#=================================================================
# PART 4: REINSTATEMENT ANALYSIS
#=================================================================


# Create a dataset for the reinstatement phase
Reinstatement_data <- Regen_subjects %>%
  select(Subject, OrigSubject, BodyPart, Condition, Reinstatement) %>%
  filter(!is.na(Reinstatement)) %>%
  mutate(
    TotalTrials = 3,
    InactiveCount = TotalTrials - Reinstatement,
    ActiveProportion = Reinstatement / TotalTrials
  )

# Test if body part affects reinstatement
reinstatement_model <- glmer(cbind(Reinstatement, InactiveCount) ~ Condition * BodyPart + (1|OrigSubject), 
                             data = Reinstatement_data, 
                             family = "binomial",
                             na.action = na.omit)

# Print reinstatement model summary
summary(reinstatement_model)

# Get overall effects for reinstatement model
reinstatement_model_output <- car::Anova(reinstatement_model, type = "III")
print(reinstatement_model_output)

# Get estimated means and comparisons
# Body part comparisons within each condition
reinstatement_bodypart_comparisons <- emmeans(reinstatement_model, 
                                              pairwise ~ BodyPart | Condition,
                                              adjust = "bonferroni",
                                              type = "response")

# Condition comparisons within each body part
reinstatement_condition_comparisons <- emmeans(reinstatement_model, 
                                               pairwise ~ Condition | BodyPart,
                                               adjust = "bonferroni",
                                               type = "response")

# Print reinstatement comparison results
print(summary(reinstatement_bodypart_comparisons$emmeans))
print(summary(reinstatement_bodypart_comparisons$contrasts))
print(summary(reinstatement_condition_comparisons$emmeans))
print(summary(reinstatement_condition_comparisons$contrasts))

# Create a similar dataset for reinstatement
Reinstatement_comparison <- bind_rows(
  # Original baseline data - convert Subject to character first
  Baseline_scores %>%
    select(Subject, Condition, BaselineProportion) %>%
    rename(
      OrigSubject = Subject,
      ActiveProportion = BaselineProportion
    ) %>%
    mutate(
      OrigSubject = as.character(OrigSubject), # Ensure it's character type
      BodyPart = "Original",
      TestPhase = "Baseline"
    ),
  
  # Reinstatement data for regenerated parts - make sure OrigSubject is character
  Reinstatement_data %>%
    select(OrigSubject, BodyPart, Condition, ActiveProportion) %>%
    mutate(
      OrigSubject = as.character(OrigSubject), # Convert to character
      TestPhase = "Reinstatement"
    )
) %>%
  mutate(
    BodyPart = factor(BodyPart, levels = c("Original", "Head", "Tail")),
    BodyStatus = factor(BodyPart, levels = c("Original", "Head", "Tail"), 
                        labels = c("Original", "Head", "Tail")), # For x-axis label
    TestPhase = factor(TestPhase, levels = c("Baseline", "Reinstatement")),
    Condition = factor(Condition)
  )


# Create the updated reinstatement plot
reinstatement_plot <- Reinstatement_comparison %>%
  group_by(Condition, BodyPart) %>%
  summarise(
    mean_prop = mean(ActiveProportion, na.rm = TRUE),
    se = sd(ActiveProportion, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = BodyPart, y = mean_prop, color = Condition, shape = BodyPart)) +
  # Use points instead of bars
  geom_point(size = 5) +
  # Error bars
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
                width = 0.2, linewidth = 1) +
  # Color and shape scales to match the figure
  scale_color_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  scale_shape_manual(values = c("Original" = 16, "Head" = 17, "Tail" = 25)) +
  # Update labels
  labs(
    title = "Reinstatement Active Arm Preference Compared to Baseline",
    y = "Proportion of Active Arm Choices",
    x = "Body Status"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  updated_theme() +
  # Add facet by condition
  facet_wrap(~ Condition, scales = "free_x", nrow = 1) +
  # Customize the legend
  guides(
    color = guide_legend(title = "Condition"),
    shape = guide_legend(title = "Legend", 
                         override.aes = list(
                           shape = c(16, 17, 25),
                           color = rep("black", 3)
                         ))
  )


print(reinstatement_plot)


#=================================================================
# PART 6: COMBINING VISUALIZATIONS
#=================================================================

# Create a common legend for regeneration and reinstatement plots
common_legend <- get_legend(
  memory_plot + 
    guides(
      color = guide_legend(title = "Condition"),
      shape = guide_legend(title = "Legend",
                           override.aes = list(
                             shape = c(16, 17, 25),
                             color = rep("black", 3)
                           ))
    ) +
    scale_shape_manual(
      values = c("Original" = 16, "Head" = 17, "Tail" = 25),
      labels = c("Intact Baseline", "Head Regenerates", "Tail Regenerates")
    )
)

# Remove legends from individual plots for combined figure
plot_A <- memory_plot + theme(legend.position = "none")
plot_B <- reinstatement_plot + theme(legend.position = "none")

# Combine plots with A/B labels - Option 1 using patchwork
combined_regen_figure <- (plot_A / plot_B) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(face = "bold", size = 16, family = "Times New Roman")
    )
  )

# Option 2: Add the legend to the right using cowplot
regen_figure_with_legend <- plot_grid(
  combined_regen_figure, common_legend,
  rel_widths = c(3, 1),
  ncol = 2
)

# Option 3: Alternative approach using patchwork with a single shared legend
combined_regen_figure_alt <- (plot_A + plot_B + 
                                plot_layout(ncol = 1) +
                                plot_annotation(tag_levels = 'A')) &
  theme(legend.position = "right") &
  guides(
    color = guide_legend(title = "Condition"),
    shape = guide_legend(title = "Legend",
                         labels = c("Intact Baseline", "Head Regenerates", "Tail Regenerates"))
  )

# Combine all plots into one figure (original code from your script)
combined_figure <- (learning_plot + intact_plot) / (memory_plot + reinstatement_plot) +
  plot_layout(heights = c(1.2, 1.2), guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(
    plot.tag = element_text(face = "bold", size = 16, family = "Times New Roman"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

print(combined_regen_figure_alt) # Print the regeneration/reinstatement combined figure
print(combined_figure) # Print the overall combined figure


#=================================================================
# PART 8: COMPARING REGENERATED PARTS TO ORIGINAL BASELINE SCORES
#=================================================================


# Reshape data correctly for paired analysis
paired_data <- Individual_comparison %>%
  filter(BodyPart %in% c("Original", "Head")) %>%
  select(OrigSubject, Condition, BodyPart, ActiveProportion) %>%
  pivot_wider(
    names_from = BodyPart,
    values_from = ActiveProportion
  )

# Check the structure of the paired data
print(head(paired_data, 10))

# Perform paired t-tests for each condition
# Control condition
control_paired_test <- paired_data %>%
  filter(Condition == "Control") %>%
  with(t.test(Head, Original, paired = TRUE))

# Treatment condition
treatment_paired_test <- paired_data %>%
  filter(Condition == "Treatment") %>%
  with(t.test(Head, Original, paired = TRUE))

# Print results
print("Control Condition: Head vs Original Baseline (Paired t-test)")
print(control_paired_test)

print("Treatment Condition: Head vs Original Baseline (Paired t-test)")
print(treatment_paired_test)

# Calculate effect sizes (Cohen's d for paired data)
# Control effect size
control_d <- paired_data %>%
  filter(Condition == "Control") %>%
  summarise(
    mean_diff = mean(Head - Original, na.rm = TRUE),
    sd_diff = sd(Head - Original, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

# Treatment effect size
treatment_d <- paired_data %>%
  filter(Condition == "Treatment") %>%
  summarise(
    mean_diff = mean(Head - Original, na.rm = TRUE),
    sd_diff = sd(Head - Original, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

print("Effect sizes (Cohen's d):")
print(control_d)
print(treatment_d)

# Create a paired data visualization that shows the individual connections
paired_viz <- Individual_comparison %>%
  filter(BodyPart %in% c("Original", "Head")) %>%
  ggplot(aes(x = BodyPart, y = ActiveProportion)) +
  # Add individual subject lines
  geom_line(aes(group = OrigSubject, color = Condition), alpha = 0.3) +
  # Add points
  geom_point(aes(color = Condition), size = 3) +
  # Add means and error bars
  stat_summary(
    aes(group = Condition, color = Condition),
    fun = mean, geom = "point", size = 5, shape = 18
  ) +
  stat_summary(
    aes(group = Condition, color = Condition),
    fun.data = function(x) {
      return(c(y = mean(x), ymin = mean(x) - sd(x)/sqrt(length(x)), 
               ymax = mean(x) + sd(x)/sqrt(length(x))))
    },
    geom = "errorbar", width = 0.2, linewidth = 1
  ) +
  # Add p-value annotations
  annotate("text", x = 1.5, y = 0.9, 
           label = paste0("Control: p = ", format(control_paired_test$p.value, digits = 3)), 
           color = "#159090", fontface = "bold", size = 3.5) +
  annotate("text", x = 1.5, y = 0.95, 
           label = paste0("Treatment: p = ", format(treatment_paired_test$p.value, digits = 3)), 
           color = "#FF8C00", fontface = "bold", size = 3.5) +
  # Customize aesthetics
  scale_color_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  labs(
    title = "Paired Comparison: Head Regeneration vs Original Baseline",
    subtitle = "Lines connect individual subjects across measurements",
    y = "Proportion of Active Arm Choices",
    x = "Body Status",
    color = "Condition"
  ) +
  ylim(0, 1) +
  consistent_theme() +
  facet_grid(. ~ Condition)

print(paired_viz)

# Create a clean bar chart for publication
clean_bar_chart <- paired_data %>%
  group_by(Condition) %>%
  summarise(
    Original_Mean = mean(Original, na.rm = TRUE),
    Original_SE = sd(Original, na.rm = TRUE)/sqrt(n()),
    Head_Mean = mean(Head, na.rm = TRUE),
    Head_SE = sd(Head, na.rm = TRUE)/sqrt(n()),
    .groups = 'drop'
  ) %>%
  pivot_longer(
    cols = c(Original_Mean, Head_Mean),
    names_to = "Measurement",
    values_to = "Mean"
  ) %>%
  mutate(
    BodyPart = ifelse(Measurement == "Original_Mean", "Original", "Head"),
    SE = ifelse(Measurement == "Original_Mean", Original_SE, Head_SE)
  ) %>%
  ggplot(aes(x = BodyPart, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SE, ymax = Mean + SE),
    position = position_dodge(0.9), width = 0.2
  ) +
  # Add significance indicators
  annotate("text", x = 1.5, y = 0.8, label = "***", color = "#159090", size = 6) +
  annotate("text", x = 2.5, y = 0.8, label = "***", color = "#FF8C00", size = 6) +
  scale_fill_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  labs(
    title = "Head Regeneration Memory Retention",
    y = "Proportion of Active Arm Choices",
    x = "Body Status",
    fill = "Condition"
  ) +
  ylim(0, 1) +
  consistent_theme()

print(clean_bar_chart)


##### reinstatement vs baseline analysis #####

# Create a dataset that links regenerated parts to their original baseline scores
Regen_vs_baseline <- Memory_test_data %>%
  select(Subject, OrigSubject, BodyPart, Condition, Test, ActiveProportion) %>%
  mutate(
    OrigSubject = as.character(OrigSubject) # Convert to character to match with Subject from Baseline_scores
  ) %>%
  # Join with baseline scores
  left_join(
    Baseline_scores %>% 
      select(Subject, BaselineProportion) %>%
      rename(OrigSubject = Subject),
    by = "OrigSubject"
  ) %>%
  # Create a change score (difference from baseline)
  mutate(
    ChangeFromBaseline = ActiveProportion - BaselineProportion,
    TestType = "Memory"
  )

# Add the reinstatement data with the same structure
Reinstatement_vs_baseline <- Reinstatement_data %>%
  select(Subject, OrigSubject, BodyPart, Condition, Reinstatement, ActiveProportion) %>%
  mutate(
    OrigSubject = as.character(OrigSubject) # Convert to character to match with Subject from Baseline_scores
  ) %>%
  # Join with baseline scores
  left_join(
    Baseline_scores %>% 
      select(Subject, BaselineProportion) %>%
      rename(OrigSubject = Subject),
    by = "OrigSubject"
  ) %>%
  # Create a change score (difference from baseline)
  mutate(
    ChangeFromBaseline = ActiveProportion - BaselineProportion,
    TestType = "Reinstatement"
  )

# Combine both datasets for unified analysis
Combined_vs_baseline <- bind_rows(
  Regen_vs_baseline,
  Reinstatement_vs_baseline
) %>%
  mutate(
    TestType = factor(TestType, levels = c("Memory", "Reinstatement")),
    BodyPart = factor(BodyPart, levels = c("Head", "Tail"))
  )

# Create a parallel dataset structure like Individual_comparison but for Reinstatement
Reinstatement_comparison <- bind_rows(
  # Original baseline data (same as in Individual_comparison)
  Baseline_scores %>%
    select(Subject, Condition, BaselineProportion) %>%
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
    select(OrigSubject, BodyPart, Condition, ActiveProportion) %>%
    mutate(TestPhase = "Reinstatement")
) %>%
  mutate(
    OrigSubject = as.character(OrigSubject), # Ensure consistent type
    BodyPart = factor(BodyPart, levels = c("Original", "Head", "Tail")),
    TestPhase = factor(TestPhase, levels = c("Baseline", "Reinstatement")),
    Condition = factor(Condition)
  )

# Check the structure of reinstatement comparison data
print("Structure of Reinstatement comparison data:")
print(head(Reinstatement_comparison, 10))

# Check completeness of data pairs for reinstatement
reinstatement_completeness <- Reinstatement_comparison %>%
  filter(BodyPart %in% c("Original", "Head", "Tail")) %>%
  group_by(Condition, OrigSubject) %>%
  summarise(
    has_original = any(BodyPart == "Original"),
    has_head = any(BodyPart == "Head"),
    has_tail = any(BodyPart == "Tail"),
    complete_head_pair = has_original & has_head,
    complete_tail_pair = has_original & has_tail,
    .groups = 'drop'
  )

print("Completeness of reinstatement data pairs:")
print(table(reinstatement_completeness$Condition, reinstatement_completeness$complete_head_pair))
print(table(reinstatement_completeness$Condition, reinstatement_completeness$complete_tail_pair))

# Reshape data for paired analysis - separate for Head and Tail
head_reinstatement_paired <- Reinstatement_comparison %>%
  filter(BodyPart %in% c("Original", "Head")) %>%
  select(OrigSubject, Condition, BodyPart, ActiveProportion) %>%
  pivot_wider(
    names_from = BodyPart,
    values_from = ActiveProportion
  )

tail_reinstatement_paired <- Reinstatement_comparison %>%
  filter(BodyPart %in% c("Original", "Tail")) %>%
  select(OrigSubject, Condition, BodyPart, ActiveProportion) %>%
  pivot_wider(
    names_from = BodyPart,
    values_from = ActiveProportion
  )

# Print sample of the paired data
print("Head reinstatement paired data sample:")
print(head(head_reinstatement_paired, 5))
print("Tail reinstatement paired data sample:")
print(head(tail_reinstatement_paired, 5))

# Perform paired t-tests for Head vs Original
# Control condition
control_head_reinstatement <- head_reinstatement_paired %>%
  filter(Condition == "Control") %>%
  with(t.test(Head, Original, paired = TRUE))

# Treatment condition
treatment_head_reinstatement <- head_reinstatement_paired %>%
  filter(Condition == "Treatment") %>%
  with(t.test(Head, Original, paired = TRUE))

# Perform paired t-tests for Tail vs Original
# Control condition
control_tail_reinstatement <- tail_reinstatement_paired %>%
  filter(Condition == "Control") %>%
  with(t.test(Tail, Original, paired = TRUE))

# Treatment condition
treatment_tail_reinstatement <- tail_reinstatement_paired %>%
  filter(Condition == "Treatment") %>%
  with(t.test(Tail, Original, paired = TRUE))

# Print results
print("Head Reinstatement vs Original Baseline:")
print("Control Condition:")
print(control_head_reinstatement)
print("Treatment Condition:")
print(treatment_head_reinstatement)

print("Tail Reinstatement vs Original Baseline:")
print("Control Condition:")
print(control_tail_reinstatement)
print("Treatment Condition:")
print(treatment_tail_reinstatement)

# Calculate effect sizes (Cohen's d)
# For Head Reinstatement
control_head_d <- head_reinstatement_paired %>%
  filter(Condition == "Control") %>%
  summarise(
    mean_diff = mean(Head - Original, na.rm = TRUE),
    sd_diff = sd(Head - Original, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

treatment_head_d <- head_reinstatement_paired %>%
  filter(Condition == "Treatment") %>%
  summarise(
    mean_diff = mean(Head - Original, na.rm = TRUE),
    sd_diff = sd(Head - Original, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

# For Tail Reinstatement
control_tail_d <- tail_reinstatement_paired %>%
  filter(Condition == "Control") %>%
  summarise(
    mean_diff = mean(Tail - Original, na.rm = TRUE),
    sd_diff = sd(Tail - Original, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

treatment_tail_d <- tail_reinstatement_paired %>%
  filter(Condition == "Treatment") %>%
  summarise(
    mean_diff = mean(Tail - Original, na.rm = TRUE),
    sd_diff = sd(Tail - Original, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

print("Effect sizes for Head Reinstatement (Cohen's d):")
print(control_head_d)
print(treatment_head_d)

print("Effect sizes for Tail Reinstatement (Cohen's d):")
print(control_tail_d)
print(treatment_tail_d)

# Create visualizations

# 1. Head Reinstatement Paired Visualization
head_reinstatement_viz <- Reinstatement_comparison %>%
  filter(BodyPart %in% c("Original", "Head")) %>%
  ggplot(aes(x = BodyPart, y = ActiveProportion)) +
  # Add individual subject lines
  geom_line(aes(group = OrigSubject, color = Condition), alpha = 0.3) +
  # Add points
  geom_point(aes(color = Condition), size = 3) +
  # Add means and error bars
  stat_summary(
    aes(group = Condition, color = Condition),
    fun = mean, geom = "point", size = 5, shape = 18
  ) +
  stat_summary(
    aes(group = Condition, color = Condition),
    fun.data = function(x) {
      return(c(y = mean(x), ymin = mean(x) - sd(x)/sqrt(length(x)), 
               ymax = mean(x) + sd(x)/sqrt(length(x))))
    },
    geom = "errorbar", width = 0.2, linewidth = 1
  ) +
  # Add p-value annotations
  annotate("text", x = 1.5, y = 0.9, 
           label = paste0("Control: p = ", format(control_head_reinstatement$p.value, digits = 3)), 
           color = "#159090", fontface = "bold", size = 3.5) +
  annotate("text", x = 1.5, y = 0.95, 
           label = paste0("Treatment: p = ", format(treatment_head_reinstatement$p.value, digits = 3)), 
           color = "#FF8C00", fontface = "bold", size = 3.5) +
  # Customize aesthetics
  scale_color_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  labs(
    title = "Head Reinstatement vs Original Baseline",
    subtitle = "Lines connect individual subjects across measurements",
    y = "Proportion of Active Arm Choices",
    x = "Body Status",
    color = "Condition"
  ) +
  ylim(0, 1) +
  consistent_theme() +
  facet_grid(. ~ Condition)

# 2. Tail Reinstatement Paired Visualization
tail_reinstatement_viz <- Reinstatement_comparison %>%
  filter(BodyPart %in% c("Original", "Tail")) %>%
  ggplot(aes(x = BodyPart, y = ActiveProportion)) +
  # Add individual subject lines
  geom_line(aes(group = OrigSubject, color = Condition), alpha = 0.3) +
  # Add points
  geom_point(aes(color = Condition), size = 3) +
  # Add means and error bars
  stat_summary(
    aes(group = Condition, color = Condition),
    fun = mean, geom = "point", size = 5, shape = 18
  ) +
  stat_summary(
    aes(group = Condition, color = Condition),
    fun.data = function(x) {
      return(c(y = mean(x), ymin = mean(x) - sd(x)/sqrt(length(x)), 
               ymax = mean(x) + sd(x)/sqrt(length(x))))
    },
    geom = "errorbar", width = 0.2, linewidth = 1
  ) +
  # Add p-value annotations
  annotate("text", x = 1.5, y = 0.9, 
           label = paste0("Control: p = ", format(control_tail_reinstatement$p.value, digits = 3)), 
           color = "#159090", fontface = "bold", size = 3.5) +
  annotate("text", x = 1.5, y = 0.95, 
           label = paste0("Treatment: p = ", format(treatment_tail_reinstatement$p.value, digits = 3)), 
           color = "#FF8C00", fontface = "bold", size = 3.5) +
  # Customize aesthetics
  scale_color_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  labs(
    title = "Tail Reinstatement vs Original Baseline",
    subtitle = "Lines connect individual subjects across measurements",
    y = "Proportion of Active Arm Choices",
    x = "Body Status",
    color = "Condition"
  ) +
  ylim(0, 1) +
  consistent_theme() +
  facet_grid(. ~ Condition)

# 3. Combined bar chart for publication
combined_means <- bind_rows(
  # Head data
  head_reinstatement_paired %>%
    group_by(Condition) %>%
    summarise(
      BodyPart = "Head",
      Original_Mean = mean(Original, na.rm = TRUE),
      Original_SE = sd(Original, na.rm = TRUE)/sqrt(n()),
      Regen_Mean = mean(Head, na.rm = TRUE),
      Regen_SE = sd(Head, na.rm = TRUE)/sqrt(n()),
      .groups = 'drop'
    ),
  # Tail data
  tail_reinstatement_paired %>%
    group_by(Condition) %>%
    summarise(
      BodyPart = "Tail",
      Original_Mean = mean(Original, na.rm = TRUE),
      Original_SE = sd(Original, na.rm = TRUE)/sqrt(n()),
      Regen_Mean = mean(Tail, na.rm = TRUE),
      Regen_SE = sd(Tail, na.rm = TRUE)/sqrt(n()),
      .groups = 'drop'
    )
) %>%
  # Create a measurement variable to distinguish between baseline and reinstatement
  mutate(
    Baseline_label = paste(BodyPart, "Baseline"),
    Reinstatement_label = paste(BodyPart, "Reinstatement")
  ) %>%
  # Convert to long format for plotting
  pivot_longer(
    cols = c(Original_Mean, Regen_Mean),
    names_to = "Measurement_type",
    values_to = "Mean"
  ) %>%
  mutate(
    Measurement = ifelse(Measurement_type == "Original_Mean", Baseline_label, Reinstatement_label),
    SE = ifelse(Measurement_type == "Original_Mean", Original_SE, Regen_SE)
  )

# Create the bar chart
reinstatement_bar_chart <- combined_means %>%
  ggplot(aes(x = Measurement, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SE, ymax = Mean + SE),
    position = position_dodge(0.9), width = 0.2
  ) +
  # Add significance stars if applicable (based on p-values)
  scale_fill_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  labs(
    title = "Reinstatement Performance vs Original Baseline",
    y = "Proportion of Active Arm Choices",
    x = "",
    fill = "Condition"
  ) +
  ylim(0, 1) +
  consistent_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Print the visualizations
print(head_reinstatement_viz)
print(tail_reinstatement_viz)
print(reinstatement_bar_chart)


#=================================================================
# PART 9: VISUALIZATION OF BASELINE VS. REGENERATtion and Reinstatement
#=================================================================

###### Initial Regeneration Plot

# Create a long format dataset for individual subject visualization
Individual_comparison <- bind_rows(
  # Original baseline data
  Baseline_scores %>%
    select(Subject, Condition, BaselineProportion) %>%
    rename(
      OrigSubject = Subject,
      ActiveProportion = BaselineProportion
    ) %>%
    mutate(
      BodyPart = "Original",
      TestPhase = "Baseline"
    ),
  
  # Memory test data for regenerated parts
  Regen_vs_baseline %>%
    select(OrigSubject, BodyPart, Condition, ActiveProportion) %>%
    mutate(TestPhase = "Memory")
) %>%
  mutate(
    OrigSubject = as.character(OrigSubject), # Ensure consistent type
    BodyPart = factor(BodyPart, levels = c("Original", "Head", "Tail")),
    TestPhase = factor(TestPhase, levels = c("Baseline", "Memory")),
    Condition = factor(Condition)
  )
# Create a cleaner version of the plot with only means and SE bars
# First, calculate group means and standard errors
group_means <- Individual_comparison %>%
  group_by(Condition, BodyPart, TestPhase) %>%
  summarise(
    Mean = mean(ActiveProportion, na.rm = TRUE),
    SE = sd(ActiveProportion, na.rm = TRUE)/sqrt(n()),
    .groups = 'drop'
  )

# Create the final plot with no connecting lines
regeneration_plot <- ggplot(
  group_means,
  aes(x = BodyPart, y = Mean)
) +
  # Add mean points with different shapes by TestPhase and BodyPart
  geom_point(
    aes(shape = interaction(TestPhase, BodyPart), fill = Condition),
    size = 5, stroke = 0.7 # Thin stroke as requested
  ) +
  # Add error bars
  geom_errorbar(
    aes(ymin = Mean - SE, ymax = Mean + SE, color = Condition),
    width = 0.2, linewidth = 0.7 # Thin error bars
  ) +
  # NO connecting lines in this version
  
  # Customize aesthetics
  scale_fill_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  scale_color_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  # Custom shape scale for TestPhase and BodyPart combinations
  scale_shape_manual(
    name = "Regeneration Test",
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
  # Use your consistent theme
  consistent_theme() +
  theme(
    legend.position = "right",
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  facet_grid(. ~ Condition) +
  guides(
    shape = guide_legend(title = "Legend", override.aes = list(stroke = 0.7)),
    fill = guide_legend(title = "Condition", override.aes = list(shape = 21, stroke = 0.7)),
    color = "none"  # Hide color legend since it's redundant
  )

# Print the final plot
print(regeneration_plot)


###### Reinstatement plot

reinstatement_group_means <- Reinstatement_comparison %>%
  group_by(Condition, BodyPart, TestPhase) %>%
  summarise(
    Mean = mean(ActiveProportion, na.rm = TRUE),
    SD = sd(ActiveProportion, na.rm = TRUE),
    N = n(),
    SE = SD/sqrt(N),
    .groups = 'drop'
  )

# Print the means to verify
print(reinstatement_group_means)


reinstatement_plot <- ggplot(
  reinstatement_group_means,
  aes(x = BodyPart, y = Mean)
) +
  # Add mean points with different shapes by TestPhase and BodyPart
  geom_point(
    aes(shape = interaction(TestPhase, BodyPart), fill = Condition),
    size = 5, stroke = 0.7 # Thin stroke as requested
  ) +
  # Add error bars
  geom_errorbar(
    aes(ymin = Mean - SE, ymax = Mean + SE, color = Condition),
    width = 0.2, linewidth = 0.7 # Thin error bars
  ) +
  
  # Customize aesthetics
  scale_fill_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  scale_color_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  # Custom shape scale for TestPhase and BodyPart combinations
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
      "Baseline.Original" = "Intact Bseline",
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
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  # Use your consistent theme
  consistent_theme() +
  theme(
    legend.position = "right",
  ) +
  facet_grid(. ~ Condition) +
  guides(
    shape = guide_legend(title = "Legend", override.aes = list(stroke = 0.7)),
    fill = guide_legend(title = "Condition", override.aes = list(shape = 21, stroke = 0.7)),
    color = "none"  # Hide color legend since it's redundant
  )

# Print the reinstatement plot
print(reinstatement_plot)

# Save the reinstatement plot
ggsave(
  "Reinstatement_memory_retention.pdf", 
  reinstatement_plot, 
  width = 12,
  height = 8,
  dpi = 300,
  device = cairo_pdf
)



###### the retention and reinstatement plots for regenerates

# Combine the regeneration and reinstatement plots

# Remvoe legend from one of the grapghs
reinstatement_plot_no_legend <- reinstatement_plot +
  theme(legend.position = "none")

# Combine the plots, keeping only the regeneration_plot legend
Regen_combined_figure <- (regeneration_plot + reinstatement_plot_no_legend) +
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

regeneration_plot_no_legend <- regeneration_plot + 
  theme(legend.position = "none")

##### Ordering the legend

reinstatement_plot_ordered <- reinstatement_plot +
  guides(
    # Keep the fill legend with its colors for the Condition
    fill = guide_legend(
      order = 1, 
      title = "Condition",
      override.aes = list(shape = 21, size = 5, stroke = 0.7)
    ),
    # Set up the shape legend
    shape = guide_legend(
      order = 2, 
      title = "Legend",
      override.aes = list(stroke = 0.7, size = 5)
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
      plot.tag = element_text(face = "bold", size = 16, family = "Times New Roman"),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
  )


print(combined_figure)


# Save the panel
ggsave("Exp7_combined_figure.png", combined_figure, 
       width = 16, height = 14, dpi = 300) 
