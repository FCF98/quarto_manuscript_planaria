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
library(ggsignif)
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

########################################
#formatting results in apa7 
######################################


# Format p-value according to APA style
format_p_value <- function(p) {
  if (p < .001) {
    return("< .001")
  } else {
    return(paste0("= ", gsub("0\\.", ".", round(p, 3))))
  }
}

# Create formatted results for the intact GLMM analysis
format_intact_glmm_results <- function(model, anova_output) {
  # Extract model summary
  model_summary <- summary(model)
  
  # Create a dataframe for the main effects and interaction
  anova_df <- as.data.frame(anova_output)
  
  # Format the results for Condition, Time and their interaction
  formatted_results <- data.frame(
    effect = rownames(anova_df),
    chisq = anova_df$Chisq,
    df = anova_df$Df,
    p.value = anova_df$`Pr(>Chisq)`,
    apa_result = paste0(
      "χ²(", anova_df$Df, ") = ", round(anova_df$Chisq, 3),
      ", *p* ", sapply(anova_df$`Pr(>Chisq)`, format_p_value)
    )
  )
  
  return(formatted_results)
}

# Format emmeans pairwise comparisons for within-group effects
format_within_group_comparisons <- function(emmeans_object) {
  # Extract contrasts
  contrasts_df <- as.data.frame(emmeans_object$contrasts)
  
  condition_names <- c("1" = "Control", "2" = "Treatment")
  
  # Create formatted results dataframe
  formatted_results <- data.frame(
    contrast = contrasts_df$contrast,
    condition = sub(":.*", "", rownames(contrasts_df)),
    odds_ratio = contrasts_df$odds.ratio,
    p.value = contrasts_df$p.value,
    apa_result = paste0(
      "OR = ", round(contrasts_df$odds.ratio, 2),
      ", *z* = ", round(contrasts_df$z.ratio, 2),
      ", *p* ", sapply(contrasts_df$p.value, format_p_value)
    )
  )
  
  # Apply the condition mapping to transform numeric codes to names
  formatted_results$condition <- condition_names[formatted_results$condition]
  
  return(formatted_results)
}

# Format emmeans pairwise comparisons for between-group effects
format_between_group_comparisons <- function(emmeans_object) {
  # Extract contrasts
  contrasts_df <- as.data.frame(emmeans_object$contrasts)
  
  # Create formatted results dataframe
  formatted_results <- data.frame(
    contrast = contrasts_df$contrast,
    time_point = sub(":.*", "", rownames(contrasts_df)),
    odds_ratio = contrasts_df$odds.ratio,
    p.value = contrasts_df$p.value,
    apa_result = paste0(
      "OR = ", round(contrasts_df$odds.ratio, 2),
      ", *z* = ", round(contrasts_df$z.ratio, 2),
      ", *p* ", sapply(contrasts_df$p.value, format_p_value)
    )
  )
  
  return(formatted_results)
}

# Format coefficients from the model
format_model_coefficients <- function(model) {
  # Extract coefficient table
  coef_table <- as.data.frame(summary(model)$coefficients)
  
  # Create formatted results dataframe
  formatted_results <- data.frame(
    term = rownames(coef_table),
    estimate = coef_table$Estimate,
    se = coef_table$`Std. Error`,
    z_value = coef_table$`z value`,
    p.value = coef_table$`Pr(>|z|)`,
    apa_result = paste0(
      "β = ", round(coef_table$Estimate, 2),
      ", SE = ", round(coef_table$`Std. Error`, 2),
      ", *z* = ", round(coef_table$`z value`, 2),
      ", *p* ", sapply(coef_table$`Pr(>|z|)`, format_p_value)
    )
  )
  
  return(formatted_results)
}

# Format estimated means for reporting
format_estimated_means <- function(emmeans_object) {
  # Extract emmeans
  means_df <- as.data.frame(emmeans_object$emmeans)
  
  # Create formatted results dataframe
  formatted_results <- data.frame(
    group = rownames(means_df),
    condition = means_df$Condition,
    time = means_df$Time,
    probability = means_df$prob,
    se = means_df$SE,
    lcl = means_df$asymp.LCL,
    ucl = means_df$asymp.UCL,
    apa_result = paste0(
      "*M* = ", round(means_df$prob, 2),
      ", *SE* = ", round(means_df$SE, 2),
      ", 95% CI [", round(means_df$asymp.LCL, 2), ", ", round(means_df$asymp.UCL, 2), "]"
    )
  )
  
  return(formatted_results)
}

# Run all formatting functions to create results tables
# This is the main function you'll call
format_all_glmm_results <- function(model, anova_output, within_comparisons, between_comparisons) {
  # Create a list to store all formatted results
  results <- list()
  
  # Format all the different result components
  results$anova <- format_intact_glmm_results(model, anova_output)
  results$coefficients <- format_model_coefficients(model)
  results$within_comparisons <- format_within_group_comparisons(within_comparisons)
  results$between_comparisons <- format_between_group_comparisons(between_comparisons)
  results$estimated_means <- format_estimated_means(within_comparisons)
  
  return(results)
}

# Helper function to make retrieving results easier for inline code
get_result <- function(results, type, filter_value = NULL, filter_col = NULL, 
                       condition = NULL, time = NULL) {
  
  # Handle main effects and interactions (anova)
  if (type == "anova") {
    # Get specific effect (Condition, Time, or Condition:Time)
    if (!is.null(filter_value)) {
      return(results$anova$apa_result[results$anova$effect == filter_value])
    } else {
      # Return all effects as a named vector
      effects <- results$anova$apa_result
      names(effects) <- results$anova$effect
      return(effects)
    }
  }
  
  # Handle within-group comparisons
  if (type == "within") {
    # If condition specified, filter by it
    if (!is.null(condition)) {
      # Find the row that matches the condition
      matching_row <- which(results$within_comparisons$condition == condition)
      
      # Check if we found a match
      if (length(matching_row) > 0) {
        return(results$within_comparisons$apa_result[matching_row])
      } else {
        # For debugging - you can remove this in the final version
        return(paste0("(No match for '", condition, "')"))
      }
    } else {
      # Return all as named vector
      comps <- results$within_comparisons$apa_result
      names(comps) <- results$within_comparisons$condition
      return(comps)
    }
  }
  
  # Handle between-group comparisons
  if (type == "between") {
    # If time specified, filter by it
    if (!is.null(time)) {
      time_filter <- paste0("Time = ", time)
      return(results$between_comparisons$apa_result[
        results$between_comparisons$time_point == time_filter])
    } else {
      # Return all as named vector
      comps <- results$between_comparisons$apa_result
      names(comps) <- results$between_comparisons$time_point
      return(comps)
    }
  }
  
  # Handle estimated means
  if (type == "means") {
    # Filter by condition and/or time if specified
    if (!is.null(condition) && !is.null(time)) {
      return(results$estimated_means$apa_result[
        results$estimated_means$condition == condition & 
          results$estimated_means$time == time])
    } else if (!is.null(condition)) {
      return(results$estimated_means$apa_result[
        results$estimated_means$condition == condition])
    } else if (!is.null(time)) {
      return(results$estimated_means$apa_result[
        results$estimated_means$time == time])
    } else {
      # Return all as named vector
      means <- results$estimated_means$apa_result
      names(means) <- paste0(results$estimated_means$condition, " - ", 
                             results$estimated_means$time)
      return(means)
    }
  }
  
  # Handle model coefficients
  if (type == "coef") {
    if (!is.null(filter_value)) {
      return(results$coefficients$apa_result[
        results$coefficients$term == filter_value])
    } else {
      # Return all as named vector
      coefs <- results$coefficients$apa_result
      names(coefs) <- results$coefficients$term
      return(coefs)
    }
  }
  
  # Return NULL if no match
  return(NULL)
}

# Apply the formatting to your existing analysis
# Just run this after your model and emmeans calculations
glmm_results <- format_all_glmm_results(
  Intact_model,
  intact_model_output,
  intact_within_group_comparisons,
  intact_between_group_comparisons
)


######### creating plot of baseline vs endpoint ###########


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
    title = "Active Arm Entries Before and After Conditioning",
    y = "Proportion of Active Arm Entries",
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
                                  "Conditioning_day5"),
                       labels = c("BL1", "BL2", "CD1", "CD2", "CD3", "CD4", "CD5")),
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
    title = "Active Arm Entries Throughout Conditionings",
    y = "Proportion of Active Arm Entries",
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
    title = "Reinstatement Compared to Baseline",
    y = "Proportion of Active Arm Entries",
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


# Create a paired data structure for tails similar to what you did for heads
tail_paired_data <- Individual_comparison %>%
  filter(BodyPart %in% c("Original", "Tail")) %>%
  select(OrigSubject, Condition, BodyPart, ActiveProportion) %>%
  pivot_wider(
    names_from = BodyPart,
    values_from = ActiveProportion
  )

# Check the structure of the paired data
print(head(tail_paired_data, 10))

# Perform paired t-tests for each condition (Tail vs Original)
# Control condition
control_tail_test <- tail_paired_data %>%
  filter(Condition == "Control") %>%
  with(t.test(Tail, Original, paired = TRUE))

# Treatment condition
treatment_tail_test <- tail_paired_data %>%
  filter(Condition == "Treatment") %>%
  with(t.test(Tail, Original, paired = TRUE))

# Print results
print("Control Condition: Tail vs Original Baseline (Paired t-test)")
print(control_tail_test)

print("Treatment Condition: Tail vs Original Baseline (Paired t-test)")
print(treatment_tail_test)

# Calculate effect sizes (Cohen's d for paired data)
# Control effect size
control_tail_d <- tail_paired_data %>%
  filter(Condition == "Control") %>%
  summarise(
    mean_diff = mean(Tail - Original, na.rm = TRUE),
    sd_diff = sd(Tail - Original, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

# Treatment effect size
treatment_tail_d <- tail_paired_data %>%
  filter(Condition == "Treatment") %>%
  summarise(
    mean_diff = mean(Tail - Original, na.rm = TRUE),
    sd_diff = sd(Tail - Original, na.rm = TRUE),
    cohen_d = mean_diff / sd_diff,
    n = n()
  )

print("Effect sizes for Tail Memory Test (Cohen's d):")
print(control_tail_d)
print(treatment_tail_d)



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


###################################################################
# formatting results in apa7
#################################################################

# Format p-value according to APA style
format_p_value <- function(p) {
  if (p < .001) {
    return("< .001")
  } else {
    return(paste0("= ", gsub("0\\.", ".", round(p, 3))))
  }
}

# Function to format a t-test result in APA style
format_t_test <- function(t_test_result, d_value) {
  # Extract values from t.test object
  t_value <- t_test_result$statistic
  df <- t_test_result$parameter
  p_value <- t_test_result$p.value
  
  # Format the APA style result
  formatted_result <- paste0(
    "*t*(", round(df, 2), ") = ", round(t_value, 2),
    ", *d* = ", round(abs(d_value), 2),
    ", *p* ", format_p_value(p_value)
  )
  
  return(formatted_result)
}

# Function to format means and SDs in APA style
format_descriptives <- function(mean_value, se_value, n_value = NULL) {
  if (is.null(n_value)) {
    return(paste0("*M* = ", round(mean_value, 2), ", *SE* = ", round(se_value, 2)))
  } else {
    return(paste0("*M* = ", round(mean_value, 2), ", *SE* = ", round(se_value, 2), ", *n* = ", n_value))
  }
}

# Function to format memory retention results for both head and tail
format_memory_retention_results <- function(
    control_head_test, treatment_head_test,
    control_tail_test, treatment_tail_test,
    control_head_d, treatment_head_d,
    control_tail_d, treatment_tail_d,
    head_paired_data, tail_paired_data) {
  
  # Calculate descriptive statistics for head
  head_descriptives <- head_paired_data %>%
    group_by(Condition) %>%
    summarise(
      n = n(),
      Original_Mean = mean(Original, na.rm = TRUE),
      Original_SD = sd(Original, na.rm = TRUE),
      Original_SE = Original_SD / sqrt(n),
      Head_Mean = mean(Head, na.rm = TRUE),
      Head_SD = sd(Head, na.rm = TRUE),
      Head_SE = Head_SD / sqrt(n),
      .groups = 'drop'
    )
  
  # Calculate descriptive statistics for tail if tail data exists
  if(!is.null(tail_paired_data) && nrow(tail_paired_data) > 0) {
    tail_descriptives <- tail_paired_data %>%
      group_by(Condition) %>%
      summarise(
        n = n(),
        Original_Mean = mean(Original, na.rm = TRUE),
        Original_SD = sd(Original, na.rm = TRUE),
        Original_SE = Original_SD / sqrt(n),
        Tail_Mean = mean(Tail, na.rm = TRUE),
        Tail_SD = sd(Tail, na.rm = TRUE),
        Tail_SE = Tail_SD / sqrt(n),
        .groups = 'drop'
      )
  } else {
    tail_descriptives <- NULL
  }
  
  # Format t-test results for head
  control_head_result <- format_t_test(control_head_test, control_head_d$cohen_d)
  treatment_head_result <- format_t_test(treatment_head_test, treatment_head_d$cohen_d)
  
  # Format t-test results for tail if tail tests exist
  if(!is.null(control_tail_test) && !is.null(treatment_tail_test) && 
     !is.null(control_tail_d) && !is.null(treatment_tail_d)) {
    control_tail_result <- format_t_test(control_tail_test, control_tail_d$cohen_d)
    treatment_tail_result <- format_t_test(treatment_tail_test, treatment_tail_d$cohen_d)
  } else {
    control_tail_result <- NULL
    treatment_tail_result <- NULL
  }
  
  # Format descriptive statistics for head
  control_original_head_desc <- format_descriptives(
    head_descriptives$Original_Mean[head_descriptives$Condition == "Control"],
    head_descriptives$Original_SE[head_descriptives$Condition == "Control"],
    head_descriptives$n[head_descriptives$Condition == "Control"]
  )
  
  control_head_desc <- format_descriptives(
    head_descriptives$Head_Mean[head_descriptives$Condition == "Control"],
    head_descriptives$Head_SE[head_descriptives$Condition == "Control"]
  )
  
  treatment_original_head_desc <- format_descriptives(
    head_descriptives$Original_Mean[head_descriptives$Condition == "Treatment"],
    head_descriptives$Original_SE[head_descriptives$Condition == "Treatment"],
    head_descriptives$n[head_descriptives$Condition == "Treatment"]
  )
  
  treatment_head_desc <- format_descriptives(
    head_descriptives$Head_Mean[head_descriptives$Condition == "Treatment"],
    head_descriptives$Head_SE[head_descriptives$Condition == "Treatment"]
  )
  
  # Format descriptive statistics for tail if tail data exists
  if(!is.null(tail_descriptives)) {
    control_original_tail_desc <- format_descriptives(
      tail_descriptives$Original_Mean[tail_descriptives$Condition == "Control"],
      tail_descriptives$Original_SE[tail_descriptives$Condition == "Control"],
      tail_descriptives$n[tail_descriptives$Condition == "Control"]
    )
    
    control_tail_desc <- format_descriptives(
      tail_descriptives$Tail_Mean[tail_descriptives$Condition == "Control"],
      tail_descriptives$Tail_SE[tail_descriptives$Condition == "Control"]
    )
    
    treatment_original_tail_desc <- format_descriptives(
      tail_descriptives$Original_Mean[tail_descriptives$Condition == "Treatment"],
      tail_descriptives$Original_SE[tail_descriptives$Condition == "Treatment"],
      tail_descriptives$n[tail_descriptives$Condition == "Treatment"]
    )
    
    treatment_tail_desc <- format_descriptives(
      tail_descriptives$Tail_Mean[tail_descriptives$Condition == "Treatment"],
      tail_descriptives$Tail_SE[tail_descriptives$Condition == "Treatment"]
    )
  } else {
    control_original_tail_desc <- NULL
    control_tail_desc <- NULL
    treatment_original_tail_desc <- NULL
    treatment_tail_desc <- NULL
  }
  
  # Create a list of all formatted results
  results <- list(
    # Head retention results
    control_head_t_test = control_head_result,
    treatment_head_t_test = treatment_head_result,
    control_original_head_desc = control_original_head_desc,
    control_head_desc = control_head_desc,
    treatment_original_head_desc = treatment_original_head_desc,
    treatment_head_desc = treatment_head_desc
  )
  
  # Add tail results if they exist
  if(!is.null(control_tail_result)) {
    results$control_tail_t_test <- control_tail_result
    results$treatment_tail_t_test <- treatment_tail_result
    results$control_original_tail_desc <- control_original_tail_desc
    results$control_tail_desc <- control_tail_desc
    results$treatment_original_tail_desc <- treatment_original_tail_desc
    results$treatment_tail_desc <- treatment_tail_desc
  }
  
  return(results)
}

# Function to format the planaria reinstatement results
format_reinstatement_results <- function(control_head_test, treatment_head_test, 
                                         control_tail_test, treatment_tail_test,
                                         control_head_d, treatment_head_d,
                                         control_tail_d, treatment_tail_d,
                                         head_paired_data, tail_paired_data) {
  
  # Calculate descriptive statistics for head
  head_descriptives <- head_paired_data %>%
    group_by(Condition) %>%
    summarise(
      n = n(),
      Original_Mean = mean(Original, na.rm = TRUE),
      Original_SD = sd(Original, na.rm = TRUE),
      Original_SE = Original_SD / sqrt(n),
      Head_Mean = mean(Head, na.rm = TRUE),
      Head_SD = sd(Head, na.rm = TRUE),
      Head_SE = Head_SD / sqrt(n),
      .groups = 'drop'
    )
  
  # Calculate descriptive statistics for tail
  tail_descriptives <- tail_paired_data %>%
    group_by(Condition) %>%
    summarise(
      n = n(),
      Original_Mean = mean(Original, na.rm = TRUE),
      Original_SD = sd(Original, na.rm = TRUE),
      Original_SE = Original_SD / sqrt(n),
      Tail_Mean = mean(Tail, na.rm = TRUE),
      Tail_SD = sd(Tail, na.rm = TRUE),
      Tail_SE = Tail_SD / sqrt(n),
      .groups = 'drop'
    )
  
  # Format t-test results for head
  control_head_result <- format_t_test(control_head_test, control_head_d$cohen_d)
  treatment_head_result <- format_t_test(treatment_head_test, treatment_head_d$cohen_d)
  
  # Format t-test results for tail
  control_tail_result <- format_t_test(control_tail_test, control_tail_d$cohen_d)
  treatment_tail_result <- format_t_test(treatment_tail_test, treatment_tail_d$cohen_d)
  
  # Format descriptive statistics for head
  control_original_head_desc <- format_descriptives(
    head_descriptives$Original_Mean[head_descriptives$Condition == "Control"],
    head_descriptives$Original_SE[head_descriptives$Condition == "Control"],
    head_descriptives$n[head_descriptives$Condition == "Control"]
  )
  
  control_head_desc <- format_descriptives(
    head_descriptives$Head_Mean[head_descriptives$Condition == "Control"],
    head_descriptives$Head_SE[head_descriptives$Condition == "Control"]
  )
  
  treatment_original_head_desc <- format_descriptives(
    head_descriptives$Original_Mean[head_descriptives$Condition == "Treatment"],
    head_descriptives$Original_SE[head_descriptives$Condition == "Treatment"],
    head_descriptives$n[head_descriptives$Condition == "Treatment"]
  )
  
  treatment_head_desc <- format_descriptives(
    head_descriptives$Head_Mean[head_descriptives$Condition == "Treatment"],
    head_descriptives$Head_SE[head_descriptives$Condition == "Treatment"]
  )
  
  # Format descriptive statistics for tail
  control_original_tail_desc <- format_descriptives(
    tail_descriptives$Original_Mean[tail_descriptives$Condition == "Control"],
    tail_descriptives$Original_SE[tail_descriptives$Condition == "Control"],
    tail_descriptives$n[tail_descriptives$Condition == "Control"]
  )
  
  control_tail_desc <- format_descriptives(
    tail_descriptives$Tail_Mean[tail_descriptives$Condition == "Control"],
    tail_descriptives$Tail_SE[tail_descriptives$Condition == "Control"]
  )
  
  treatment_original_tail_desc <- format_descriptives(
    tail_descriptives$Original_Mean[tail_descriptives$Condition == "Treatment"],
    tail_descriptives$Original_SE[tail_descriptives$Condition == "Treatment"],
    tail_descriptives$n[tail_descriptives$Condition == "Treatment"]
  )
  
  treatment_tail_desc <- format_descriptives(
    tail_descriptives$Tail_Mean[tail_descriptives$Condition == "Treatment"],
    tail_descriptives$Tail_SE[tail_descriptives$Condition == "Treatment"]
  )
  
  # Create a list of all formatted results
  results <- list(
    # Head reinstatement results
    control_head_t_test = control_head_result,
    treatment_head_t_test = treatment_head_result,
    control_original_head_desc = control_original_head_desc,
    control_head_desc = control_head_desc,
    treatment_original_head_desc = treatment_original_head_desc,
    treatment_head_desc = treatment_head_desc,
    
    # Tail reinstatement results
    control_tail_t_test = control_tail_result,
    treatment_tail_t_test = treatment_tail_result,
    control_original_tail_desc = control_original_tail_desc,
    control_tail_desc = control_tail_desc,
    treatment_original_tail_desc = treatment_original_tail_desc,
    treatment_tail_desc = treatment_tail_desc
  )
  
  return(results)
}

# Helper function to make getting results easier
get_planaria_result <- function(results, part = NULL, condition = NULL, measurement = NULL, bodypart = NULL) {
  # For memory retention results
  if (is.null(part) || part == "memory") {
    if (!is.null(condition) && !is.null(measurement)) {
      # For head measurements
      if (is.null(bodypart) || bodypart == "head") {
        # Return specific descriptive
        if (condition == "Control" && measurement == "Original") {
          return(results$memory$control_original_head_desc)
        } else if (condition == "Control" && measurement == "Head") {
          return(results$memory$control_head_desc)
        } else if (condition == "Treatment" && measurement == "Original") {
          return(results$memory$treatment_original_head_desc)
        } else if (condition == "Treatment" && measurement == "Head") {
          return(results$memory$treatment_head_desc)
        }
      }
      # For tail measurements, if they exist
      else if (bodypart == "tail" && !is.null(results$memory$control_tail_t_test)) {
        if (condition == "Control" && measurement == "Original") {
          return(results$memory$control_original_tail_desc)
        } else if (condition == "Control" && measurement == "Tail") {
          return(results$memory$control_tail_desc)
        } else if (condition == "Treatment" && measurement == "Original") {
          return(results$memory$treatment_original_tail_desc)
        } else if (condition == "Treatment" && measurement == "Tail") {
          return(results$memory$treatment_tail_desc)
        }
      }
    } else if (!is.null(condition) && is.null(measurement)) {
      # Return t-test result
      if (is.null(bodypart) || bodypart == "head") {
        if (condition == "Control") {
          return(results$memory$control_head_t_test)
        } else if (condition == "Treatment") {
          return(results$memory$treatment_head_t_test)
        }
      } else if (bodypart == "tail" && !is.null(results$memory$control_tail_t_test)) {
        if (condition == "Control") {
          return(results$memory$control_tail_t_test)
        } else if (condition == "Treatment") {
          return(results$memory$treatment_tail_t_test)
        }
      }
    }
  }
  
  # For reinstatement results
  else if (part == "reinstatement") {
    if (!is.null(condition) && !is.null(measurement)) {
      # Return specific descriptive for head
      if (is.null(bodypart) || bodypart == "head") {
        if (condition == "Control" && measurement == "Original") {
          return(results$reinstatement$control_original_head_desc)
        } else if (condition == "Control" && measurement == "Head") {
          return(results$reinstatement$control_head_desc)
        } else if (condition == "Treatment" && measurement == "Original") {
          return(results$reinstatement$treatment_original_head_desc)
        } else if (condition == "Treatment" && measurement == "Head") {
          return(results$reinstatement$treatment_head_desc)
        }
      }
      # Return specific descriptive for tail
      else if (bodypart == "tail") {
        if (condition == "Control" && measurement == "Original") {
          return(results$reinstatement$control_original_tail_desc)
        } else if (condition == "Control" && measurement == "Tail") {
          return(results$reinstatement$control_tail_desc)
        } else if (condition == "Treatment" && measurement == "Original") {
          return(results$reinstatement$treatment_original_tail_desc)
        } else if (condition == "Treatment" && measurement == "Tail") {
          return(results$reinstatement$treatment_tail_desc)
        }
      }
    } else if (!is.null(condition) && is.null(measurement)) {
      # Return t-test result
      if (is.null(bodypart) || bodypart == "head") {
        if (condition == "Control") {
          return(results$reinstatement$control_head_t_test)
        } else if (condition == "Treatment") {
          return(results$reinstatement$treatment_head_t_test)
        }
      } else if (bodypart == "tail") {
        if (condition == "Control") {
          return(results$reinstatement$control_tail_t_test)
        } else if (condition == "Treatment") {
          return(results$reinstatement$treatment_tail_t_test)
        }
      }
    }
  }
  
  # Return NULL if no match
  return(NULL)
}

# Main function to run all the formatting
format_planaria_analyses <- function(
    # Memory retention t-tests and data
  control_head_memory_test, treatment_head_memory_test, 
  control_head_memory_d, treatment_head_memory_d, 
  head_memory_paired_data,
  
  # Optional: tail memory retention tests (may be NULL if not available)
  control_tail_memory_test = NULL, treatment_tail_memory_test = NULL,
  control_tail_memory_d = NULL, treatment_tail_memory_d = NULL,
  tail_memory_paired_data = NULL,
  
  # Reinstatement t-tests and data
  control_head_reinstatement, treatment_head_reinstatement,
  control_tail_reinstatement, treatment_tail_reinstatement,
  control_head_reinstatement_d, treatment_head_reinstatement_d,
  control_tail_reinstatement_d, treatment_tail_reinstatement_d,
  head_reinstatement_paired, tail_reinstatement_paired
) {
  
  # Format memory retention results (for both head and tail if available)
  memory_results <- format_memory_retention_results(
    control_head_memory_test, treatment_head_memory_test,
    control_tail_memory_test, treatment_tail_memory_test,
    control_head_memory_d, treatment_head_memory_d,
    control_tail_memory_d, treatment_tail_memory_d,
    head_memory_paired_data, tail_memory_paired_data
  )
  
  # Format reinstatement results
  reinstatement_results <- format_memory_retention_results(
    control_head_reinstatement, treatment_head_reinstatement,
    control_tail_reinstatement, treatment_tail_reinstatement,
    control_head_reinstatement_d, treatment_head_reinstatement_d,
    control_tail_reinstatement_d, treatment_tail_reinstatement_d,
    head_reinstatement_paired, tail_reinstatement_paired
  )
  
  # Combine all results
  results <- list(
    memory = memory_results,
    reinstatement = reinstatement_results
  )
  
  return(results)
}


Exp7_regen_reinstate_results <- format_planaria_analyses(
  # Memory tests (head)
  control_paired_test, treatment_paired_test,
  control_d, treatment_d, paired_data,
  
  # Memory tests (tail) - now we have these variables
  control_tail_test, treatment_tail_test,
  control_tail_d, treatment_tail_d, tail_paired_data,
  
  # Reinstatement tests (both head and tail)
  control_head_reinstatement, treatment_head_reinstatement,
  control_tail_reinstatement, treatment_tail_reinstatement,
  control_head_d, treatment_head_d,
  control_tail_d, treatment_tail_d,
  head_reinstatement_paired, tail_reinstatement_paired
)




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
  
  # Add significance bars
  # Control head vs original (p < 0.01)
  geom_signif(
    data = subset(group_means, Condition == "Control"),
    annotations = "**",
    xmin = "Original", xmax = "Head",
    y_position = 0.8,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  # Control tail vs original (p < 0.05)
  geom_signif(
    data = subset(group_means, Condition == "Control"),
    annotations = "*",
    xmin = "Original", xmax = "Tail",
    y_position = 0.7,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  # Treatment head vs original (p < 0.001)
  geom_signif(
    data = subset(group_means, Condition == "Treatment"),
    annotations = "***",
    xmin = "Original", xmax = "Head",
    y_position = 0.8,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  
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
    title = "Regeneration Compared to Baseline",
    y = "Proportion of Active Arm Entries",
    x = "Body Status"
  ) +
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

# Display the plot
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
    title = "Reinstatement Compared to Baseline",
    y = "Proportion of Active Arm Entries",
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

# Apply theme to each individual plot
learning_plot_no_legend <- learning_plot_no_legend +
  theme(
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 18, face = "bold", family = "Times New Roman") 
  )

intact_plot_no_legend <- intact_plot_no_legend +
  theme(
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 18, face = "bold", family = "Times New Roman", hjust = 0.5) 
  )

regeneration_plot_no_legend <- regeneration_plot_no_legend +
  theme(
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 18, face = "bold", family = "Times New Roman") 
  )



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

reinstatement_plot_ordered <- reinstatement_plot_ordered +
  theme(
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 18, face = "bold", family = "Times New Roman") 
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


#=======================================================================
# PART 11: descriptives
#=======================================================================

#Pulling baseline arm preference descriptives

Exp7_prefered_arm_count <- Exp7_data %>% filter(Subject %in% (1:30)) %>% count(Prefered_arm)

Exp7_active_arm_count <- Exp7_data %>% filter(Subject %in% (1:30)) %>% count(`Active arm`)

Exp7_left_active_arm_count <- Exp7_active_arm_count %>% filter(`Active arm` == "L") %>% pull(n)

Exp7_right_active_arm_count <- Exp7_active_arm_count %>% filter(`Active arm` == "R") %>% pull(n)

#Pulling active arm descriptives

Left_active_arm_count <- Exp7_prefered_arm_count

