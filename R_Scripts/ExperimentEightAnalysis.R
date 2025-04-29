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
library(patchwork)

# read in the data
Exp8_full_data <- read_excel("Datasets/Experiment_8_Full_Data.xlsx")


#Creating theme for plots

# Updated theme for regeneration and reinstatement plots
consistant_theme <- function() {
  theme_classic() +
    theme(
      text = element_text(family = "Times New Roman", size = 20),
      axis.title = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(margin = margin(r = 20)),
      axis.text = element_text(size = 18, color = "black"),
      legend.position = "right", # Updated to place legend on right side
      legend.title = element_text(size = 17, face = "bold"),
      legend.text = element_text(size = 16),
      axis.line = element_line(color = "black", linewidth = 0.8),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      # Add facet theming to match the figure
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 14, face = "bold", color = "black")
    )
}

#Restructure data for baseline/endpoint comparison

Exp8_data_long <- Exp8_full_data %>% mutate(across(ends_with("_arm_entries"),
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
                                  ActiveCount/ TotalTrials,
                                  NA_real_)
  ) %>% 
  filter(is.na(ActiveCount) | ActiveCount <= TotalTrials)
  

# statistical analysis of baseline vs endpoint

#set contrasts
contrasts(Exp8_data_long$Time) <- contr.sum
contrasts(Exp8_data_long$Condition) <- contr.sum

#baseline vs endpoint model
learning_model <- glmer(cbind(ActiveCount, InactiveCount) ~ Condition * Time + (1|Subject),
                        data = Exp8_data_long,
                        family = "binomial",
                        na.action = na.omit)

summary(learning_model)

#get overall effects for learning model

learning_model_output <- car::Anova(learning_model, type = "III")
print(learning_model_output)


# getting estimated means and comparisons for within and between groups
#Within-group comparisons
conditioning_within_groups_comparisons <- emmeans(learning_model,
                                                  pairwise ~ Time | Condition,
                                                  adjust = "bonferroni",
                                                  type = "response")

#between-group comparisons
conditioning_between_groups_comparisons <- emmeans(learning_model,
                                                   pairwise ~ Condition | Time,
                                                   adjust = "bonferroni",
                                                   type = "response")

#printing comparisons
print(summary(conditioning_within_groups_comparisons$emmeans))
print(summary(conditioning_within_groups_comparisons$contrasts))
print(summary(conditioning_between_groups_comparisons$emmeans))
print(summary(conditioning_between_groups_comparisons$contrasts))


#=================================
#apa formatted results
#================================

format_p_value <- function(p) {
  if (p <.001) {
    return ("< .001")
  }  else {
      return(paste0("= ", gsub("0\\.", ". ", round(p, 3))))
    }
}


#Formatting results for glmm analysis
format_conditioning_results <- function(model, anova_output) {
  Exp8_model_summary <- summary(learning_model)
  Exp8_anova_df <- as.data.frame(anova_output)
  
  formatted_results <- data.frame(
    effect = rownames(Exp8_anova_df),
    chisq = Exp8_anova_df$Chisq,
    df = Exp8_anova_df$Df,
    p.value = Exp8_anova_df$`Pr(>Chisq)`,
    apa_results = paste0(
      "χ²(", Exp8_anova_df$Df, ") = ", round(Exp8_anova_df$Chisq, 3),
      ", *p* ", sapply(Exp8_anova_df$`Pr(>Chisq)`, format_p_value) 
    )
  )
  return(formatted_results)
  }

#Format within-group emmeans comparisons

format_within_group_comparisons <- function(emmeans_object) {
  # Extract contrasts
  contrasts_df <- as.data.frame(emmeans_object$contrasts)
  
  # Extract condition more reliably, or set explicit condition labels
  # Make sure they're in the correct order based on your graph
  condition_labels <- c("Control", "Treatment")  # Assuming this order matches your data
  condition_index <- 1:length(contrasts_df$contrast)
  
  # Create formatted results dataframe with INVERTED odds ratios
  formatted_results <- data.frame(
    contrast = contrasts_df$contrast,
    condition = condition_labels[condition_index],
    # Invert the odds ratio
    odds_ratio = 1/contrasts_df$odds.ratio,  # This inverts baseline→endpoint to endpoint→baseline
    p.value = contrasts_df$p.value,
    apa_result = paste0(
      "*OR* = ", round(contrasts_df$odds.ratio, 2),  # Invert odds ratio
      ", *z* = ", round(-contrasts_df$z.ratio, 2),   # Flip sign of z-ratio
      ", *p* ", sapply(contrasts_df$p.value, format_p_value)
    )
  )
  
  return(formatted_results)
}

#format estimated means for my manuscript

format_estimated_means <- function(emmeans_object) {
  exp8_means_df <- as.data.frame(emmeans_object$emmeans)
  
  formatted_results <- data.frame(
    group = rownames(exp8_means_df),
    condition = exp8_means_df$Condition,
    time = exp8_means_df$Time,
    probablity = exp8_means_df$prob,
    se = exp8_means_df$SE,
    lcl = exp8_means_df$asymp.LCL,
    ucl = exp8_means_df$asymp.UCL,
    apa_results = paste0(
      "*M* = ", round(exp8_means_df$prob, 2),
      ", *SE* = ", round(exp8_means_df$SE, 2),
      ", 95% CI [", round(exp8_means_df$asymp.LCL, 2), ", ", round(exp8_means_df$asymp.UCL, 2), "]"
    )
  )
  return(formatted_results)
}


# Format emmeans pairwise comparisons for between-group effects
# Format emmeans pairwise comparisons for between-group effects
format_between_group_comparisons <- function(emmeans_object) {
  # Extract contrasts
  contrasts_df <- as.data.frame(emmeans_object$contrasts)
  
  # Use explicit time point labels
  time_labels <- c("Baseline", "Endpoint")
  time_index <- 1:length(contrasts_df$contrast)
  
  # Create formatted results dataframe
  formatted_results <- data.frame(
    contrast = contrasts_df$contrast,
    time_point = time_labels[time_index],  # Use explicit labels
    odds_ratio = contrasts_df$odds.ratio,
    p.value = contrasts_df$p.value,
    apa_result = paste0(
      "*OR* = ", round(contrasts_df$odds.ratio, 2),
      ", *z* = ", round(contrasts_df$z.ratio, 2),
      ", *p* ", sapply(contrasts_df$p.value, format_p_value)
    )
  )
  
  return(formatted_results)
}

# Formating coefficients from our comparison model

format_model_coefficients <- function(model) {
  coef_table <- as.data.frame(summary(model)$coefficients)
  
  formatted_results <- data.frame(
    term = rownames(coef_table),
    estimate = coef_table$Estimate,
    se = coef_table$`Std. Error`,
    z_value = coef_table$`z value`,
    p.value = coef_table$`Pr(>|z|)`,
    apa_result = paste0(
      "β = ", round(coef_table$Estimate, 2),
      ", *SE* = ", round(coef_table$`Std. Error`, 2),
      ", *z* = ", round(coef_table$`z value`, 2),
      ", *p* = ", sapply(coef_table$`Pr(>|z|)`, format_p_value)
    )
  )
  
  return(formatted_results)
}


#primary function for inline code

format_all_glmm_results <- function(model, anova_output, within_comparisons, between_comparisons) {
  results <- list()
  
  results$anova <- format_conditioning_results(model, anova_output)
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
      # Direct match with condition label
      match_idx <- which(results$within_comparisons$condition == condition)
      
      if (length(match_idx) > 0) {
        return(results$within_comparisons$apa_result[match_idx[1]])
      } else {
        return(paste0("(no match found for condition '", condition, "')"))
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
      # Direct match with time label
      match_idx <- which(results$between_comparisons$time_point == time)
      
      if (length(match_idx) > 0) {
        return(results$between_comparisons$apa_result[match_idx[1]])
      } else {
        return(paste0("(no match found for time '", time, "')"))
      }
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
  
  # Return informative message if no match
  return("(function type not recognized)")
}

# Apply the formatting to your existing analysis
glmm_results <- format_all_glmm_results(
  learning_model,
  learning_model_output,
  conditioning_within_groups_comparisons,
  conditioning_between_groups_comparisons
)




####################################################################
#Plotting baseline to endpoint comparison
####################################################################
Exp8_Baseline_endpoint_comparison <- Exp8_data_long %>%
  group_by(Condition, Time) %>%
  summarise(
    mean_proportion = mean(ActiveArmProportion, na.rm = TRUE),
    se = sd(ActiveArmProportion, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>% 
  ggplot(aes(x = Time, y = mean_proportion, fill = Condition))+
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean_proportion - se, ymax = mean_proportion + se),
                position = position_dodge(0.7), width = 0.2, linewidth = 0.7) +
  geom_signif( #for treatment group
    annotations = "*",
    xmin = 1.2, xmax = 2.2,
    y_position = 0.85,
    color = "black",
    size = 0.6,
    textsize = 6 ) +
  geom_signif( #  for control group
    annotations = "***",
    xmin = 0.8, xmax = 1.8,
    y_position = 0.7,
    color = "black",
    size = 0.6,
    textsize = 6 ) +
  scale_fill_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  labs( 
    title = "Active Arm Entries Before and After Conditioning",
    y = "Proportion of Active Arm Entries",
    x = "Time point",
    fill = "Condition") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = (seq(0, 1, 0.1)
    )) +
  consistant_theme()




print(Exp8_Baseline_endpoint_comparison)


#=================================================================
# PART 2: LEARNING ACROSS CONDITIONING DAYS
#=================================================================

Exp8_data_long_days <- Exp8_full_data %>%
  select(Subject, Condition,
         Baseline_day1,
         Baseline_day2,
         Conditioning_day1,
         Conditioning_day2,
         Conditioning_day3, 
         Conditioning_day4) %>%
  pivot_longer(
    cols = c(Baseline_day1:Conditioning_day4),
    names_to = "TimePoint",
    values_to = "ActiveArmChoices"
  ) %>%
  mutate(
    TImePoint = factor(TimePoint,
                       levels = c("Baseline_day1", "Baseline_day2", 
                                  "Conditioning_day1", "Conditioning_day2", 
                                  "Conditioning_day3", "Conditioning_day4")),
    Subject = factor(Subject),
    Condition = factor(Condition),
    Proportion = ActiveArmChoices / 3
  )

#reshaping data for plotting

Exp8_conditioning_plot <- Exp8_data_long_days %>%
  group_by(Condition, TimePoint) %>%
  summarise(
    mean_proportion = mean(Proportion, na.rm = TRUE),
    se = sd(Proportion, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = TimePoint, y = mean_proportion, color = Condition, group = Condition)) +
  geom_line(linewidth = 1.5)+
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_proportion - se, ymax = mean_proportion + se),
                width = 0.2, linewidth = 0.8) +
  scale_color_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  # Add this scale_x_discrete function to rename the x-axis labels
  scale_x_discrete(
    labels = c(
      "Baseline_day1" = "BL1", 
      "Baseline_day2" = "BL2", 
      "Conditioning_day1" = "CD1", 
      "Conditioning_day2" = "CD2", 
      "Conditioning_day3" = "CD3", 
      "Conditioning_day4" = "CD4"
    )
  ) +
  labs(
    title = "Active Arm Entries Throughout Conditioning",
    y = "Proportion of Active Arm Entries",
    x = "Day"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  consistant_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none")

print(Exp8_conditioning_plot)

#=================================================================
# PART 3: Combinig figures for presentation
#=================================================================

Exp8_combined_figure <- Exp8_conditioning_plot / Exp8_Baseline_endpoint_comparison +
  plot_layout(
    widths = c(1, 1), 
    guides = "collect",
    ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 16, family = "Times New Roman"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
  ) 

print(Exp8_combined_figure)

ggsave("Figures/Exp8_combined_figure.png", Exp8_combined_figure, 
       width = 14, 
       height = 7, 
       dpi = 300)

#=================================================================
# PART 4: Descriptives
#=================================================================

Exp8_prefered_arm_count <- Exp8_full_data %>% count(Prefered_arm)

Exp8_active_arm_count <- Exp8_full_data %>% count(`Active arm`)

Exp8_left_active_arm_count = Exp8_active_arm_count %>% filter(`Active arm` == "L") %>% pull(n)

Exp8_right_active_arm_count = Exp8_active_arm_count %>% filter(`Active arm` == "R") %>% pull(n)

print(Exp8_left_active_arm_count)
print(Exp8_right_active_arm_count)



