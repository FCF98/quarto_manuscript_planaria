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

data <- read_excel("Datasets/ExperimentTwoSummaryData.xlsx")



#########structuring data correctly############


data_long <- data %>%
  # First ensure all percentage columns are numeric
  mutate(across(ends_with("_arm_%"), 
                ~as.numeric(as.character(.)))) %>%
  pivot_longer(
    cols = c(`baseline_active_arm_%`, `endpoint_active_arm_%`, `test_active_arm_%`, `reinstatement_active_arm_%`),
    names_to = "Time",
    values_to = "ActiveCount"
  ) %>%
  mutate(
    Time = factor(Time, 
                  levels = c("baseline_active_arm_%", "endpoint_active_arm_%", "test_active_arm_%", "reinstatement_active_arm_%"),
                  labels = c("Baseline", "Endpoint", "Test", "Reinstatement")),
    Condition = factor(Condition),
    active_arm = factor(active_arm),
    Subject = factor(Subject),
    # Convert to integer only for non-NA values
    ActiveCount = if_else(!is.na(ActiveCount), 
                          as.integer(ActiveCount), 
                          NA_integer_)
  ) %>%
  mutate(
    TotalTrials = case_when(
      Time %in% c("Baseline", "Endpoint") ~ 6,
      Time %in% c("Test", "Reinstatement") ~ 4
    ),
    # Handle NAs in calculations
    InactiveCount = if_else(!is.na(ActiveCount),
                            TotalTrials - ActiveCount,
                            NA_integer_),
    ActiveArmProportion = if_else(!is.na(ActiveCount),
                                  ActiveCount / TotalTrials,
                                  NA_real_)
  ) %>%
  # Only filter out impossible data points where ActiveCount is not NA
  dplyr::filter(is.na(ActiveCount) | ActiveCount <= TotalTrials)




###### model to test for main effects and interactions across 4 time points #####

contrasts(data_long$Condition) <- contr.sum

contrasts(data_long$Time) <- contr.sum

m1 <- glmer(cbind(ActiveCount, InactiveCount) ~ Condition * Time + (1|Subject), 
            data = data_long, 
            family = "binomial",
            na.action = na.omit)


summary(m1)

exp2_decisions_model_output <- car::Anova(m1, type = "III")

exp2_decision_comparisons_within_group <-  emmeans(m1,pairwise~Time|Condition,adjust="bonferroni",type="response")

exp2_decision_comparisons_between_group <- emmeans(m1,pairwise~Condition|Time,adjust="bonferroni",type="response")

exp2_decision_comparisons_within_group <- summary(exp2_decision_comparisons_within_group)$contrasts
 
exp2_decision_comparisons_between_group <- summary(exp2_decision_comparisons_between_group)$contrasts


##### calculating cohens h values for effect size reporting #####

### Getting the within group effect sizes

# Function for Cohen's h calculation
calc_cohens_h <- function(p1, p2) {
  2 * (asin(sqrt(p1)) - asin(sqrt(p2)))
}

# Get mean proportions for each time point
time_props <- data_long %>%
  group_by(Time, Condition) %>%
  summarize(
    prop = mean(ActiveArmProportion, na.rm = TRUE),
    .groups = "drop"
  )

# Define contrast levels explicitly
contrast_levels <- c(
  "Baseline / Endpoint",
  "Baseline / Test",
  "Baseline / Reinstatement",
  "Endpoint / Test",
  "Endpoint / Reinstatement",
  "Test / Reinstatement"
)

# Get p-values from emmeans results
p_values_df <- exp2_decision_comparisons_within_group %>%
  as.data.frame() %>%
  select(Condition, contrast, p.value)

# Get all unique times
times <- levels(data_long$Time)
# Create all possible pairs
pairs <- t(combn(times, 2))

# Create the final dataframe with all comparisons
within_group_h <- do.call(rbind, lapply(unique(data_long$Condition), function(cond) {
  # Get props for this condition
  condition_props <- time_props %>% filter(Condition == cond)
  
  # Create dataframe of all comparisons
  data.frame(
    Condition = cond,
    contrast = factor(
      ifelse(grepl("Reinstatement", paste(pairs[,1], "/", pairs[,2])),
             gsub("Reinstatement", "Reinstatement", paste(pairs[,1], "/", pairs[,2])),
             paste(pairs[,1], "/", pairs[,2])),
      levels = contrast_levels
    ),
    prop1 = condition_props$prop[match(pairs[,1], condition_props$Time)],
    prop2 = condition_props$prop[match(pairs[,2], condition_props$Time)],
    stringsAsFactors = FALSE
  )
})) %>%
  # Calculate Cohen's h
  mutate(
    cohens_h = abs(calc_cohens_h(prop1, prop2))
  ) %>%
  # Join with p-values
  left_join(p_values_df, by = c("Condition", "contrast")) %>%
  # Format results
  mutate(
    apa_result = paste0("*h* = ", round(cohens_h, 2), 
                        ", *p* ", ifelse(p.value < .001, "< .001",
                                       paste0("= ",gsub("0\\.", ".", round(p.value, 3)))))
  )


# Create separate dataframes for control and treatment results
control_h_values_exp2_decisions <- within_group_h %>%
  filter(Condition == "Control") %>%
  select(contrast, cohens_h, p.value, apa_result)

treatment_h_values_exp2_decisions <- within_group_h %>%
  filter(Condition == "Treatment") %>%
  select(contrast, cohens_h, p.value, apa_result)



### Getting the between group effect sizes


p_values_between_df <- exp2_decision_comparisons_between_group %>%
  as.data.frame()

# Calculate between-group h values
between_group_h <- time_props %>%
  pivot_wider(
    names_from = Condition,
    values_from = prop
  ) %>%
  mutate(
    cohens_h = abs(calc_cohens_h(Control, Treatment))
  )

# Create final between-group results dataframe
between_group_h_final <- between_group_h %>%
  mutate(
    cohens_h = abs(calc_cohens_h(Control, Treatment))
  ) %>%
  left_join(
    p_values_between_df %>% 
      select(Time, p.value),
    by = "Time"
  ) %>%
  mutate(
    contrast = factor(Time, levels = times),  # Use the times object we defined earlier
    apa_result = paste0("*h* = ", round(cohens_h, 2), 
                        ", *p* ", ifelse(p.value < .001, "< .001",
                                         paste0("= ", round(p.value, 3))))
  ) %>%
  select(Time, contrast, cohens_h, p.value, apa_result)


############# plotting the grouped results   ############


# First, create a data frame for the between-group significance annotations
between_signif <- data.frame(
  x = 3,  # Test time point
  y = 0.45,  # Position between the two points
  label = "**"  # Two asterisks for p < 0.01
)

# Create the summary data for graphing

summary_data_w_sample <- data_long %>%
  group_by(Condition, Time) %>%
  summarise(
    mean_proportion = mean(ActiveArmProportion, na.rm = TRUE),
    se_proportion = sd(ActiveArmProportion, na.rm = TRUE) / sqrt(sum(!is.na(ActiveArmProportion))),
    n = sum(!is.na(ActiveArmProportion)),
    .groups = 'drop'
  ) %>%
  mutate(label = sprintf("n=%d", n))



# Create APA-style plots, and arranging in a panel of 3 plots

# 1. Control Group Plot with colored sample sizes
control_plot <- summary_data_w_sample %>%
  filter(Condition == "Control") %>%
  ggplot(aes(x = Time, y = mean_proportion, group = 1)) +
  geom_line(linewidth = 1.5, color = "#159090") +
  geom_point(size = 4, color = "#159090") +
  geom_errorbar(aes(ymin = mean_proportion - se_proportion, 
                    ymax = mean_proportion + se_proportion), 
                width = 0.2, linewidth = 1, color = "#159090") +
  # Add sample size labels just above the x-axis with matching color
  geom_text(aes(label = label, y = 0.03), color = "#159090", size = 4.5) +
  # Control group within comparisons
  geom_signif(
    annotations = "***",
    xmin = 1, xmax = 2,
    y_position = 0.63,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  geom_signif(
    annotations = "***",
    xmin = 2, xmax = 3,
    y_position = 0.70,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  labs(
    title = "Control Group",
    x = "Time Point",
    y = "Mean Proportion of Active Arm Entries"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 16)),
    axis.text = element_text(size = 16, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(t = 30, r = 30, b = 50, l = 30, unit = "pt"),
    panel.border = element_blank()
  )+
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
  ) +
  scale_color_manual(values = c("Control" = "#159090", "Treatment" = "#FF8C00"))

# 2. Treatment Group Plot with colored sample sizes
treatment_plot <- summary_data_w_sample %>%
  filter(Condition == "Treatment") %>%
  ggplot(aes(x = Time, y = mean_proportion, group = 1)) +
  geom_line(linewidth = 1.5, color = "#FF8C00") +
  geom_point(size = 4, color = "#FF8C00") +
  geom_errorbar(aes(ymin = mean_proportion - se_proportion, 
                    ymax = mean_proportion + se_proportion), 
                width = 0.2, linewidth = 1, color = "#FF8C00") +
  # Add sample size labels just above the x-axis with matching color
  geom_text(aes(label = label, y = 0.03), color = "#FF8C00", size = 4.5) +
  # Treatment group within comparisons
  geom_signif(
    annotations = "**",
    xmin = 1, xmax = 2,
    y_position = 0.6,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  geom_signif(
    annotations = "*",
    xmin = 1, xmax = 3,
    y_position = 0.7,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  geom_signif(
    annotations = "*",
    xmin = 2, xmax = 4,
    y_position = 0.85,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  labs(
    title = "Treatment Group",
    x = "Time Point",
    y = "Mean Proportion of Active Arm Enties"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 16)),
    axis.text = element_text(size = 16, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(t = 30, r = 30, b = 50, l = 30, unit = "pt"),
    panel.border = element_blank()
  )+
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
  ) +
  scale_color_manual(values = c("Control" = "#159090", "Treatment" = "#FF8C00"))

# 3. Combined Plot with closer-spaced sample sizes
# Create a data frame for label positioning in the combined plot with smaller offsets
combined_labels <- summary_data_w_sample %>%
  # Calculate small x-offsets to keep labels closer together
  mutate(
    x_pos = as.numeric(Time) + ifelse(Condition == "Control", -0.1, 0.1),
    y_pos = 0.03
  )

combined_plot <- ggplot(summary_data_w_sample, 
                        aes(x = Time, y = mean_proportion, 
                            color = Condition, group = Condition)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_proportion - se_proportion, 
                    ymax = mean_proportion + se_proportion), 
                width = 0.1, linewidth = 1) +
  # Add sample size labels for each condition at each time point with smaller offset
  geom_text(
    data = combined_labels,
    aes(x = x_pos, y = y_pos, label = label, color = Condition),
    size = 4.5
  ) +
  # Only between group comparison
  annotate("text", 
           x = 3, 
           y = 0.53, 
           label = "*", 
           size = 6,
           color = "black")+
  labs(
    title = "Between-Group Comparison",
    x = "Time Point",
    y = "Mean Proportion of Active Arm Entries",
    color = "Condition"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 16)),
    axis.text = element_text(size = 16, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(t = 30, r = 30, b = 50, l = 30, unit = "pt"),
    panel.border = element_blank()
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
  ) +
  scale_color_manual(values = c("Control" = "#159090", "Treatment" = "#FF8C00"))

# Arrange the plots in a panel
combined_panel <- (control_plot + treatment_plot) / (combined_plot) +
  plot_layout(heights = c(1.2, 1.2)) +
  plot_annotation(tag_levels = 'A')

# Display the panel
combined_panel

# Save the panel
ggsave("panel_plot.png", combined_panel, 
       width = 12, height = 12, dpi = 300) 



###############################################
# Descriptives
###############################################

Exp2_fulldata <- read_excel("Datasets/ExperimentTwoFullData.xlsx", sheet = "Data for analysis")

Exp2_data_w_exclusions <- Exp2_fulldata %>%
  filter(!Subject %in% c(3, 4, 14, 22, 27, 35, 40, 44, 47, 50, 52, 55, 56))

#Pulling baseline arm preference descriptives

Exp2_prefered_arm_count <- Exp2_data_w_exclusions %>% filter(Subject %in% (1:60)) %>% count(Prefered_arm)

Exp2_active_arm_count <- Exp2_data_w_exclusions %>% filter(Subject %in% (1:60)) %>% count(active_arm)

Exp2_left_active_arm_count <- Exp2_active_arm_count %>% filter(active_arm == "L") %>% pull(n)

Exp2_right_active_arm_count <- Exp2_active_arm_count %>% filter(active_arm == "R") %>% pull(n)

# Calculate summary statistics for plotting
Exp2_summary_Stats <- data_long %>% filter(Time == "Test") %>%
  group_by(Condition) %>%
  summarise(
    mean_prop = mean(ActiveArmProportion, na.rm = TRUE),
    sd = sd(ActiveArmProportion, na.rm = TRUE)
  )

############## descriptives for excluded subjects ######################

Exp2_excluded_subjects_data <- Exp2_fulldata %>%
  filter(Subject %in% c(3, 4, 14, 22, 27, 35, 40, 44, 47, 50, 52, 55, 56))

Exp2_excluded_subjects_prefered_arm <- Exp2_excluded_subjects_data %>% count(Prefered_arm)

############ checking subjects avaliable for between groups comparisons

active_arm_sample_sizes <- data_long %>%
  group_by(Condition, Time) %>%
  summarise(
    n = sum(!is.na(ActiveArmProportion)),
    .groups = 'drop'
  )


############ checking subjects available for within subjects comparisons

# Create a data frame showing which subjects have data at each timepoint
active_arm_availability <- data_long %>%
  group_by(Subject, Condition, Time) %>%
  summarize(has_data = !is.na(ActiveArmProportion), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(Subject, Condition),
    names_from = Time,
    values_from = has_data
  )

# Count subjects with data at specific pairs of timepoints
baseline_test_subjects <- active_arm_availability %>%
  filter(Baseline == TRUE & Test == TRUE) %>%
  group_by(Condition) %>%
  summarize(count = n())

baseline_endpoint_subjects <- active_arm_availability %>%
  filter(Baseline == TRUE & Endpoint == TRUE) %>%
  group_by(Condition) %>%
  summarize(count = n())

endpoint_test_subjects <- active_arm_availability %>%
  filter(Endpoint == TRUE & Test == TRUE) %>%
  group_by(Condition) %>%
  summarize(count = n())

test_reinstatement_subjects <- active_arm_availability %>%
  filter(Test == TRUE & Reinstatement == TRUE) %>%
  group_by(Condition) %>%
  summarize(count = n())


# For the most accurate view of which subjects were used in the GLMER model's
# pairwise comparisons via emmeans:
# Define the pairs you want to check
time_pairs <- list(
  c("Baseline", "Endpoint"),
  c("Baseline", "Test"),
  c("Baseline", "Reinstatement"),
  c("Endpoint", "Test"),
  c("Endpoint", "Reinstatement"),
  c("Test", "Reinstatement")
)

# Get reference grid from the model
ref_grid <- ref_grid(m1)

