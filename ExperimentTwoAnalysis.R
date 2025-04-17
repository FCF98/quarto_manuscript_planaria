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

# 1. Control Group Plot
control_plot <- summary_data_w_sample %>%
  filter(Condition == "Control") %>%
  ggplot(aes(x = Time, y = mean_proportion, group = 1)) +
  geom_line(linewidth = 1.5, color = "#159090") +
  geom_point(size = 4, color = "#159090") +
  geom_errorbar(aes(ymin = mean_proportion - se_proportion, 
                    ymax = mean_proportion + se_proportion), 
                width = 0.2, linewidth = 1, color = "#159090") +
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
    y = "Mean Proportion of Active Arm Choices"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 14),  # Reduced from 16
    axis.title = element_text(size = 16, face = "bold"),  # Reduced from 18
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.text = element_text(size = 12, color = "black"),  # Reduced from 14
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),  # Reduced from 16
    legend.text = element_text(size = 12),  # Reduced from 14
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


# 2. Treatment Group Plot
treatment_plot <- summary_data_w_sample %>%
  filter(Condition == "Treatment") %>%
  ggplot(aes(x = Time, y = mean_proportion, group = 1)) +
  geom_line(linewidth = 1.5, color = "#FF8C00") +
  geom_point(size = 4, color = "#FF8C00") +
  geom_errorbar(aes(ymin = mean_proportion - se_proportion, 
                    ymax = mean_proportion + se_proportion), 
                width = 0.2, linewidth = 1, color = "#FF8C00") +
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
    y_position = 0.8,
    color = "black",
    size = 0.6,
    textsize = 6
  ) +
  labs(
    title = "Treatment Group",
    x = "Time Point",
    y = "Mean Proportion of Active Arm Choices"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 14),  # Reduced from 16
    axis.title = element_text(size = 16, face = "bold"),  # Reduced from 18
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.text = element_text(size = 12, color = "black"),  # Reduced from 14
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),  # Reduced from 16
    legend.text = element_text(size = 12),  # Reduced from 14
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

# 3. Combined Plot (modify your original plot to show only between-group differences)
combined_plot <- ggplot(summary_data_w_sample, 
                        aes(x = Time, y = mean_proportion, 
                            color = Condition, group = Condition)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_proportion - se_proportion, 
                    ymax = mean_proportion + se_proportion), 
                width = 0.1, linewidth = 1) +
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
    y = "Mean Proportion of Active Arm Choices",
    color = "Condition"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 14),  # Reduced from 16
    axis.title = element_text(size = 16, face = "bold"),  # Reduced from 18
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.text = element_text(size = 12, color = "black"),  # Reduced from 14
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),  # Reduced from 16
    legend.text = element_text(size = 12),  # Reduced from 14
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

# Save the panel
ggsave("panel_plot.png", combined_panel, 
       width = 12, height = 12, dpi = 300) 




############ Creating individual plots  #########


# First, create an ordered factor for subjects based on their condition
data_long <- data_long %>%
  mutate(Subject_ordered = factor(Subject, 
                                  levels = data_long %>%
                                    select(Subject, Condition) %>%
                                    distinct() %>%
                                    arrange(desc(Condition == "Control"), Subject) %>%
                                    pull(Subject)))

# Then modify your plot to use the new ordered factor
individual_plots <- ggplot(data_long, aes(x = Time, y = ActiveArmProportion, group = 1, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~Subject_ordered, ncol = 6) +  # Use Subject_ordered instead of Subject
  labs(
    title = "Individual Subject Active Arm Choices Over Time",
    x = "Time Point",
    y = "Proportion of Active Arm Choices"
  ) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    panel.border = element_rect(color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = c("Control" = "#159090", "Treatment" = "#FF8C00"))

individual_plots

# Save the plot
ggsave("individual_subject_plots.png", individual_plots, 
       width = 20, height = 24, dpi = 300, units = "in",
       bg = "white")






###### lollipop chart ######


# Prepare the data for the lollipop chart
lollipop_chart_data <- data_long %>%
  filter(Time %in% c("Baseline", "Endpoint")) %>%
  group_by(Subject, Condition) %>%
  summarize(change = ActiveArmProportion[Time == "Endpoint"] - ActiveArmProportion[Time == "Baseline"],
            .groups = "drop")

# Create the lollipop chart
ggplot(lollipop_chart_data, aes(x = as.factor(Subject), y = change, color = Condition)) +
  geom_segment(aes(x = as.factor(Subject), xend = as.factor(Subject), y = 0, yend = change), 
               color = "grey") +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Treatment" = rgb(0.2,0.7,0.1,0.5), "Control" = rgb(0.7,0.2,0.1,0.5))) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = "Subject",
    y = "Change in proportion for active arm (Endpoint - Baseline)",
    title = "Change active arm preference from Baseline to Endpoint",
    subtitle = "Positive values indicate an increase, negative values indicate a decrease"
  )



######### Plotting learning across conditioning days ##########


Exp2_fulldata <- read_excel("Datasets/ExperimentTwoFullData.xlsx", sheet = "Data for analysis")

Exp2_data_w_exclusions <- Exp2_fulldata %>%
  filter(!Subject %in% c(3, 4, 14, 22, 27, 35, 40, 44, 47, 50, 52, 55, 56))


Data_long_days <- Exp2_data_w_exclusions %>%
  select(Subject, Condition,
         Baseline = `Baseline%`,
         Day1 = `Day1%`,
         Day2 = `Day2%`,
         Day3 = `Day3%`,
         Day4 = `Day4%`,
         Day5 = `Day5%`) %>%
  pivot_longer(
    cols = c(Baseline:Day5),
    names_to = "TimePoint",
    values_to = "Proportion"
  ) %>%
  mutate(
    TimePoint = factor(TimePoint,
                       levels = c("Baseline", "Day1", "Day2", "Day3", "Day4", "Day5")),
    Condition = factor(Condition)
  ) %>%
  # Calculate mean and SE for each condition and timepoint
  group_by(Condition, TimePoint) %>%
  summarise(
    mean_prop = mean(Proportion, na.rm = TRUE),
    se = sd(Proportion, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Create the visualization
ggplot(Data_long_days, aes(x = TimePoint, y = mean_prop,
                      color = Condition, group = Condition)) +
  # Add mean lines
  geom_line(linewidth = 1.5) +
  # Add mean points
  geom_point(size = 3) +
  # Add error bars
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
                width = 0.2) +
  # Customize the theme and labels
  theme_minimal() +
  labs(
    title = "Average Proportion of Active Arm Choices Over Time",
    y = "Proportion of Active Arm Choices",
    x = "Time Point"
  ) +
  scale_color_manual(values = c("Control" = "#0072B2", "Treatment" = "#D55E00")) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylim(0, 0.75)  # Set y-axis limits from 0 to 1



###############################################
# Descriptives
###############################################

Exp2_data_w_exclusions <- Exp2_fulldata %>%
  filter(!Subject %in% c(3, 4, 14, 22, 27, 35, 40, 44, 47, 50, 52, 55, 56))

#Pulling baseline arm preference descriptives

Exp2_prefered_arm_count <- Exp2_data_w_exclusions %>% filter(Subject %in% (1:60)) %>% count(Prefered_arm)

Exp2_active_arm_count <- Exp2_data_w_exclusions %>% filter(Subject %in% (1:60)) %>% count(active_arm)

Exp2_left_active_arm_count <- Exp2_active_arm_count %>% filter(active_arm == "L") %>% pull(n)

Exp2_right_active_arm_count <- Exp2_active_arm_count %>% filter(active_arm == "R") %>% pull(n)

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

print(active_arm_sample_sizes)


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

# Print counts for each comparison
print("Subjects available for Baseline vs Test comparison:")
print(baseline_test_subjects)

print("Subjects available for Baseline vs Endpoint comparison:")
print(baseline_endpoint_subjects)

print("Subjects available for Endpoint vs Test comparison:")
print(endpoint_test_subjects)

print("Subjects available for Test vs Reinstatement comparison:")
print(test_reinstatement_subjects)

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
print(summary(ref_grid))

# Check samples used in each emmeans comparison
for (pair in time_pairs) {
  emm_pair <- emmeans(m1, specs = "Time", at = list(Time = pair))
  # Attempt to get subjects used in this comparison
  subjects_used <- tryCatch({
    recover_data(m1, specs = emm_pair)$Subject
  }, error = function(e) {
    # Handle the case where recover_data might not work as expected
    return(NULL)
  })
  
  if (!is.null(subjects_used)) {
    unique_subjects <- unique(subjects_used)
    cat(sprintf("\nFor %s vs %s comparison, %d unique subjects were used\n", 
                pair[1], pair[2], length(unique_subjects)))
  } else {
    cat(sprintf("\nCouldn't extract subjects for %s vs %s comparison\n", 
                pair[1], pair[2]))
  }
}
