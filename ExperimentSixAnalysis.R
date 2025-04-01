# R Code for Planaria CPP Analysis with Movement Filtering
# This analysis explores whether drug conditioning affects planaria preference,
# and whether this memory persists in regenerated head and tail segments

# Load required libraries
library(tidyverse)
library(lme4)
library(emmeans)
library(ggplot2)
library(effsize)
library(patchwork) # Added for combining plots

# Read in your original data file
data_original <- read.csv("Datasets/ExperimentSixFullData.csv")

# First, filter out subjects with low movement during regeneration test
# Calculate group medians for regenerate movement by condition and body part
print("Calculating movement thresholds...")

# Identify whole animals vs head/tail regenerates
data_original$is_whole <- !grepl("[HT]", data_original$Subject)
data_original$body_part <- "Whole"
data_original$body_part[grepl("H", data_original$Subject)] <- "Head"
data_original$body_part[grepl("T", data_original$Subject)] <- "Tail"

# Calculate group medians for movement during regeneration test
group_medians <- data_original %>%
  filter(!is_whole) %>%  # Only regenerates have distance data
  group_by(Condition, body_part) %>%
  summarise(median_distance = median(Regeneration_test_distance, na.rm = TRUE),
            .groups = "drop")

# Calculate the 25% threshold for each group
group_medians$threshold <- group_medians$median_distance * 0.25

# Print the thresholds
print("Movement thresholds (25% of group medians):")
print(group_medians)

# Identify subjects with movement less than 25% of group median
to_exclude <- data_original %>%
  filter(!is_whole) %>%
  left_join(group_medians, by = c("Condition", "body_part")) %>%
  filter(Regeneration_test_distance < threshold) %>%
  pull(Subject)

print(paste("Number of excluded subjects due to low movement:", length(to_exclude)))
if(length(to_exclude) > 0) {
  print("Excluded subjects:")
  print(to_exclude)
}

# Filter the dataset to remove excluded subjects
data_filtered <- data_original %>%
  filter(!(Subject %in% to_exclude))

print(paste("Original dataset had", nrow(data_original), "rows"))
print(paste("After filtering, dataset has", nrow(data_filtered), "rows"))

# Continue with the analysis using data_filtered instead of data_original
# Extract original subject number for regenerates
data_filtered$original_subject <- data_filtered$Subject
data_filtered$original_subject[!data_filtered$is_whole] <- 
  as.numeric(gsub("[HT]", "", data_filtered$Subject[!data_filtered$is_whole]))

# Create a long format dataset for all phases
create_analysis_data <- function() {
  # Create empty dataframe
  analysis_data <- data.frame()
  
  # Process whole animals
  whole_animals <- data_filtered[data_filtered$is_whole, ]
  
  # Add baseline data for whole animals
  baseline_whole <- data.frame(
    Subject = whole_animals$Subject,
    OriginalSubject = whole_animals$Subject,
    BodyPart = "Whole",
    Phase = "Baseline",
    Condition = whole_animals$Condition,
    ActiveSurface = whole_animals$Active_surface,
    PropActiveArm = whole_animals$Baseline_active_proportion
  )
  analysis_data <- rbind(analysis_data, baseline_whole)
  
  # Add test data for whole animals
  test_whole <- data.frame(
    Subject = whole_animals$Subject,
    OriginalSubject = whole_animals$Subject,
    BodyPart = "Whole",
    Phase = "Test",
    Condition = whole_animals$Condition,
    ActiveSurface = whole_animals$Active_surface,
    PropActiveArm = whole_animals$Test_active_proportion
  )
  analysis_data <- rbind(analysis_data, test_whole)
  
  # Process head regenerates
  head_regenerates <- data_filtered[data_filtered$body_part == "Head", ]
  
  # Add baseline data for heads (using baseline from whole animal)
  baseline_head <- data.frame(
    Subject = head_regenerates$Subject,
    OriginalSubject = head_regenerates$original_subject,
    BodyPart = "Head",
    Phase = "Baseline",
    Condition = head_regenerates$Condition,
    ActiveSurface = head_regenerates$Active_surface,
    PropActiveArm = NA
  )
  
  # Match baseline values from whole animals
  for (i in 1:nrow(baseline_head)) {
    original_id <- baseline_head$OriginalSubject[i]
    original_row <- which(whole_animals$Subject == original_id)
    if (length(original_row) > 0) {
      baseline_head$PropActiveArm[i] <- whole_animals$Baseline_active_proportion[original_row]
    }
  }
  
  analysis_data <- rbind(analysis_data, baseline_head)
  
  # Add regeneration test data for heads
  regen_head <- data.frame(
    Subject = head_regenerates$Subject,
    OriginalSubject = head_regenerates$original_subject,
    BodyPart = "Head",
    Phase = "Regeneration",
    Condition = head_regenerates$Condition,
    ActiveSurface = head_regenerates$Active_surface,
    PropActiveArm = head_regenerates$Regeneration_test_active_proportion
  )
  analysis_data <- rbind(analysis_data, regen_head)
  
  # Add reinstatement data for heads if available
  if (any(!is.na(head_regenerates$Reinstatement_test_active_proportion))) {
    reinst_head <- data.frame(
      Subject = head_regenerates$Subject,
      OriginalSubject = head_regenerates$original_subject,
      BodyPart = "Head",
      Phase = "Reinstatement",
      Condition = head_regenerates$Condition,
      ActiveSurface = head_regenerates$Active_surface,
      PropActiveArm = head_regenerates$Reinstatement_test_active_proportion
    )
    analysis_data <- rbind(analysis_data, reinst_head)
  }
  
  # Process tail regenerates
  tail_regenerates <- data_filtered[data_filtered$body_part == "Tail", ]
  
  # Add baseline data for tails (using baseline from whole animal)
  baseline_tail <- data.frame(
    Subject = tail_regenerates$Subject,
    OriginalSubject = tail_regenerates$original_subject,
    BodyPart = "Tail",
    Phase = "Baseline",
    Condition = tail_regenerates$Condition,
    ActiveSurface = tail_regenerates$Active_surface,
    PropActiveArm = NA
  )
  
  # Match baseline values from whole animals
  for (i in 1:nrow(baseline_tail)) {
    original_id <- baseline_tail$OriginalSubject[i]
    original_row <- which(whole_animals$Subject == original_id)
    if (length(original_row) > 0) {
      baseline_tail$PropActiveArm[i] <- whole_animals$Baseline_active_proportion[original_row]
    }
  }
  
  analysis_data <- rbind(analysis_data, baseline_tail)
  
  # Add regeneration test data for tails
  regen_tail <- data.frame(
    Subject = tail_regenerates$Subject,
    OriginalSubject = tail_regenerates$original_subject,
    BodyPart = "Tail",
    Phase = "Regeneration",
    Condition = tail_regenerates$Condition,
    ActiveSurface = tail_regenerates$Active_surface,
    PropActiveArm = tail_regenerates$Regeneration_test_active_proportion
  )
  analysis_data <- rbind(analysis_data, regen_tail)
  
  # Add reinstatement data for tails if available
  if (any(!is.na(tail_regenerates$Reinstatement_test_active_proportion))) {
    reinst_tail <- data.frame(
      Subject = tail_regenerates$Subject,
      OriginalSubject = tail_regenerates$original_subject,
      BodyPart = "Tail",
      Phase = "Reinstatement",
      Condition = tail_regenerates$Condition,
      ActiveSurface = tail_regenerates$Active_surface,
      PropActiveArm = tail_regenerates$Reinstatement_test_active_proportion
    )
    analysis_data <- rbind(analysis_data, reinst_tail)
  }
  
  # Convert factors
  analysis_data$Subject <- as.factor(analysis_data$Subject)
  analysis_data$OriginalSubject <- as.factor(analysis_data$OriginalSubject)
  analysis_data$BodyPart <- as.factor(analysis_data$BodyPart)
  analysis_data$Phase <- as.factor(analysis_data$Phase)
  analysis_data$Condition <- as.factor(analysis_data$Condition)
  analysis_data$ActiveSurface <- as.factor(analysis_data$ActiveSurface)
  
  return(analysis_data)
}

# Create the analysis dataset
data <- create_analysis_data()

# Check the structure of the data
str(data)
summary(data)

# Calculate difference scores (Test - Baseline)
calculate_diff_scores <- function(data) {
  # Wide format for baseline and test phases
  baseline_data <- data %>%
    filter(Phase == "Baseline") %>%
    select(OriginalSubject, BodyPart, Condition, ActiveSurface, PropActiveArm) %>%
    rename(Baseline = PropActiveArm)
  
  # Get all post-baseline phases
  post_phases <- data %>%
    filter(Phase != "Baseline") %>%
    select(OriginalSubject, BodyPart, Phase, Condition, ActiveSurface, PropActiveArm)
  
  # Calculate differences by joining
  diff_scores <- post_phases %>%
    left_join(baseline_data, 
              by = c("OriginalSubject", "BodyPart", "Condition", "ActiveSurface")) %>%
    mutate(DiffScore = PropActiveArm - Baseline)
  
  return(diff_scores)
}

diff_scores <- calculate_diff_scores(data)

#########################
# ANALYSIS 1: Did conditioning change preference?
# Compare Baseline vs Test for Whole subjects
#########################

# Filter data for just whole subjects, baseline and test phases
whole_baseline_test <- data %>% 
  filter(BodyPart == "Whole", Phase %in% c("Baseline", "Test"))

# Descriptive statistics
whole_summary <- whole_baseline_test %>%
  group_by(Condition, Phase) %>%
  summarise(
    mean_prop = mean(PropActiveArm, na.rm = TRUE),
    sd_prop = sd(PropActiveArm, na.rm = TRUE),
    n = n(),
    se_prop = sd_prop / sqrt(n),
    .groups = "drop"
  )

print("Summary statistics for whole animals - Baseline vs Test:")
print(whole_summary)

# Create bar plot for comparing Baseline to Test (Figure A)
# Rename Phase levels to match the desired figure
whole_summary <- whole_summary %>%
  mutate(Phase = factor(Phase, levels = c("Baseline", "Test"),
                        labels = c("Baseline", "Endpoint")))

# Make sure conditions are in the desired order (Control, Cocaine, Methamphetamine)
whole_summary$Condition <- factor(whole_summary$Condition, 
                                  levels = c("Control", "Cocaine", "Methamphetamine"))

# Create bar plot
figA <- ggplot(whole_summary, aes(x = Phase, y = mean_prop, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_errorbar(aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop),
                position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c("Control" = "#159090", "Cocaine" = "#FF8C00", "Methamphetamine" = "#D81B60"),
                    name = "Condition") +
  labs(title = "Time on Active Surface Before and After Conditioning",
       y = "Proportion of Time on Active Surface",
       x = "Time") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 10),
    legend.position = "right",
    axis.title = element_text(size = 11),
    plot.title = element_text(size = 12, hjust = 0.5),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  ylim(0, 1)  # Match y-axis scale in your figure

print(figA)

# Statistical test - mixed effects model
model_whole <- lmer(PropActiveArm ~ Phase * Condition + (1|Subject), data = whole_baseline_test)
print(summary(model_whole))
print(anova(model_whole))

# Post-hoc tests
emm_whole <- emmeans(model_whole, ~ Phase | Condition)
print(pairs(emm_whole))

# Effect size - Cohen's d for each condition
whole_baseline_test %>%
  group_by(Condition) %>%
  summarize(
    d = cohen.d(PropActiveArm[Phase == "Test"], PropActiveArm[Phase == "Baseline"])$estimate,
    .groups = "drop"
  )

#########################
# ANALYSIS 2: Did regenerated parts retain conditioning?
# Compare Baseline vs Regeneration for Head and Tail
#########################

# Filter for baseline and regeneration phases
regen_data <- data %>%
  filter(Phase %in% c("Baseline", "Regeneration"))

# Descriptive statistics by body part and condition
regen_summary <- regen_data %>%
  group_by(Condition, BodyPart, Phase) %>%
  summarise(
    mean_prop = mean(PropActiveArm, na.rm = TRUE),
    sd_prop = sd(PropActiveArm, na.rm = TRUE),
    n = n(),
    se_prop = sd_prop / sqrt(n),
    .groups = "drop"
  )

print("Summary statistics for regenerated parts:")
print(regen_summary)

# Create Figure C: Regeneration comparison plot
# Filter for baseline (whole animals) and regeneration data (head and tail)
regen_comparison <- data %>%
  filter((BodyPart == "Whole" & Phase == "Baseline") | 
           (BodyPart %in% c("Head", "Tail") & Phase == "Regeneration"))

# Create a new variable for the x-axis
regen_comparison <- regen_comparison %>%
  mutate(BodyStatus = case_when(
    BodyPart == "Whole" ~ "Original",
    BodyPart == "Head" ~ "Head",
    BodyPart == "Tail" ~ "Tail"
  ))

# Calculate group means and standard errors
regen_summary <- regen_comparison %>%
  group_by(Condition, BodyStatus) %>%
  summarise(
    mean_prop = mean(PropActiveArm, na.rm = TRUE),
    se_prop = sd(PropActiveArm, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Make sure conditions are in the desired order (Control, Cocaine, Methamphetamine)
regen_summary$Condition <- factor(regen_summary$Condition, 
                                  levels = c("Control", "Cocaine", "Methamphetamine"))

# Make sure body status is in the desired order (Original, Head, Tail)
regen_summary$BodyStatus <- factor(regen_summary$BodyStatus,
                                   levels = c("Original", "Head", "Tail"))

# Create point plot
figC <- ggplot(regen_summary, aes(x = BodyStatus, y = mean_prop, color = Condition, shape = BodyStatus)) +
  facet_wrap(~ Condition, scales = "free_x") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop), width = 0.2) +
  scale_color_manual(values = c("Control" = "#159090", "Cocaine" = "#FF8C00", "Methamphetamine" = "#D81B60"),
                     name = "Condition") +
  scale_shape_manual(values = c("Original" = 16, "Head" = 17, "Tail" = 25),
                     name = "Body Part") +
  labs(title = "Retention of CPP Memory in Regenerated Tissue",
       y = "Proportion of Time on Active Surface",
       x = "Body Status") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10),
    legend.position = "right",
    axis.title = element_text(size = 11),
    plot.title = element_text(size = 12, hjust = 0.5),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  ylim(0, 1)  # Match y-axis scale in your figure

print(figC)

# Statistical test - separate models for head and tail
head_regen <- regen_data %>% filter(BodyPart == "Head")
tail_regen <- regen_data %>% filter(BodyPart == "Tail")

# Head model
model_head_regen <- lmer(PropActiveArm ~ Phase * Condition + (1|OriginalSubject), data = head_regen)
print(summary(model_head_regen))
print(anova(model_head_regen))

# Tail model
model_tail_regen <- lmer(PropActiveArm ~ Phase * Condition + (1|OriginalSubject), data = tail_regen)
print(summary(model_tail_regen))
print(anova(model_tail_regen))

# Post-hoc tests
emm_head <- emmeans(model_head_regen, ~ Phase | Condition)
print(pairs(emm_head))

emm_tail <- emmeans(model_tail_regen, ~ Phase | Condition)
print(pairs(emm_tail))

# Effect sizes
head_regen %>%
  group_by(Condition) %>%
  summarize(
    d = cohen.d(PropActiveArm[Phase == "Regeneration"], PropActiveArm[Phase == "Baseline"])$estimate,
    .groups = "drop"
  )

tail_regen %>%
  group_by(Condition) %>%
  summarize(
    d = cohen.d(PropActiveArm[Phase == "Regeneration"], PropActiveArm[Phase == "Baseline"])$estimate,
    .groups = "drop"
  )

#########################
# ANALYSIS 3: Are there differences between head and tail memory?
# Compare head vs tail regenerates directly
#########################

head_tail_regen <- data %>% 
  filter(Phase == "Regeneration", BodyPart %in% c("Head", "Tail"))

# Statistical test
model_head_tail <- lmer(PropActiveArm ~ BodyPart * Condition + (1|OriginalSubject), data = head_tail_regen)
print(summary(model_head_tail))
print(anova(model_head_tail))

# Post-hoc comparisons
emm_head_tail <- emmeans(model_head_tail, ~ BodyPart | Condition)
print(pairs(emm_head_tail))

#########################
# ANALYSIS 5: If reinstatement data is available, analyze it
#########################

# Check if reinstatement data is available
has_reinstatement <- any(data$Phase == "Reinstatement", na.rm = TRUE)

if (has_reinstatement) {
  # Reinstatement data
  reinst_data <- data %>%
    filter(Phase %in% c("Baseline", "Reinstatement"))
  
  # Descriptive statistics
  reinst_summary <- reinst_data %>%
    group_by(Condition, BodyPart, Phase) %>%
    summarise(
      mean_prop = mean(PropActiveArm, na.rm = TRUE),
      sd_prop = sd(PropActiveArm, na.rm = TRUE),
      n = n(),
      se_prop = sd_prop / sqrt(n),
      .groups = "drop"
    )
  
  print("Summary statistics for reinstatement test:")
  print(reinst_summary)
  
  # Create Figure D: Reinstatement comparison plot
  # Filter for baseline (whole animals) and reinstatement data (head and tail)
  reinst_comparison <- data %>%
    filter((BodyPart == "Whole" & Phase == "Baseline") | 
             (BodyPart %in% c("Head", "Tail") & Phase == "Reinstatement"))
  
  # Create a new variable for the x-axis
  reinst_comparison <- reinst_comparison %>%
    mutate(BodyStatus = case_when(
      BodyPart == "Whole" ~ "Original",
      BodyPart == "Head" ~ "Head",
      BodyPart == "Tail" ~ "Tail"
    ))
  
  # Calculate group means and standard errors
  reinst_summary <- reinst_comparison %>%
    group_by(Condition, BodyStatus) %>%
    summarise(
      mean_prop = mean(PropActiveArm, na.rm = TRUE),
      se_prop = sd(PropActiveArm, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    )
  
  # Make sure conditions are in the desired order (Control, Cocaine, Methamphetamine)
  reinst_summary$Condition <- factor(reinst_summary$Condition, 
                                     levels = c("Control", "Cocaine", "Methamphetamine"))
  
  # Make sure body status is in the desired order (Original, Head, Tail)
  reinst_summary$BodyStatus <- factor(reinst_summary$BodyStatus,
                                      levels = c("Original", "Head", "Tail"))
  
  # Create point plot
  figD <- ggplot(reinst_summary, aes(x = BodyStatus, y = mean_prop, color = Condition, shape = BodyStatus)) +
    facet_wrap(~ Condition, scales = "free_x") +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop), width = 0.2) +
    scale_color_manual(values = c("Control" = "#159090", "Cocaine" = "#FF8C00", "Methamphetamine" = "#D81B60"),
                       name = "Condition") +
    scale_shape_manual(values = c("Original" = 16, "Head" = 17, "Tail" = 25),
                       name = "Body Part") +
    labs(title = "Reinstatement Active Arm Preference Compared to Baseline",
         y = "Proportion of Time on Active Surface",
         x = "Body Status") +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = 10),
      legend.position = "right",
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 12, hjust = 0.5),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    ) +
    ylim(0, 1)  # Match y-axis scale in your figure
  
  print(figD)
  
  # Statistical tests for head and tail
  head_reinst <- reinst_data %>% filter(BodyPart == "Head")
  tail_reinst <- reinst_data %>% filter(BodyPart == "Tail")
  
  # Only run models if there's sufficient data
  if (nrow(head_reinst) > 10) {
    model_head_reinst <- lmer(PropActiveArm ~ Phase * Condition + (1|OriginalSubject), data = head_reinst)
    print(summary(model_head_reinst))
    print(anova(model_head_reinst))
    
    emm_head_reinst <- emmeans(model_head_reinst, ~ Phase | Condition)
    print(pairs(emm_head_reinst))
  }
  
  if (nrow(tail_reinst) > 10) {
    model_tail_reinst <- lmer(PropActiveArm ~ Phase * Condition + (1|OriginalSubject), data = tail_reinst)
    print(summary(model_tail_reinst))
    print(anova(model_tail_reinst))
    
    emm_tail_reinst <- emmeans(model_tail_reinst, ~ Phase | Condition)
    print(pairs(emm_tail_reinst))
  }
  
  # Create a Combined Panel Figure with reinstatement
  combined_figure <- (figA | figA) / (figC | figD) + 
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Print combined figure
  print(combined_figure)
  
  # Save the combined figure
  ggsave("planaria_cpp_figures.png", combined_figure, width = 10, height = 8, dpi = 300)
  
} else {
  print("No reinstatement data available for analysis")
  
  # Create a Combined Panel Figure without reinstatement
  combined_figure <- (figA | figA) / (figC | plot_spacer()) + 
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Print combined figure
  print(combined_figure)
  
  # Save the combined figure
  ggsave("planaria_cpp_figures.png", combined_figure, width = 10, height = 8, dpi = 300)
}

# Create a summary table of excluded subjects
if(length(to_exclude) > 0) {
  exclusion_info <- data_original %>%
    filter(Subject %in% to_exclude) %>%
    select(Subject, Condition, body_part, Regeneration_test_distance) %>%
    left_join(group_medians, by = c("Condition", "body_part")) %>%
    mutate(percent_of_mean = (Regeneration_test_distance / median_distance) * 100)

}