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
  

#Basic plot

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
                position = position_dodge(0.8), width = 0.2, linewidth = 0.7) +
  scale_fill_manual(values = c("Treatment" = "#FF8C00", "Control" = "#159090")) +
  labs( 
    title = "Active Arm Choices Before and After Conditioning",
    y = "Proportion of Active Arm Entries",
    x = "Time point",
    fill = "Condition") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = (seq(0, 1, 0.1)
  )) +
  consistant_theme()


print(Exp8_Baseline_endpoint_comparison)
