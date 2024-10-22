library(readxl)
library(tidyr)
library(dplyr)
library(lme4)
library(ggplot2)
library(hrbrthemes)
library(effectsize)
library(emmeans)
library(DHARMa)

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
  filter(is.na(ActiveCount) | ActiveCount <= TotalTrials)



###### model to test for main effects and interactions across 4 time points #####

contrasts(data_long$Condition) <- contr.sum

contrasts(data_long$Time) <- contr.sum

m1 <- glmer(cbind(ActiveCount, InactiveCount) ~ Condition * Time + (1|Subject), 
            data = data_long, 
            family = "binomial",
            na.action = na.omit)

summary(m1)

car::Anova(m1, type = "III")

emmeans(m1,pairwise~Time|Condition,adjust="bonferroni",type="response")

emmeans(m1,pairwise~Condition|Time,adjust="bonferroni",type="response")


############# plotting the grouped results   ############


# Calculate mean and standard error for each Condition and Time, explicitly handling NAs
summary_data <- data_long %>%
  group_by(Condition, Time) %>%
  summarise(
    mean_proportion = mean(ActiveArmProportion, na.rm = TRUE),
    se_proportion = sd(ActiveArmProportion, na.rm = TRUE) / sqrt(sum(!is.na(ActiveArmProportion))),
    n = sum(!is.na(ActiveArmProportion)),  # Add sample size calculation
    .groups = 'drop'
  )

# Add sample size to the plot labels
summary_data <- summary_data %>%
  mutate(label = sprintf("n=%d", n))

# Create the line plot with sample sizes
ggplot(summary_data, aes(x = Time, y = mean_proportion, color = Condition, group = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_proportion - se_proportion, 
                    ymax = mean_proportion + se_proportion), 
                width = 0.2) +
  # Add sample size labels below the x-axis
  geom_text(aes(label = label, y = -0.02), 
            position = position_dodge(width = 0.4),
            size = 3,
            show.legend = FALSE) +
  labs(
    title = "Active Arm Choice Proportion Over Time",
    subtitle = "By Experimental Condition",
    x = "Time Point",
    y = "Proportion of Active Arm Choices",
    color = "Condition"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # Add more space at the bottom for sample size labels
    plot.margin = margin(b = 40, l = 20, r = 20, t = 20, unit = "pt")
  ) +
  scale_y_continuous(limits = c(-0.05, 0.7), breaks = seq(0, 0.7, 0.1)) +
  scale_color_brewer(palette = "Set1")

# Print the summary data to check the calculations
print(summary_data)



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
  scale_color_manual(values = c("Control" = "navy", "Treatment" = "darkred"))

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


