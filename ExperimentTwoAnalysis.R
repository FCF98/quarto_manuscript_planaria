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

exp2_model_output <- car::Anova(m1, type = "III")

emmeans(m1,pairwise~Time|Condition,adjust="bonferroni",type="response")

emmeans(m1,pairwise~Condition|Time,adjust="bonferroni",type="response")

#nice_output <- nice(m1, sig_symbols = c(" = ", "<"))


## trying afex format ##

#m1_afex <- mixed(cbind(ActiveCount, InactiveCount) ~ Condition * Time + (1|Subject), 
               #  data = data_long, 
                # family = binomial,
                # method = "LRT")

#################################


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



# Create APA-style plot with both types of comparisons
grouped_comprison_w_significance_apa<- ggplot(summary_data_w_sample, aes(x = Time, y = mean_proportion, color = Condition, group = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_proportion - se_proportion, 
                    ymax = mean_proportion + se_proportion), 
                width = 0.2) +
  # Control group within comparisons
  geom_signif(
    annotations = "***",
    xmin = 1, xmax = 2,
    y_position = 0.63,
    color = "#159090"
  ) +
  geom_signif(
    annotations = "***",
    xmin = 2, xmax = 3,
    y_position = 0.8,
    color = "#159090"
  ) +
  # Treatment group within comparisons
  geom_signif(
    annotations = "**",
    xmin = 1, xmax = 2,
    y_position = 0.7,
    color = "#FF8C00"
  ) +
  geom_signif(
    annotations = "*",
    xmin = 1, xmax = 3,
    y_position = 0.88,
    color = "#FF8C00"
  ) +
  geom_signif(
    annotations = "*",
    xmin = 2, xmax = 4,
    y_position = 0.95,
    color = "#FF8C00"
  ) +
  # Between group comparison
  geom_signif(
    annotations = "#",
    xmin = 2.9, xmax = 3.1,
    y_position = 0.6,
    color = "black"
  ) +
  # Add sample size labels below the x-axis
  geom_text(aes(label = label, y = 0.04), 
            position = position_dodge(width = 0.4),
            size = 3,
            show.legend = FALSE) +
  # APA-style labels
  labs(
    x = "Time Point",
    y = "Mean Proportion of Active Arm Choices",
    color = "Condition"
  ) +
  # APA-style theme modifications
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 12),
    axis.title = element_text(size = 12),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(t = 20, r = 20, b = 40, l = 20, unit = "pt"),
    panel.border = element_blank()
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    expand = c(0.02, 0)
  ) +
  scale_color_manual(values = c("Control" = "#159090", "Treatment" = "#FF8C00"))

grouped_comprison_w_significance_apa

#saving the plot
ggsave("grouped_comprison_w_significance_apa.png", grouped_comprison_w_significance_apa, 
       width = 12, height = 9, dpi = 300, units = "in",
       bg = "white")



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

#Ctr + Shft + c to uncomment

# Data_long_days <- data %>%
#   select(Subject, Condition,
#          Baseline = `Baseline%`,
#          Day1 = `Day1%`,
#          Day2 = `Day2%`,
#          Day3 = `Day3%`,
#          Day4 = `Day4%`,
#          Day5 = `Day5%`) %>%
#   pivot_longer(
#     cols = c(Baseline:Day5),
#     names_to = "TimePoint",
#     values_to = "Proportion"
#   ) %>%
#   mutate(
#     TimePoint = factor(TimePoint,
#                        levels = c("Baseline", "Day1", "Day2", "Day3", "Day4", "Day5")),
#     Condition = factor(Condition)
#   ) %>%
#   # Calculate mean and SE for each condition and timepoint
#   group_by(Condition, TimePoint) %>%
#   summarise(
#     mean_prop = mean(Proportion, na.rm = TRUE),
#     se = sd(Proportion, na.rm = TRUE) / sqrt(n()),
#     .groups = 'drop'
#   )
# 
# # Create the visualization
# ggplot(Data_long_days, aes(x = TimePoint, y = mean_prop,
#                       color = Condition, group = Condition)) +
#   # Add mean lines
#   geom_line(linewidth = 1.5) +
#   # Add mean points
#   geom_point(size = 3) +
#   # Add error bars
#   geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
#                 width = 0.2) +
#   # Customize the theme and labels
#   theme_minimal() +
#   labs(
#     title = "Average Proportion of Active Arm Choices Over Time",
#     y = "Proportion of Active Arm Choices",
#     x = "Time Point"
#   ) +
#   scale_color_manual(values = c("Control" = "#0072B2", "Treatment" = "#D55E00")) +
#   theme(
#     legend.position = "bottom",
#     panel.grid.minor = element_blank(),
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 12, face = "bold"),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   ) +
#   ylim(0, 1)  # Set y-axis limits from 0 to 1




