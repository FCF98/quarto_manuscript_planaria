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



# plotting residuals

simulationOutput <- simulateResiduals(fittedModel = m1)
plot(simulationOutput)



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

####### Create inidividual plots for each subject with condition-based colors ##########


# First, create an ordered factor for subjects based on their condition
data_long_factor <- data_long %>%
  mutate(Subject_ordered = factor(Subject, 
                                  levels = data_long %>%
                                    select(Subject, Condition) %>%
                                    distinct() %>%
                                    arrange(desc(Condition == "Control"), Subject) %>%
                                    pull(Subject)))

# Then modify your plot to use the new ordered factor
individual_plots <- ggplot(data_long_factor, aes(x = Time, y = ActiveArmProportion, group = 1, color = Condition)) +
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
       width = 20, height = 24, dpi = 600, units = "in",
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




# Create the line plot with significance markers


summary_data <- data_long %>%
  group_by(Condition, Time) %>%
  summarise(
    mean_proportion = mean(ActiveArmProportion, na.rm = TRUE),
    se_proportion = sd(ActiveArmProportion, na.rm = TRUE) / sqrt(sum(!is.na(ActiveArmProportion))),
    n = sum(!is.na(ActiveArmProportion)),  # Add sample size calculation
    .groups = 'drop'
  )

ggplot(summary_data, aes(x = Time, y = mean_proportion, color = Condition, group = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_proportion - se_proportion, 
                    ymax = mean_proportion + se_proportion), 
                width = 0.2) +
  # Add significance markers for Control group
  geom_signif(
    data = data.frame(x = c(1, 2), xend = c(2, 3)),
    aes(x = x, xend = xend, y = c(0.65, 0.60), yend = c(0.65, 0.60)),
    annotation = c("***", "***"),
    color = "navy",
    map_signif_level = FALSE,
    tip_length = 0.02,
    vjust = 0.5
  ) +
  # Add significance markers for Treatment group
  geom_signif(
    data = data.frame(x = c(1, 1, 2), xend = c(2, 3, 4)),
    aes(x = x, xend = xend, y = c(0.55, 0.50, 0.45), yend = c(0.55, 0.50, 0.45)),
    annotation = c("**", "*", "*"),
    color = "darkred",
    map_signif_level = FALSE,
    tip_length = 0.02,
    vjust = 0.5
  ) +
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
    plot.margin = margin(b = 40, l = 20, r = 20, t = 20, unit = "pt")
  ) +
  scale_y_continuous(limits = c(-0.05, 0.7), breaks = seq(0, 0.7, 0.1)) +
  scale_color_brewer(palette = "Set1")



###### caclulating effect sizes ########

# Function to calculate Cohen's h

cohens_h <- function(p1, p2) {
  2 * (asin(sqrt(p1)) - asin(sqrt(p2)))
}

# 1. Cohen's h for treatment vs control at each time point

cohens_h_time <- data_long %>%
  group_by(Time, Condition) %>%
  summarize(mean_proportion = mean(ActiveArmProportion), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = mean_proportion) %>%
  mutate(cohens_h = cohens_h(Treatment, Control))

print("Cohen's h for Treatment vs Control at each time point:")
print(cohens_h_time)


# 2. Cohen's h for baseline vs endpoint within each condition
cohens_h_condition <- data_long %>%
  filter(Time %in% c("Baseline", "Endpoint")) %>%
  group_by(Condition, Time) %>%
  summarize(mean_proportion = mean(ActiveArmProportion), .groups = "drop") %>%
  pivot_wider(names_from = Time, values_from = mean_proportion) %>%
  mutate(cohens_h = cohens_h(Baseline, Endpoint))

print("Cohen's h for Baseline vs Endpoint within each condition:")
print(cohens_h_condition)

# 3. Odds ratios from GLMM models
# For m1 (Baseline vs Endpoint)
or_m1 <- exp(fixef(m1))
ci_m1 <- exp(confint(m1, method="Wald"))


# For m2 (Baseline vs Test)
or_m2 <- exp(fixef(m2))
ci_m2 <- exp(confint(m2, method="Wald"))

# For m3 (Baseline vs Reinstatement)
or_m3 <- exp(fixef(m3))
ci_m3 <- exp(confint(m3, method="Wald"))


# Function to pretty print odds ratios and CIs
print_or <- function(or, ci, model_name) {
  cat(paste0("\nOdds Ratios for ", model_name, ":\n"))
  for (i in 1:length(or)) {
    cat(sprintf("%s: OR = %.2f, 95%% CI [%.2f, %.2f]\n", 
                names(or)[i], or[i], ci[i,1], ci[i,2]))
  }
}



print_or(or_m1, ci_m1, "Model 1 (Baseline vs Endpoint)")
print_or(or_m2, ci_m2, "Model 2 (Baseline vs Test)")
print_or(or_m3, ci_m3, "Model 3 (Baseline vs Reinstatement)")



#these below are a seperate way to do the adjustment 

res <- emmeans(m1,pairwise~Time|Condition,adjust="none",type="response")

res <- as.data.frame(res$contrasts)

res$p.adj <- p.adjust(res$p.value,method="bonferroni")