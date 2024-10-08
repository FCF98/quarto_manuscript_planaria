library(readxl)
library(tidyr)
library(dplyr)
library(lme4)
library(ggplot2)
library(hrbrthemes)

data <- read_excel("Datasets/PipelineAnalysisTestData2.xlsx")


#structuring data correctly

data_long <- data %>%
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
    ActiveCount = as.integer(ActiveCount)  # Ensure this is an integer
  ) %>%
  mutate(
    TotalTrials = case_when(
      Time %in% c("Baseline", "Endpoint") ~ 6,
      Time %in% c("Test", "Reinstatement") ~ 4
    ),
    InactiveCount = TotalTrials - ActiveCount,
    ActiveArmProportion = ActiveCount / TotalTrials
  ) %>%
  filter(ActiveCount <= TotalTrials)  # Remove any impossible data points


# Subset the data for Baseline and Endpoint
data_subset_BE <- subset(data_long, Time %in% c("Baseline", "Endpoint"))

###### model comparing baseline to endpoint using subsetted data #####
m1 <- glmer(cbind(ActiveCount, InactiveCount) ~ Condition * Time + (1|Subject), 
            data = data_subset_BE, 
            family = "binomial")

summary(m1)

car::Anova(m1, type = "III")


# Subset the data for Baseline and test
data_subset_BT <- subset(data_long, Time %in% c("Baseline", "Test"))

m2 <- glmer(cbind(ActiveCount, InactiveCount) ~ Condition * Time + (1|Subject), 
            data = data_subset_BT, 
            family = "binomial")

summary(m2)

# Subset the data for Baseline and reinstatement
data_subset_BR <- subset(data_long, Time %in% c("Baseline", "Reinstatement"))

m3 <- glmer(cbind(ActiveCount, InactiveCount) ~ Condition * Time + (1|Subject), 
            data = data_subset_BR, 
            family = "binomial")

summary(m3)



############# plotting the results   ############


# Calculate mean and standard error for each Condition and Time
summary_data <- data_long %>%
  group_by(Condition, Time) %>%
  summarise(
    mean_proportion = mean(ActiveArmProportion),
    se_proportion = sd(ActiveArmProportion) / sqrt(n()),
    .groups = 'drop'
  )


# Create the line plot
ggplot(summary_data, aes(x = Time, y = mean_proportion, color = Condition, group = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_proportion - se_proportion, 
                    ymax = mean_proportion + se_proportion), 
                width = 0.2) +
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
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  scale_color_brewer(palette = "Set1")


# creating a interaction plot with by-participant data

ggplot(data_long, aes(x = Time, y = ActiveArmProportion, 
                     group = Condition, shape = Condition)) +
  # adds raw data points in each condition
  geom_point(aes(colour = Condition),alpha = .2) +
  # add lines to connect each participant's data points across conditions
  geom_line(aes(group = Subject, colour = Condition), alpha = .2) +
  # add data points representing cell means
  stat_summary(fun = "mean", geom = "point", size = 2, colour = "black") +
  # add lines connecting cell means by condition
  stat_summary(fun = "mean", geom = "line", colour = "black") +
  # add errorbars to cell means
  stat_summary(fun.data = "mean_se", geom = "errorbar", 
               width = .2, colour = "black") +
  # change colours and theme
  scale_color_brewer(palette = "Dark2") +
  theme_minimal()



###### lollipop chart ####

# Change baseline
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


