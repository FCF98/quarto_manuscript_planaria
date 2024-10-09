library(readxl)
library(tidyr)
library(dplyr)
library(lme4)
library(ggplot2)
library(hrbrthemes)
library(effectsize)

data <- read_excel("Datasets/PipelineAnalysisTestData2.xlsx")


#########structuring data correctly############

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


###### For model 2 and 3, need to decide whether to test cocaine group (Treatment) only

# Subset the data for Baseline and test
data_subset_BT <- subset(data_long, Time %in% c("Baseline", "Test"))

m2 <- glmer(cbind(ActiveCount, InactiveCount) ~ Condition * Time + (1|Subject), 
            data = data_subset_BT, 
            family = "binomial")

summary(m2)

# Model 2: Baseline vs Test (only for treatment group)
data_subset_BT_treatment <- subset(data_subset_BT, Condition == "Treatment")
m2 <- glmer(cbind(ActiveCount, InactiveCount) ~ Time + (1|Subject), 
            data = data_subset_BT_treatment, 
            family = "binomial")

# Subset the data for Baseline and reinstatement
data_subset_BR <- subset(data_long, Time %in% c("Baseline", "Reinstatement"))

m3 <- glmer(cbind(ActiveCount, InactiveCount) ~ Condition * Time + (1|Subject), 
            data = data_subset_BR, 
            family = "binomial")

summary(m3)

# Model 3: Baseline vs Reinstatement (only for treatment group)
data_subset_BR_treatment <- subset(data_subset_BR, Condition == "Treatment")
m3 <- glmer(cbind(ActiveCount, InactiveCount) ~ Time + (1|Subject), 
            data = data_subset_BR_treatment, 
            family = "binomial")




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




########### Correcting for multiple comparisons #########

# Function to extract p-values from model summary
extract_p_values <- function(model) {
  summary_data <- summary(model)$coefficients
  p_values <- summary_data[, "Pr(>|z|)"]
  return(p_values[-1])  # Exclude intercept
}

# Extract p-values from all models
p_values_m1 <- extract_p_values(m1)
p_values_m2 <- extract_p_values(m2)
p_values_m3 <- extract_p_values(m3)

# Combine all p-values
all_p_values <- c(p_values_m1, p_values_m2, p_values_m3)

# Apply Bonferroni correction
n_tests <- length(all_p_values)
adjusted_p_values <- p.adjust(all_p_values, method = "holm")

# Function to print results with adjusted p-values
print_results <- function(model, adjusted_p, model_name) {
  cat(paste0("\nResults for ", model_name, ":\n"))
  summary_data <- summary(model)$coefficients
  coefficients <- summary_data[, "Estimate"]
  std_errors <- summary_data[, "Std. Error"]
  
  for (i in 2:nrow(summary_data)) {  # Start from 2 to skip intercept
    cat(sprintf("%s: Estimate = %.4f, SE = %.4f, Adjusted p-value = %.4f\n",
                rownames(summary_data)[i], coefficients[i], std_errors[i], 
                adjusted_p[i-1]))
  }
}

# Print results for each model
print_results(m1, adjusted_p_values[1:3], "Model 1 (Baseline vs Endpoint)")
print_results(m2, adjusted_p_values[4:6], "Model 2 (Baseline vs Test)")
print_results(m3, adjusted_p_values[7:9], "Model 3 (Baseline vs Reinstatement)")



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


