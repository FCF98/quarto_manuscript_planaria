library(readxl)
library(tidyr)
library(dplyr)
library(lme4)

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
data_subset <- subset(data_long, Time %in% c("Baseline", "Endpoint"))

# # model comparing baseline to endpoint using subsetted data
m1 <- glmer(cbind(ActiveCount, InactiveCount) ~ Condition * Time + (1|Subject), 
            data = data_subset, 
            family = "binomial")

summary(m1)

car::Anova(m1, type = "III")