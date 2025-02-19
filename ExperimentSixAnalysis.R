library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(effsize)

# Read in the data
data <- read_excel("Exp 6 CPP v2 data for test analysis.xlsx")

# Fix the column name with the space
names(data) <- gsub(" ", "_", names(data))

# Clean the initial data
clean_data <- data %>%
  filter(!is.na(Subject)) %>%
  select(Subject, Condition, Baseline_active_proportion, Test_active_proportion)

# Create long format data
long_data <- clean_data %>%
  pivot_longer(
    cols = c("Baseline_active_proportion", "Test_active_proportion"),
    names_to = "Time",
    values_to = "Active_Preference"
  ) %>%
  mutate(
    Time = factor(ifelse(Time == "Baseline_active_proportion", "Baseline", "Test")),
    Condition = factor(Condition),
    Subject = factor(Subject)  # Convert Subject to factor
  )

# Fit mixed effects model
model <- lmer(Active_Preference ~ Condition * Time + (1|Subject), 
              data = long_data)

# Create visualization
plot <- ggplot(long_data, aes(x = Time, y = Active_Preference, 
                              color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_point(alpha = 0.3, position = position_jitter(width = 0.1)) +
  theme_classic() +
  labs(
    title = "Change in Active Surface Preference: Baseline vs Test",
    y = "Proportion of Time in Active Surface",
    x = "Time Point"
  ) +
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)

# Print plot
print(plot)

# Print model results using anova
print("ANOVA results from mixed model:")
print(anova(model))

# Calculate descriptive statistics
aggregate_stats <- long_data %>%
  group_by(Condition, Time) %>%
  summarise(
    mean = mean(Active_Preference, na.rm = TRUE),
    se = sd(Active_Preference, na.rm = TRUE)/sqrt(n()),
    n = n(),
    .groups = 'drop'
  )

print("Descriptive Statistics:")
print(aggregate_stats)

# Add post-hoc comparisons
# Within-group comparisons (Baseline vs Test for each condition)
pairwise_tests <- long_data %>%
  group_by(Condition) %>%
  summarise(
    t_test = list(t.test(Active_Preference ~ Time, paired = TRUE)),
    p_value = t_test[[1]]$p.value,
    t_stat = t_test[[1]]$statistic,
    .groups = 'drop'
  )

print("\nPairwise comparisons (Baseline vs Test) within each condition:")
print(pairwise_tests %>% select(Condition, p_value, t_stat))

# Between-group comparisons at each timepoint
between_group_tests <- long_data %>%
  group_by(Time) %>%
  summarise(
    aov_result = list(aov(Active_Preference ~ Condition)),
    p_value = summary(aov_result[[1]])[[1]][["Pr(>F)"]][1],
    f_stat = summary(aov_result[[1]])[[1]][["F value"]][1],
    .groups = 'drop'
  )

print("\nBetween-group comparisons at each timepoint:")
print(between_group_tests %>% select(Time, p_value, f_stat))


# Add effect size calculations
cohen_d <- long_data %>%
  group_by(Condition) %>%
  summarise(
    effect_size = cohen.d(Active_Preference ~ Time, paired = TRUE)$estimate
  )

# More detailed descriptive statistics
detailed_stats <- long_data %>%
  group_by(Condition, Time) %>%
  summarise(
    mean = mean(Active_Preference),
    sd = sd(Active_Preference),
    n = n(),
    sem = sd/sqrt(n)
  )

print("Effect sizes:")
print(cohen_d)
print("\nDetailed statistics:")
print(detailed_stats)