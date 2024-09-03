library(ggplot2)
library(RColorBrewer)

dose_response_data <- read.csv('Datasets/DoseResponseStatistics.csv')

# I need my conditions as a factor

dose_response_data <- dose_response_data %>%
  mutate(Condition = as.factor(Condition))

# Violin plot

dose_response_violin <- ggplot(dose_response_data, aes(x = factor(Condition, levels = c("Control", "1uM", "5uM", "10uM", "20uM")), y = Distance, fill = Condition)) +
  geom_violin(trim = FALSE) +  # trim = FALSE to show full range of data
  geom_jitter(width = 0.2, alpha = 0.5) +  # Add data points with some jitter
  stat_summary(fun = mean, geom = "crossbar",color = "black", width = 0.5) +  # Add mean bars
  coord_cartesian(ylim = c(0, NA)) +  # Set y-axis limits with lower bound at 0
  theme_minimal() +  
  labs(title = "Distribution of Distance by Condition",
       x = "Condition",
       y = "Distance (cm)") +
  scale_fill_brewer(palette = "Pastel1")