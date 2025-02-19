library(ggplot2)
library(dplyr)

# Create data for single cohort with paired training
start_date <- as.Date("2024-10-29")  # arbitrary start date

tasks <- data.frame(
  Task = c("Baseline (all subjects - 40 x 15min)", 
           "Training Days (20 pairs x 20min/day)",
           "First Test (40 subjects x 15min)",
           "Processing",
           "Second Test (40 subjects x 15min)",
           "Final Test (40 subjects x 20min)"),
  Start = c(
    start_date,
    start_date + 1,
    start_date + 11,
    start_date + 12,
    start_date + 26,
    start_date + 27
  ),
  Duration = c(1, 10, 1, 1, 1, 1)
)

# Calculate end dates
tasks$End <- tasks$Start + tasks$Duration

# Create the Gantt chart
ggplot(tasks, aes(y = reorder(Task, desc(Start)), x = Start, xend = End)) +
  geom_segment(size = 10, aes(yend = Task), color = "skyblue") +
  labs(title = "Experimental Schedule - Single Cohort with Paired Training",
       subtitle = "Running control + treatment pairs during training days",
       x = "Date",
       y = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.text = element_text(size = 10)) +
  scale_x_date(date_breaks = "2 days", date_labels = "%b %d")


















