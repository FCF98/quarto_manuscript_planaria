# library(ggplot2)
# library(RColorBrewer)
# library(dplyr)
# library(ggpubr)
# library(tidyr)
# library(ggridges)
# library(viridis)
# library(ggdist)
# library(ghibli)
# library(agridat)
# library(psych)

dose_response_data <- read.csv('Datasets/DoseResponseStatistics.csv')

time_binned_dose_response_data <- read.csv('Datasets/dose_response_time_binned.csv')




######### ridgeplot of average distance across time for whole sample ############




# Define the order of 'Minute' levels with highest values at the top

minute_order <- as.character(15:1)  # Adjust this to reverse the order, e.g., 15, 14, ..., 1

# Reorder the 'Minute' factor according to the specified order
time_binned_dose_response_data$Minute <- factor(time_binned_dose_response_data$Minute, levels = minute_order)

time_binned_dose_response_data$Distance <- as.numeric(time_binned_dose_response_data$Distance)

# Compute the means for each 'Minute'
time_binned_mean_data <- time_binned_dose_response_data %>%
  group_by(Minute) %>%
  summarise(mean_distance = mean(Distance, na.rm = TRUE))

# Create the ridgeline plot with mean lines
ggplot(time_binned_dose_response_data, aes(x = Distance, y = Minute, fill = as.factor(Minute))) +
  geom_density_ridges(alpha = 0.7, scale = 1.5) +  # Ridgeline part with transparency
  geom_point(aes(color = as.factor(Minute)), size = 3, 
             position = position_jitter(width = 0, height = 0.3), alpha = 0.6) +  # Spread data points vertically within the ridge
  geom_segment(data = time_binned_mean_data, aes(x = mean_distance, xend = mean_distance, 
                                                 y = as.numeric(Minute) - 0.3, yend = as.numeric(Minute) + 0.3), 
               color = "black", size = 1) +  # Add segment lines for mean in black
  scale_fill_viridis_d(guide = guide_legend(reverse = TRUE)) +  # Reverse legend order
  scale_color_viridis_d(guide = guide_legend(reverse = TRUE)) +  # Reverse legend order
  labs(x = "Distance Moved (cm)",
       y = "Time interval (minutes)",
       fill = "Minute",
       color = "Minute") +
  theme_minimal() +
  theme(legend.position = "right",  # Position legend on the right
        axis.title.y = element_text(margin = margin(r = 10), size = 16),  # Add space to the right of the y-axis label
        axis.title.x = element_text(margin = margin(t = 10), size = 16))  # Add space above the x-axis label


# t-test to compare minute 1 to minute 15 

minute_comparison_data <- time_binned_dose_response_data %>%
  filter(Minute %in% c("1", "15"))

minute_1 <- minute_comparison_data$Distance[minute_comparison_data$Minute == 1]
minute_15 <- minute_comparison_data$Distance[minute_comparison_data$Minute == 15]

#replacing NAs with 0
minute_15 <- replace_na(minute_15, 0)
minute_1 <- replace_na(minute_1, 0)

time_binned_t_test_result <-  t.test(minute_1, minute_15, na.rm = TRUE, paired = TRUE)

# Print the t-test result
print(time_binned_t_test_result)

# get the mean of min 1 and 15

minute1_descriptives <- as.data.frame(describe(minute_1))

minute15_descriptives <- as.data.frame(describe(minute_15), na.rm = TRUE)



############# Creating rain plot for distance ################



dose_rainplot_data <- read.csv('Datasets/dose_response_50subjects_trial.csv')

dose_rainplot_data <- mutate(dose_rainplot_data,
                             Subject = as.numeric(Subject),
                             Distance = as.numeric(Distance))

dose_rainplot_data$Condition <- factor(dose_rainplot_data$Condition,
                                       levels = c('Control', '1uM', '5uM', '10uM', '20uM', '100uM'))



# Create the plot
rainplot <- ggplot(dose_rainplot_data, aes(x = Condition, y = Distance, fill = Condition)) +
  scale_fill_ghibli_d("SpiritedMedium", direction = -1) +
  geom_boxplot(width = 0.1) +
  xlab('Condition') +
  ylab('Distance moved (cm)') +
  theme_classic(base_size = 18, base_family = "serif") +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 0, hjust = .5, vjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(r = 10)),  # Add space to the right of the y-axis label
        axis.title.x = element_text(margin = margin(t = 10))) +  # Add space above the x-axis label
  scale_y_continuous(breaks = seq(0, 180, by = 20), limits = c(0, 180), expand = c(0, 0)) +
  # Line below adds dot plots from {ggdist} package
  stat_dots(side = "left", justification = 1.12, binwidth = 1.9)  
# Line below adds half-violin from {ggdist} package 
 # stat_halfeye(adjust = .5, width = .6, justification = -.2, .width = 0, point_colour = NA)

print(rainplot)




############# Creating boxplot for distance ################

data_summary <- dose_response_data %>%
  group_by(Condition) %>%
  summarise(
    mean = mean(Distance),
    sem = sd(Distance) / sqrt(n())
  )

dose_response_boxplot <- ggplot(dose_response_data, aes(x = Condition, y = Distance, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +  # Add jitter points
  scale_fill_brewer(palette = "Pastel1") +  # Use Pastel1 color scheme
  theme_minimal() +
  labs(title = "Planaria Cocaine Dose-Response",
       x = "Condition",
       y = "Distance moved (cm)") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # Remove legend
  )

dose_response_boxplot



########## creating line grpagh for time binned data ##############

# Load required libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Reorder Condition factor levels for the legend
time_binned_dose_response_data$Condition <- factor(time_binned_dose_response_data$Condition,
                                                   levels = c("control", "1uM", "5uM", "10uM", "20uM", "100uM"))

# Create summary statistics for plotting (Mean and SEM for each Condition and Minute)
plot_data <- time_binned_dose_response_data %>%
  group_by(Condition, Minute) %>%
  summarize(MeanDistance = mean(Distance, na.rm = TRUE),
            SEMDistance = sd(Distance, na.rm = TRUE) / sqrt(n()))  # Calculate SEM

# Use a bolder color palette (modified manually for better visibility)
strong_colors <- c("control" = "#377eb8",  # Blue
                   "1uM" = "#e41a1c",     # Red
                   "5uM" = "#4daf4a",     # Green
                   "10uM" = "#ff7f00",    # Orange
                   "20uM" = "#984ea3", # Purple
                  "100uM" = "#A3984E") #gold

# Create the plot with SEM error bars and custom strong colors

ggplot(plot_data, aes(x = Minute, y = MeanDistance, color = Condition)) +
  geom_line(size = 1.2) +  # Line plot with thicker lines
  geom_point(size = 2.5) +  # Larger points
  geom_errorbar(aes(ymin = MeanDistance - SEMDistance, ymax = MeanDistance + SEMDistance), 
                width = 0.2) +  # Error bars with SEM
  labs(x = "Time bins (minutes)",
       y = "Total Distance (cm)") +  # No title
  scale_color_manual(values = strong_colors) +  # Apply custom strong colors
  theme_minimal() +  # Minimal theme
  theme(
    legend.title = element_blank(),  # Remove legend title
    legend.position = "right",  # Keep legend on the right
    axis.title.x = element_text(margin = margin(t = 15)),  # Add space between x-axis label and plot
    axis.title.y = element_text(margin = margin(r = 15)),  # Add space between y-axis label and plot
    panel.grid.major = element_line(color = "gray", linetype = "dashed"),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_blank()  # Remove the plot title
  )

########## Boxplot of distance moved in first  5 mins ###########


# Ensure 'Minute' is numeric for filtering
time_binned_dose_response_data$Minute <- as.numeric(as.character(time_binned_dose_response_data$Minute))

# Filter data for minutes 1 through 5
filtered_data <- time_binned_dose_response_data %>%
  filter(Minute >= 1 & Minute <= 5)

# Summarize total distance for each condition and subject
total_distance_data <- filtered_data %>%
  group_by(Condition, Subject) %>%
  summarise(TotalDistance = sum(Distance, na.rm = TRUE), .groups = 'drop')

# Create the boxplot
boxplot_graph <- ggplot(total_distance_data, aes(x = Condition, y = TotalDistance, fill = Condition)) +
  geom_boxplot() +
  labs(title = "Total Distance by Condition (Minutes 1-5)",
       x = "Condition",
       y = "Total Distance") +
  theme_minimal()

# Display the boxplot
print(boxplot_graph)




