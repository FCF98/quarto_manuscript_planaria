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
ridgeplot <- ggplot(time_binned_dose_response_data, aes(x = Distance, y = Minute, fill = as.factor(Minute))) +
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
