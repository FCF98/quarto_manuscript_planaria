#library(ggplot2)
#library(RColorBrewer)
#library(dplyr)
#library(ggpubr)
#library(tidyr)
#library(ggridges)
#library(viridis)
#library(ggdist)
#library(ghibli)
#library(agridat)

dose_rainplot_data <- read.csv('Datasets/dose_response_50subjects_trial.csv')

dose_rainplot_data <- mutate(dose_rainplot_data,
                             Subject = as.numeric(Subject),
                             Distance = as.numeric(Distance))

dose_rainplot_data$Condition <- factor(dose_rainplot_data$Condition,
                                       levels = c('Control', '1uM', '5uM', '10uM', '20uM'))



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
        axis.title.y = element_text(margin = margin(r = 20)),  # Add space to the right of the y-axis label
        axis.title.x = element_text(margin = margin(t = 20))) +  # Add space above the x-axis label
  scale_y_continuous(breaks = seq(0, 180, by = 20), limits = c(0, 180), expand = c(0, 0)) +
  # Line below adds dot plots from {ggdist} package
  stat_dots(side = "left", justification = 1.12, binwidth = 1.9)  
# Line below adds half-violin from {ggdist} package 
# stat_halfeye(adjust = .5, width = .6, justification = -.2, .width = 0, point_colour = NA)


