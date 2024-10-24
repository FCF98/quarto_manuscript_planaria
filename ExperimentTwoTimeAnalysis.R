library(readxl)
library(tidyverse)
library(lme4)
library(emmeans)
library(ggsignif)
library(ggplot2)
library(gghalves)

Exp2_Time_Data <- read_excel("Datasets/ExpeimentTwoTimeData.xlsx")


########### filtering out subjects that met exclusion criteria ###########

Exp2_Time_Data_Clean <- Exp2_Time_Data %>%
  filter(!Subject %in% c(3, 4, 14, 22, 27, 35, 40, 44, 47, 50, 52, 55, 56))


########## Converting time to seconds #############

# Vectorized function to convert time values
convert_time <- function(x) {
  # Convert to numeric, handling NAs automatically
  x_num <- as.numeric(x)
  
  # Create result vector
  result <- ifelse(is.na(x_num), NA,
                   ifelse(x_num > 1, x_num, x_num * 1440))
  
  return(result)
}

# Get all time column names
time_cols <- names(Exp2_Time_Data_Clean)[grep("Time", names(Exp2_Time_Data_Clean))]

# Apply conversion to all time columns
Exp2_Time_Data_converted <- Exp2_Time_Data_Clean %>%
  mutate(across(all_of(time_cols), convert_time))



############ Creating mean completion time per day #############

# Calculate means for different phases
Exp2_Time_Data_Clean <- Exp2_Time_Data_converted %>%
  rowwise() %>%
  mutate(
    # Baseline mean (B1-B6)
    Baseline_Mean = mean(c(B1_Time, B2_Time, B3_Time, B4_Time, B5_Time, B6_Time), na.rm = TRUE),
    
    # Conditioning means by day (4 trials per day)
    Cond_Day1_Mean = mean(c(Time_C1, Time_C2, Time_C3, Time_C4), na.rm = TRUE),
    Cond_Day2_Mean = mean(c(Time_C5, Time_C6, Time_C7, Time_C8), na.rm = TRUE),
    Cond_Day3_Mean = mean(c(Time_C9, Time_C10, Time_C11, Time_C12), na.rm = TRUE),
    Cond_Day4_Mean = mean(c(Time_C13, Time_C14, Time_C15, Time_C16), na.rm = TRUE),
    Cond_Day5_Mean = mean(c(Time_C17, Time_C18, Time_C19, Time_C20), na.rm = TRUE),
    
    # Cacluating endpoint value (last 6 conditioning trials)
    Endpoint = mean(c(Time_C15, Time_C16, Time_C17, Time_C18, Time_C19, Time_C20), na.rm = TRUE),
    
    # Test mean (T1-T4)
    Test_Mean = mean(c(Time_T1, Time_T2, Time_T3, Time_T4), na.rm = TRUE),
    
    # Reinstatement mean (R1-R4)
    Reinst_Mean = mean(c(Time_R1, Time_R2, Time_R3, Time_R4), na.rm = TRUE)
  ) %>%
  ungroup()



########### Converting data to long format ################

# Convert to long format using the existing mean columns
Exp2_Time_Long <- Exp2_Time_Data_Clean %>%
  select(Subject, Condition, 
         Baseline = Baseline_Mean,
         Day1 = Cond_Day1_Mean,
         Day2 = Cond_Day2_Mean,
         Day3 = Cond_Day3_Mean,
         Day4 = Cond_Day4_Mean,
         Day5 = Cond_Day5_Mean,
         Endpoint = Endpoint,
         Test = Test_Mean,
         Reinst = Reinst_Mean) %>%
  pivot_longer(
    cols = c(Baseline:Reinst),
    names_to = "TimePoint",
    values_to = "Time"
  ) %>%
  mutate(TimePoint = factor(TimePoint, 
                            levels = c("Baseline", "Day1", "Day2", "Day3", "Day4", "Day5", "Endpoint", "Test", "Reinst")))


######### Creating a model and testing change in completion time across time points between conditions ##########

# Need to remove time points for conditiong day averages

# Convert to long format using the existing mean columns
Exp2_Time_Long_Four_Points <- Exp2_Time_Data_Clean %>%
  select(Subject, Condition, 
         Baseline = Baseline_Mean,
         Endpoint = Endpoint,
         Test = Test_Mean,
         Reinst = Reinst_Mean) %>%
  pivot_longer(
    cols = c(Baseline:Reinst),
    names_to = "TimePoint",
    values_to = "Time"
  ) %>%
  mutate(TimePoint = factor(TimePoint, 
                            levels = c("Baseline", "Endpoint", "Test", "Reinst")))



m2 <- lmer(log(Time) ~ Condition * TimePoint + (1|Subject), 
           data = Exp2_Time_Long_Four_Points, 
           na.action = na.omit)

summary(m2)

car::Anova(m2, type = "III")

emmeans(m2,pairwise~TimePoint|Condition,adjust="bonferroni",type="response")

emmeans(m2,pairwise~Condition|TimePoint,adjust="bonferroni",type="response")






######### Plotting the grouped data by condition across days #############


grouped_time_plots <- ggplot(Exp2_Time_Long_Four_Points, aes(x = TimePoint, y = Time, color = Condition, group = Condition)) +
  # Add mean lines
  stat_summary(fun = mean, geom = "line", size = 1) +
  # Add error bars (standard error)
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  # Add individual subject points (optional - comment out if too cluttered)
  geom_point(alpha = 0.2) +
  # Customize appearance
  theme_classic() +
  labs(
    x = "Time Point",
    y = "Time to Decision (seconds)",
    title = "Decision Time Across Experimental Phases",
    color = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) 

grouped_time_plots

# save the plot
ggsave("grouped_time_plots.png", grouped_time_plots, 
       width = 8, height = 10, dpi = 300, units = "in",
       bg = "white")




########### Trying to plot with significance bars ########

grouped_time_plots_with_sig <- ggplot(Exp2_Time_Long_Four_Points, aes(x = TimePoint, y = Time, color = Condition, group = Condition)) +
  # Add individual subject points with jitter and dodge
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.4),
             alpha = 0.3,  # slightly increased transparency
             size = 2) +   # increased point size
  
  # Add mean lines - dodged to align with points
  stat_summary(fun = mean, 
               geom = "line", 
               linewidth = 1.2,    # slightly thicker lines
               position = position_dodge(0.4)) +
  
  # Add error bars (standard error) - dodged to align with points
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               width = 0.2,
               linewidth = 1,      # thicker error bars
               position = position_dodge(0.4)) +
  
  # Control group significance markers
  geom_signif(
    annotations = "***",
    xmin = "Baseline", xmax = "Test",
    y_position = 200,
    color = "#377EB8"
  ) +
  geom_signif(
    annotations = "***",
    xmin = "Endpoint", xmax = "Test",
    y_position = 300,
    color = "#377EB8"
  ) +
  
  # Treatment group significance markers
  geom_signif(
    annotations = "***",
    xmin = "Baseline", xmax = "Test",
    y_position = 230,
    color = "#E60012"
  ) +
  geom_signif(
    annotations = "***",
    xmin = "Baseline", xmax = "Reinst",
    y_position = 260,    # fixed the y-position from 2600
    color = "#E60012"
  ) +
  geom_signif(
    annotations = "**",
    xmin = "Endpoint", xmax = "Test",
    y_position = 330,
    color = "#E60012"
  ) +
  geom_signif(
    annotations = "**",
    xmin = "Endpoint", xmax = "Reinst",
    y_position = 360,
    color = "#E60012"
  ) +
  
  # APA-style theme modifications
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 12),
    axis.title = element_text(size = 12),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(t = 20, r = 20, b = 40, l = 20, unit = "pt"),
    panel.border = element_blank()
  ) +
  scale_y_continuous(
    limits = c(0, 400),
    breaks = seq(0, 300, 50)
  ) +
  scale_color_manual(values = c("Control" = "#377EB8", "Treatment" = "#E60012")) +
  labs(
    x = "Time Point",
    y = "Time to Decision (seconds)",
    color = "Condition"
  )

grouped_time_plots_with_sig

ggplot(Exp2_Time_Long_Four_Points, aes(x = TimePoint, y = Time, color = Condition, fill = Condition)) +
  # Add boxplots with transparency
  geom_boxplot(position = position_dodge(0.8), 
               alpha = 0.2, 
               width = 0.5) +
  # Add individual points with jitter
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
             alpha = 0.3, 
             size = 2) +
  # Add mean lines
  stat_summary(fun = mean, 
               geom = "line", 
               linewidth = 1.2,
               position = position_dodge(0.8),
               aes(group = Condition)) +
  # Add your existing significance markers and theme
  # Control group significance markers
  geom_signif(
    annotations = "***",
    xmin = "Baseline", xmax = "Test",
    y_position = 200,
    color = "#377EB8"
  ) +
  geom_signif(
    annotations = "***",
    xmin = "Endpoint", xmax = "Test",
    y_position = 300,
    color = "#377EB8"
  ) +
  
  # Treatment group significance markers
  geom_signif(
    annotations = "***",
    xmin = "Baseline", xmax = "Test",
    y_position = 230,
    color = "#E60012"
  ) +
  geom_signif(
    annotations = "***",
    xmin = "Baseline", xmax = "Reinst",
    y_position = 260,    # fixed the y-position from 2600
    color = "#E60012"
  ) +
  geom_signif(
    annotations = "**",
    xmin = "Endpoint", xmax = "Test",
    y_position = 330,
    color = "#E60012"
  ) +
  geom_signif(
    annotations = "**",
    xmin = "Endpoint", xmax = "Reinst",
    y_position = 360,
    color = "#E60012"
  ) +
  
  # APA-style theme modifications
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 12),
    axis.title = element_text(size = 12),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(t = 20, r = 20, b = 40, l = 20, unit = "pt"),
    panel.border = element_blank()
  ) +
  scale_color_manual(values = c("Control" = "#377EB8", "Treatment" = "#E60012")) +
  scale_fill_manual(values = c("Control" = "#377EB8", "Treatment" = "#E60012"))




########## Raincloud plot ############


library(ggridges)

ggplot(Exp2_Time_Long_Four_Points, 
       aes(x = Time, y = TimePoint, fill = Condition)) +
  geom_density_ridges(alpha = 0.5, 
                      jittered_points = TRUE,
                      point_alpha = 0.3,
                      position = position_points_jitter(width = 0.05, height = 0)) +
  theme_ridges()

ggplot(Exp2_Time_Long_Four_Points, aes(x = TimePoint, y = Time, color = Condition)) +
  # Individual trajectories (light)
  geom_line(aes(group = interaction(Subject, Condition)), alpha = 0.1) +
  # Group means (bold)
  stat_summary(fun = mean, geom = "line", size = 2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  theme_minimal()

############ plotting individual data ################

# First, create an ordered factor for subjects based on their condition
Exp2_Time_Long <- Exp2_Time_Long %>%
  mutate(Subject_ordered = factor(Subject, 
                                  levels = Exp2_Time_Long %>%
                                    select(Subject, Condition) %>%
                                    distinct() %>%
                                    arrange(desc(Condition == "Control"), Subject) %>%
                                    pull(Subject)))

# Create individual plots
individual_time_plots <- ggplot(Exp2_Time_Long, aes(x = TimePoint, y = Time, group = 1, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~Subject_ordered, ncol = 6) +
  labs(
    title = "Individual Subject Decision Times Across Experimental Phases",
    x = "Time Point",
    y = "Time to Decision (seconds)"
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
  scale_color_manual(values = c("Control" = "navy", "Treatment" = "darkred"))


individual_time_plots

# save the plot
ggsave("individual_subject_time_plots.png", individual_time_plots, 
       width = 20, height = 24, dpi = 300, units = "in",
       bg = "white")
