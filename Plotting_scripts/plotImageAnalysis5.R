# RelA and cRel nuclear level, batch corrected

# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)

# Specify the Excel file path
excel_file <- "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/Microscopy/VK014/Vk014_analysis_100cells_v3.xlsx"
folder <- "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/Microscopy/VK014/"

# Read all sheet names from the Excel file
sheet_names <- c("CD40_low_BCR_high", "CD40_low", "CD40_high_BCR_high", "CD40_high")

bc_table <- read_excel(excel_file, sheet = "Batch_norm")

# Read data from all sheets into a list and add a 'Group' column
data_list <- lapply(seq_along(sheet_names), function(i) {
  data <- read_excel(excel_file, sheet = sheet_names[i])
  data$Condition <- sheet_names[i]
  return(data)
})

# Function to process and plot combined data for a specified intensity variable
process_and_plot_combined <- function(data_indices, plot_colors, plot_filename, channel, group, intensity_var, title_suffix, y_lim, normalization, Bcorrect) {
  # Combine data frames for the specified indices
  combined_data <- bind_rows(data_list[data_indices])
  
  ### Edit #1: NEW filtering step based on the new "Ch" column.
  # Filter to only include rows for the desired channel
  combined_data <- combined_data %>% filter(Ch == channel) # RelA or cRel
  combined_data <- combined_data %>% filter(Group == group) # nucleus or whole cell, cytoplasm
  ### End Edit #1
  
  # Convert Time to numeric (you can remove the commented-out factor line if not needed)
  # combined_data$Time <- factor(combined_data$Time, levels = c(0, 7, 24, 48, 72))
  combined_data$Time <- as.numeric(as.character(combined_data$Time))
  
  # Compute the corrected intensity
  corrected_intensity_var <- paste0(channel, "_Corrected")
  
  if (Bcorrect) {
    # join the lookup table 'bc_table' based on Time, Ch, and Condition
    combined_data <- combined_data %>%
      left_join(bc_table, by = c("Time", "Ch", "Condition")) %>%
      # Then calculate the corrected intensity by dividing the measured intensity
      # by the lookup value from 'bc_table' (named 'bc_value').
      mutate(!!corrected_intensity_var := .data[[intensity_var]] / bc_value)
  }
  else {
    combined_data <- combined_data %>%
      mutate(!!corrected_intensity_var := .data[[intensity_var]])
  }
  # Calculate mean and standard error for each Time and Group
  summary_data <- combined_data %>%
    group_by(Time, Condition) %>%
    summarise(
      mean = mean(.data[[corrected_intensity_var]], na.rm = TRUE) / normalization,
      sd = sd(.data[[corrected_intensity_var]], na.rm = TRUE) / normalization,
      n = n(),
      se = sd / sqrt(n)
    ) %>%
    ungroup()
  
  # Plot the trajectory with error bars
  p <- ggplot(summary_data, aes(x = Time, y = mean, color = Condition, group = Condition)) +
    geom_line(size = 1, show.legend = FALSE) +
    geom_point(size = 2, show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), size = 1, width = 3, show.legend = FALSE) +
    scale_color_manual(values = plot_colors) +
    scale_x_continuous(breaks = c(0, 7, 24, 48, 72)) +
    scale_y_continuous(limits = y_lim) +
    labs(
      title = paste(title_suffix),
      y = paste(intensity_var, " (AU)", sep = ""),
      x = "Time (hours)"
    ) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(
      legend.title = element_blank(),
      text = element_text(size = 16)
    )
  
  # Save the plot as a PNG file
  ggsave(filename = plot_filename, plot = p, width = 6.4, height = 4.8, units="in")
  
  # Display the plot
  print(p)
}

custom_colors <- c("#005866", "#8e04b9", "#069ab6", "#e977ea")

# Plot all 4 sheets for area
process_and_plot_combined(
  data_indices = 1:4,
  plot_colors = custom_colors[1:4],
  plot_filename = paste0(folder, "Trajectory_Area", postfix, ".png"),
  channel = "cRel",
  group = "whole cell",
  intensity_var = "Area",
  title_suffix = "Cell Area",
  y_lim = c(0.5, 1.8),
  normalization = 1/8,
  Bcorrect = FALSE
)