# RelA and cRel nuclear level, batch corrected

# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(introdataviz)

# Specify the Excel file path
excel_file <- "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/Microscopy/VK014/Vk014_analysis_100cells_v3.xlsx"
folder <- "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/Microscopy/VK014/"

Westerns <- data.frame(
  timepoints = c(0, 8, 24, 48, 72),
  expRelA_CD40h = c(1, 9.570835624, 6.647391667, 8.812192916, 9.290313426),
  expRelA_CD40l = c(1,5.464728524,1.064973305,0.402735158,0.135871011),
  expRelA_BCR_CD40h = c(1,8.662575164,8.825318836,9.3780265,10.15210103),
  
  expcRel_CD40h = c(1,12.76960509,4.205811206,4.813371962,7.088068182),
  expcRel_CD40l = c(1,4.263613861,0.441887939,0,0.018986274),
  expcRel_BCR_CD40h = c(1,9.470775203,10.41122862,13.80555243,20.19675405)
)

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
  
  # Filter to only include rows for the desired channel
  combined_data <- combined_data %>% filter(Ch == channel) # RelA or cRel
  combined_data <- combined_data %>% filter(Group == group) # nucleus or whole cell, cytoplasm

  # Convert Time to numeric (you can remove the commented-out factor line if not needed)
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
  
  # Determine the two groups (assumes exactly 2 groups)
  group_levels <- unique(combined_data$Condition)
  if(length(group_levels) != 2) {
    stop("There must be exactly 2 groups for the split violin plot.")
  }
  combined_data$Condition <- factor(combined_data$Condition, levels = rev(group_levels))

  # Create the plot with two layers: one for each group
  p <- ggplot(combined_data, aes(x = Time, y = .data[[corrected_intensity_var]], fill = Condition, group = interaction(Condition, Time))) +
    introdataviz::geom_split_violin(alpha = 0.3, width = 12, trim = FALSE, color = scales::alpha("black", 0.1)) +
    # scale_y_continuous(limits = y_lim) +
    coord_cartesian(ylim = y_lim) +
    scale_x_continuous(breaks = c(0, 7, 24, 48, 72)) +
    labs(
      title = paste(title_suffix),
      y = paste(intensity_var, " (Corrected, AU)", sep = ""),
      x = "Time (hours)"
    ) +
    stat_summary(fun = mean, geom = "point",
                 aes(color = Condition), 
                 position = position_dodge(width = 0), 
                 size = 2) +
    stat_summary(fun = mean, geom = "line",
                 aes(color = Condition, group = Condition), 
                 position = position_dodge(width = 0), 
                 size = 0.75) +
    scale_fill_manual(values = plot_colors) +
    scale_color_manual(values = plot_colors) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          text = element_text(size = 16))
  
  # Add Western blot timepoints
  # transform_to_primary <- function(y) { y / 80 + 0.18 } # RelA transformation
  transform_to_primary <- function(y) { y / 90 + 0.2 } # cRel transformation
  
  # p <- p + geom_point(
  #   data = Westerns,
  #   aes(x = timepoints, y = transform_to_primary(expRelA_CD40l)),
  #   inherit.aes = FALSE,  # disables inheritance from the main plot
  #   shape = 24, fill = "#069ab6", alpha = 0.5, color="black", stroke=1, size = 3
  # )
  
  # p <- p + geom_point(
  #   data = Westerns,
  #   aes(x = timepoints, y = transform_to_primary(expcRel_CD40h)),
  #   inherit.aes = FALSE,  # disables inheritance from the main plot
  #   shape = 24, fill = "#005866", alpha = 0.5, color="black", stroke=1, size = 3
  # )
  # 
  # p <- p + geom_point(
  #   data = Westerns,
  #   aes(x = timepoints, y = transform_to_primary(expcRel_BCR_CD40h)),
  #   inherit.aes = FALSE,  # disables inheritance from the main plot
  #   shape = 24, fill = "#8e04b9", alpha = 0.5, color="black", stroke=1, size = 3
  # )
  
  # Save the plot as a PNG file
  ggsave(filename = plot_filename, plot = p, width = 6.4, height = 4.8, units="in")
  
  # Display the plot
  print(p)
}

# Custom colors
custom_colors <- c("#069ab6", "#e977ea", "#005866", "#8e04b9")
Int_var <- "Adj.Median"
BC_var <- "Batch correction"

# RelA_max <- c(0.03, 0.6)
# cRel_max <- c(0.05, 0.7)
RelA_max <- c(0.13, 0.43)
cRel_max <- c(0.1, 0.5)
postfix <- "_min"

normalization <- 1

# Plot the first 2 sheets on the same plot for RelA_Intensity
process_and_plot_combined(
  data_indices = 1:2,
  plot_colors = custom_colors[1:2],
  plot_filename = paste0(folder, "SplitViolin_FirstTwoSheets_RelA", postfix, ".png"),
  channel = "RelA",
  group = "nucleus",
  intensity_var = Int_var,
  title_suffix = "RelA Intensity",
  y_lim = RelA_max,
  normalization = normalization,
  Bcorrect = TRUE
)

# Plot the last 2 sheets on the same plot for RelA_Intensity
process_and_plot_combined(
  data_indices = 3:4,
  plot_colors = custom_colors[3:4],
  plot_filename = paste0(folder, "SplitViolin_LastTwoSheets_RelA", postfix, ".png"),
  channel = "RelA",
  group = "nucleus",
  intensity_var = Int_var,
  title_suffix = "RelA Intensity",
  y_lim = RelA_max,
  normalization = normalization,
  Bcorrect = TRUE
)

# Plot the first 2 sheets on the same plot for cRel_Intensity
process_and_plot_combined(
  data_indices = 1:2,
  plot_colors = custom_colors[1:2],
  plot_filename = paste0(folder, "SplitViolin_FirstTwoSheets_cRel", postfix, ".png"),
  channel = "cRel",
  group = "nucleus",
  intensity_var = Int_var,
  title_suffix = "cRel Intensity",
  y_lim = cRel_max,
  normalization = normalization,
  Bcorrect = TRUE
)

# Plot the last 2 sheets on the same plot for cRel_Intensity
process_and_plot_combined(
  data_indices = 3:4,
  plot_colors = custom_colors[3:4],
  plot_filename = paste0(folder, "SplitViolin_LastTwoSheets_cRel", postfix, ".png"),
  channel = "cRel",
  group = "nucleus",
  intensity_var = Int_var,
  title_suffix = "cRel Intensity",
  y_lim = cRel_max,
  normalization = normalization,
  Bcorrect = TRUE
)
