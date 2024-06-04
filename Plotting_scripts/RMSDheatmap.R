# Install and load necessary packages
library(readxl)
library(ggplot2)
library(cetcolor)
library(scales)

# Read data from Excel file
# xlsx_sheet = "HX007 fit w mismatch (96h)"
# xlsx_sheet = "8P 96h"
xlsx_sheet = "HX007 fit w mismatch (before96)"
excel_data <- read_excel("/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/BCR-CD40\ crosstalk\ project/RMSD\ fit.xlsx", 
                         sheet = xlsx_sheet)
dose_order <- c("none", "low", "med", "high")
# dose_order <- c("5P1", "5P3", "5P5", "5P8")
# dose_order <- c("9P1", "9P3", "9P5", "9P8")
# dose_order <- c("8P1", "8P3", "8P5", "8P8")
# dose_label_order <- c('1hr', '3hr', '5hr', '8hr')
dose_label_order <- dose_order

excel_data$exp_dose <- factor(excel_data$exp_dose, levels = dose_order)
excel_data$model_dose <- factor(excel_data$model_dose, levels = rev(dose_order))

# Create a heatmap
ggplot(data = excel_data, aes(x = exp_dose, y = model_dose, fill = generation, label = substring(generation, 1, 4))) +
  geom_tile() +
  geom_text(size = 6) +
  # scale_fill_distiller(palette = "RdYlBu") +
  scale_x_discrete(labels=dose_label_order) +
  scale_y_discrete(labels=rev(dose_label_order)) +
  scale_fill_gradientn(colours = cet_pal(100, name = "d1"), guide = "colourbar",
                       oob=squish, limits = c(0.1,1.3), breaks = c(0.1, 0.4, 0.7, 1.0, 1.3)) + 
  labs(x = "Experimental Dose", y = "Model Dose", fill = "Generation") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0)) +
  theme(legend.position="top",panel.background = element_blank(), legend.key.width = unit(1.8, "cm"), legend.key.height = unit(0.7, "cm"),
        axis.text = element_text(size = 18), axis.title = element_text(size = 20, face="bold"),
        legend.text = element_text(size=18), legend.title = element_text(size=20, face="bold")) +
  # theme(plot.margin = margin(c(5.5, 70, 5.5, 5.5), "points"))
  theme(plot.margin = margin(c(5.5, 40, 5.5, 5.5), "points"))


ggplot(data = excel_data, aes(x = exp_dose, y = model_dose, fill = pop_expansion, label = substring(pop_expansion, 1, 4))) +
  geom_tile() +
  geom_text(size = 6) +
  # scale_fill_distiller(palette = "RdYlBu") +
  scale_x_discrete(labels=dose_label_order) +
  scale_y_discrete(labels=rev(dose_label_order)) +
  scale_fill_gradientn(colours = cet_pal(100, name = "d1"), guide = "colourbar",
                       oob=squish, limits = c(0.4,1.9), breaks = c(0.4, 0.7, 1.0, 1.3, 1.6, 1.9)) +
                       # oob=squish, limits = c(0.1,1.3), breaks = c(0.1, 0.4, 0.7, 1.0, 1.3)) +  # 5P
                       # oob=squish, limits = c(0.1,1.0), breaks = c(0.4, 0.7, 1.0)) +  # 9P
  labs(x = "Experimental Dose", y = "Model Dose", fill = "Population Expansion") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0)) +
  theme(legend.position="top",panel.background = element_blank(), legend.key.width = unit(1.8, "cm"), legend.key.height = unit(0.7, "cm"), # 1.5 vs. 1.7 for key width
        axis.text = element_text(size = 18), axis.title = element_text(size = 20, face="bold"),
        legend.text = element_text(size=18), legend.title = element_text(size=20, face="bold")) +
  # theme(plot.margin = margin(c(5.5, 80, 5.5, 5.5), "points"))
  theme(plot.margin = margin(c(5.5, 40, 5.5, 5.5), "points"))
