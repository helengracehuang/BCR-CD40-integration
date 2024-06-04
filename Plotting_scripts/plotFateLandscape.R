library(ggplot2)
library(latticeExtra)
library(dplyr)
library(cetcolor)

MAX_TIME <- 120
Nfounder <- 125

survivalT1 <- 24
survivalT2 <- 96
proliferateT <- 84

# COLOR_SCHEME = terrain.colors(100)
COLOR_SCHEME = cet_pal(100, name = "l8", alpha = 1)
COLOR_SCHEME_TR = cet_pal(100, name = "l8", alpha = 0)
COLOR_SCHEME_DIV = cet_pal(100, name = "d1")
SMOOTHNESS = 100

# Function to process files
process_file <- function(file_path) {
  BcellLineages <- read.table(file=file_path, sep="\t", header = TRUE)
  cellCount <- nrow(BcellLineages)
  cells <- rep(0, MAX_TIME) # live cells
  Pcells <- rep(0, MAX_TIME) # proliferating founder cells
  for (i in 1:cellCount) {
    for (j in ceiling(BcellLineages$birthday[i]):ceiling(BcellLineages$abs_fate_t[i])) {
      cells[j] = cells[j] + 1
    }
  }
  for (i in 1:Nfounder) {
    for (j in ceiling(BcellLineages$abs_fate_t[i]):MAX_TIME) {
      if (BcellLineages$fate[i] == 2) {
        Pcells[j] = Pcells[j] + 1
      }
    }
  }
  return(list("live_cells" = cells, "pf_cells" = Pcells))
}

modifier <- "withAICD"
# modifier <- "noAICD"
datadir <- file.path(paste("/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/BCR-CD40\ crosstalk\ project/", modifier, sep=""))
outdir <- file.path(paste("/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/Raw\ plots\ &\ figures/Simulation/", modifier, sep=""))
dataFiles <- dir(file.path(datadir), pattern = ".txt$",
                full.names = TRUE, recursive = TRUE)

# extract the dose number in the middle of the file names using regex
if (modifier == "withAICD") {
  doseList <- gsub('^..............................................................................................|...................$', '', dataFiles)
} else if (modifier == "noAICD") {
  doseList <- gsub('^............................................................................................|...................$', '', dataFiles)
}
  
doseList <- as.integer(doseList)
dataTable <- data.frame(doseList, dataFiles)
dataTable$CD40_delay <- ifelse(dataTable$doseList <= 25, 0,
                             ifelse(dataTable$doseList > 25 & dataTable$doseList <= 50, 1, 8))
dataTable$CD40L_dose <- ((dataTable$doseList-1) %% 25) %/% 5  # integer division
dataTable$Antigen_dose <- ((dataTable$doseList-1) %% 25) %% 5  # modulo

AllCells <- list()
LiveCells <- list()
ProlifFounderCells <- list()

liveT1 <- c()
liveT2 <- c()
prolifT <- c()

for (file_path in dataTable$dataFiles) {
  AllCells[[file_path]] <- process_file(file_path)
  LiveCells <- AllCells[[file_path]]$live_cells
  ProlifFounderCells <- AllCells[[file_path]]$pf_cells
  liveT1 <- c(liveT1, LiveCells[survivalT1])
  liveT2 <- c(liveT2, LiveCells[survivalT2])
  prolifT <- c(prolifT, ProlifFounderCells[proliferateT])
}

dataTable$liveT1 <- liveT1
dataTable$liveT2 <- liveT2
dataTable$prolifT <- prolifT

# Collapse by averaging values of column1 across rows with the same column2 and column3
collapsedDataTable <- dataTable %>%
  group_by(CD40_delay, CD40L_dose, Antigen_dose) %>%
  summarize(liveT1 = mean(liveT1), liveT2 = mean(liveT2), prolifT = mean(prolifT))

data_noDelay <- subset(collapsedDataTable, CD40_delay == 0)
data_1hrDelay <- subset(collapsedDataTable, CD40_delay == 1)
data_8hrDelay <- subset(collapsedDataTable, CD40_delay == 8)

if (modifier == "withAICD") {
  data_noDelay[10,4:6] = data_noDelay[15,4:6]
}
# =============================================================
# # plotting cell survival landscape at survivalT1 hrs
# plot <- levelplot(liveT1 / Nfounder ~ CD40L_dose * Antigen_dose, data_noDelay, col.regions = COLOR_SCHEME, main = paste("Cell survival rate at ", survivalT1, "hrs",sep=""),
#           panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
#           at=seq(0.15, 0.7, length=SMOOTHNESS), colorkey=list( at=seq(0.15, 0.7, length=SMOOTHNESS), col=COLOR_SCHEME))  +
#   layer_(panel.2dsmoother(..., n = 200))
# plot
# png(paste(outdir, "/Terrain_live24_noDelay.png", sep=""), width = 2000, height = 1900, res = 300)
# print(plot)
# dev.off()
# 
# # plotting cell proliferative landscape at prolifT hrs
# plot <- levelplot(prolifT / liveT1 ~ CD40L_dose * Antigen_dose, data_noDelay, col.regions = COLOR_SCHEME, main = paste("% Cumulative proliferative cell at ", proliferateT, "hrs",sep=""),
#                   panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
#                   at=seq(0, 0.35, length=SMOOTHNESS), colorkey=list( at=seq(0, 0.35, length=SMOOTHNESS), col=COLOR_SCHEME) ) +
# layer_(panel.2dsmoother(..., n = 500, method = "loess"))
# plot
# png(paste(outdir, "/Terrain_prolifT_noDelay.png", sep=""), width = 2000, height = 1900, res = 300)
# print(plot)
# dev.off()
# 
# # plotting cell fold change landscape at survivalT2 hrs
# plot <- levelplot(liveT2 / Nfounder ~ CD40L_dose * Antigen_dose, data_noDelay, col.regions = COLOR_SCHEME, main = paste("Population fold-change at ", survivalT2, "hrs",sep=""),
#                   panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose", 
#                   at=seq(0, 2.5, length=SMOOTHNESS), colorkey=list( at=seq(0, 2.5, length=SMOOTHNESS), col=COLOR_SCHEME) ) +
#   layer_(panel.2dsmoother(..., n = 500, method = "loess"))
# plot
# png(paste(outdir, "/Terrain_liveT2_noDelay.png", sep=""), width = 2000, height = 1900, res = 300)
# print(plot)
# dev.off()
# 
# # par(bg="transparent")
# plot <- contourplot(liveT2 / Nfounder ~ CD40L_dose * Antigen_dose, data_noDelay, col.regions = COLOR_SCHEME, 
#                     main = paste("Population fold-change at ", survivalT2, "hrs",sep=""), pretty = TRUE,
#                     labels=list(col = "white", cex = 2),label.style="align",
#                     panel = panel.levelplot.points, cex = 1.2, pch=1, xlab="CD40L dose", ylab="Antigen dose",
#                   # at=seq(0, 3.5, length=SMOOTHNESS), 
#                   alpha.regions=0, colorkey=list( at=seq(0, 2.5, length=SMOOTHNESS), col=COLOR_SCHEME_TR), 
#                   cuts = 4) + layer_(panel.2dsmoother(..., n = 500, method = "loess", lwd=2, col = "white"))
# plot
# png(paste(outdir, "/Terrain_liveT2_noDelay_contour.png", sep=""), bg="transparent", width = 2000, height = 1900, res = 300)
# print(plot)
# dev.off()



# # =============================================================
# # plotting cell survival landscape at survivalT1 hrs with 8hr delay
# plot <- levelplot(liveT1 / Nfounder ~ CD40L_dose * Antigen_dose, data_8hrDelay, col.regions = COLOR_SCHEME, main = paste("Cell survival rate at ", survivalT1, "hrs",sep=""),
#                   panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
#                   at=seq(0, 0.7, length=SMOOTHNESS), colorkey=list( at=seq(0, 0.7, length=SMOOTHNESS), col=COLOR_SCHEME))  +
#   layer_(panel.2dsmoother(..., n = 200))
# plot
# png(paste(outdir, "/Terrain_live24_8hrDelay.png", sep=""), width = 2000, height = 1900, res = 300)
# print(plot)
# dev.off()
# 
# # plotting cell proliferative landscape at prolifT hrs
# plot <- levelplot(prolifT / liveT1 ~ CD40L_dose * Antigen_dose, data_8hrDelay, col.regions = COLOR_SCHEME, main = paste("% Cumulative proliferative cell at ", proliferateT, "hrs",sep=""),
#                   panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
#                   at=seq(-0.025, 0.375, length=SMOOTHNESS), colorkey=list( at=seq(0, 0.375, length=SMOOTHNESS), col=COLOR_SCHEME) ) +
#   layer_(panel.2dsmoother(..., n = 200))
# plot
# png(paste(outdir, "/Terrain_prolifT_8hrDelay.png", sep=""), width = 2000, height = 1900, res = 300)
# print(plot)
# dev.off()
# 
# # plotting cell fold change landscape at survivalT2 hrs
# plot <- levelplot(liveT2 / Nfounder ~ CD40L_dose * Antigen_dose, data_8hrDelay, col.regions = COLOR_SCHEME, main = paste("Population fold-change at ", survivalT2, "hrs",sep=""),
#                   panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
#                   at=seq(0, 2.3, length=SMOOTHNESS), colorkey=list( at=seq(0, 2.3, length=SMOOTHNESS), col=COLOR_SCHEME) ) +
#   layer_(panel.2dsmoother(..., n = 500, method = "loess"))
# plot
# png(paste(outdir, "/Terrain_liveT2_8hrDelay.png", sep=""), width = 2000, height = 1900, res = 300)
# print(plot)
# dev.off()
# 
# 
# 
# 
# # =============================================================
# # plotting cell survival landscape at survivalT1 hrs with 1hr delay
# plot <- levelplot(liveT1 / Nfounder ~ CD40L_dose * Antigen_dose, data_1hrDelay, col.regions = COLOR_SCHEME, main = paste("Cell survival rate at ", survivalT1, "hrs",sep=""),
#                   panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
#                   at=seq(0, 0.7, length=SMOOTHNESS), colorkey=list( at=seq(0, 0.7, length=SMOOTHNESS), col=COLOR_SCHEME))  +
#   layer_(panel.2dsmoother(..., n = 200))
# plot
# png(paste(outdir, "/Terrain_live24_1hrDelay.png", sep=""), width = 2000, height = 1900, res = 300)
# print(plot)
# dev.off()
# 
# # plotting cell proliferative landscape at prolifT hrs
# plot <- levelplot(prolifT / liveT1 ~ CD40L_dose * Antigen_dose, data_1hrDelay, col.regions = COLOR_SCHEME, main = paste("% Cumulative proliferative cell at ", proliferateT, "hrs",sep=""),
#                   panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
#                   at=seq(0, 0.375, length=SMOOTHNESS), colorkey=list( at=seq(0, 0.375, length=SMOOTHNESS), col=COLOR_SCHEME))  + #expand = c(0.01, 0.01) +
#   layer_(panel.2dsmoother(..., n = 200))
# plot
# png(paste(outdir, "/Terrain_prolifT_1hrDelay.png", sep=""), width = 2000, height = 1900, res = 300)
# print(plot)
# dev.off()
# 
# 
# # plotting cell fold change landscape at survivalT2 hrs
# plot <- levelplot(liveT2 / Nfounder ~ CD40L_dose * Antigen_dose, data_1hrDelay, col.regions = COLOR_SCHEME, main = paste("Population fold-change at ", survivalT2, "hrs",sep=""),
#                   panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
#                   at=seq(0, 2.3, length=SMOOTHNESS), colorkey=list( at=seq(0, 2.3, length=SMOOTHNESS), col=COLOR_SCHEME) ) +
#   layer_(panel.2dsmoother(..., n = 500, method = "loess"))
# plot
# png(paste(outdir, "/Terrain_liveT2_1hrDelay.png", sep=""), width = 2000, height = 1900, res = 300)
# print(plot)
# dev.off()







# plotting differences between 2 plots
data_difference <- cbind(data_8hrDelay[, 1:3], (data_8hrDelay[, 4:6]-data_1hrDelay[, 4:6]))
plot <- levelplot(liveT1 / Nfounder ~ CD40L_dose * Antigen_dose, data_difference, col.regions = COLOR_SCHEME_DIV, main = paste("Cell survival rate at ", survivalT1, "hrs",sep=""),
                  panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
                  at=seq(-0.5, 0.5, length=SMOOTHNESS), colorkey=list( at=seq(-0.5, 0.5, length=SMOOTHNESS), col=COLOR_SCHEME_DIV))  +
  layer_(panel.2dsmoother(..., n = 200))
plot
png(paste(outdir, "/Terrain_live24_difference.png", sep=""), width = 2000, height = 1900, res = 300)
print(plot)
dev.off()

plot <- levelplot((data_8hrDelay$prolifT / data_8hrDelay$liveT1 - data_1hrDelay$prolifT / data_1hrDelay$liveT1) ~ CD40L_dose * Antigen_dose, data_difference, col.regions = COLOR_SCHEME_DIV, main = paste("% Cumulative proliferative cell at ", proliferateT, "hrs",sep=""),
                  panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
                  at=seq(-0.25, 0.25, length=SMOOTHNESS), colorkey=list( at=seq(-0.25, 0.25, length=SMOOTHNESS), col=COLOR_SCHEME_DIV))  + #expand = c(0.01, 0.01) +
  layer_(panel.2dsmoother(..., n = 200))
plot
png(paste(outdir, "/Terrain_prolifT_difference.png", sep=""), width = 2000, height = 1900, res = 300)
print(plot)
dev.off()
plot <- levelplot(liveT2 / Nfounder ~ CD40L_dose * Antigen_dose, data_difference, col.regions = COLOR_SCHEME_DIV, main = paste("Population fold-change at ", survivalT2, "hrs",sep=""),
                  panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
                  at=seq(-1.5, 1.5, length=SMOOTHNESS), colorkey=list( at=seq(-1.5, 1.5, length=SMOOTHNESS), col=COLOR_SCHEME_DIV) ) +
  layer_(panel.2dsmoother(..., n = 500, method = "loess"))
plot
png(paste(outdir, "/Terrain_liveT2_difference.png", sep=""), width = 2000, height = 1900, res = 300)
print(plot)
dev.off()





# plotting differences between 2 plots
data_difference <- cbind(data_8hrDelay[, 1:3], (data_8hrDelay[, 4:6]/data_1hrDelay[, 4:6]))
plot <- levelplot(log10(liveT1) ~ CD40L_dose * Antigen_dose, data_difference, col.regions = COLOR_SCHEME_DIV, main = paste("Cell survival rate at ", survivalT1, "hrs",sep=""),
                  panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
                  at=seq(-0.8, 0.8, length=SMOOTHNESS), colorkey=list( at=seq(-0.8, 0.8, length=SMOOTHNESS), col=COLOR_SCHEME_DIV))  +
  layer_(panel.2dsmoother(..., n = 200))
plot
png(paste(outdir, "/Terrain_live24_ratio.png", sep=""), width = 2000, height = 1900, res = 300)
print(plot)
dev.off()

plot <- levelplot(log10(prolifT/liveT1) ~ CD40L_dose * Antigen_dose, data_difference, col.regions = COLOR_SCHEME_DIV, main = paste("% Cumulative proliferative cell at ", proliferateT, "hrs",sep=""),
                  panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
                  at=seq(-0.8, 0.8, length=SMOOTHNESS), colorkey=list( at=seq(-0.8, 0.8, length=SMOOTHNESS), col=COLOR_SCHEME_DIV))  + #expand = c(0.01, 0.01) +
  layer_(panel.2dsmoother(..., n = 200))
plot
png(paste(outdir, "/Terrain_prolifT_ratio.png", sep=""), width = 2000, height = 1900, res = 300)
print(plot)
dev.off()
plot <- levelplot(log10(liveT2) ~ CD40L_dose * Antigen_dose, data_difference, col.regions = COLOR_SCHEME_DIV, main = paste("Population fold-change at ", survivalT2, "hrs",sep=""),
                  panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
                  at=seq(-0.8, 0.8, length=SMOOTHNESS), colorkey=list( at=seq(-0.8, 0.8, length=SMOOTHNESS), col=COLOR_SCHEME_DIV) ) +
  layer_(panel.2dsmoother(..., n = 500, method = "loess"))
plot
png(paste(outdir, "/Terrain_liveT2_ratio.png", sep=""), width = 2000, height = 1900, res = 300)
print(plot)
dev.off()
