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

modifier <- "noAICD"
datadir <- file.path(paste("/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/BCR-CD40\ crosstalk\ project/", modifier, sep=""))
outdir <- file.path(paste("/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/Raw\ plots\ &\ figures/Simulation/", modifier, sep=""))
dataFiles <- dir(file.path(datadir), pattern = ".txt$",
                 full.names = TRUE, recursive = TRUE)

# extract the dose number in the middle of the file names using regex
doseList <- gsub('^............................................................................................|...................$', '', dataFiles)

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

data_noDelay_noAICD <- subset(collapsedDataTable, CD40_delay == 0)
data_1hrDelay_noAICD <- subset(collapsedDataTable, CD40_delay == 1)
data_8hrDelay_noAICD <- subset(collapsedDataTable, CD40_delay == 8)




modifier <- "withAICD"
datadir <- file.path(paste("/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/BCR-CD40\ crosstalk\ project/", modifier, sep=""))
outdir <- file.path(paste("/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/Raw\ plots\ &\ figures/Simulation/", modifier, sep=""))
dataFiles <- dir(file.path(datadir), pattern = ".txt$",
                 full.names = TRUE, recursive = TRUE)

# extract the dose number in the middle of the file names using regex
doseList <- gsub('^..............................................................................................|...................$', '', dataFiles)


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

data_noDelay_withAICD <- subset(collapsedDataTable, CD40_delay == 0)
data_1hrDelay_withAICD <- subset(collapsedDataTable, CD40_delay == 1)
data_8hrDelay_withAICD <- subset(collapsedDataTable, CD40_delay == 8)

data_noDelay_withAICD[10,4:6] = data_noDelay_withAICD[15,4:6]


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





# plotting differences between 2 plots
data_difference <- cbind(data_noDelay_withAICD[, 1:3], (data_noDelay_withAICD[, 4:6]-data_noDelay_noAICD[, 4:6]))
plot <- levelplot(liveT1 / Nfounder ~ CD40L_dose * Antigen_dose, data_difference, col.regions = COLOR_SCHEME_DIV, main = paste("Cell survival rate at ", survivalT1, "hrs",sep=""),
                  panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
                  at=seq(-0.5, 0.5, length=SMOOTHNESS), colorkey=list( at=seq(-0.5, 0.5, length=SMOOTHNESS), col=COLOR_SCHEME_DIV))  +
  layer_(panel.2dsmoother(..., n = 200))
plot
png(paste(outdir, "/Terrain_live24_AICDdifference.png", sep=""), width = 2000, height = 1900, res = 300)
print(plot)
dev.off()

plot <- levelplot((data_noDelay_withAICD$prolifT / data_noDelay_withAICD$liveT1 - data_noDelay_noAICD$prolifT / data_noDelay_noAICD$liveT1) ~ CD40L_dose * Antigen_dose, data_difference, col.regions = COLOR_SCHEME_DIV, main = paste("% Cumulative proliferative cell at ", proliferateT, "hrs",sep=""),
                  panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
                  at=seq(-0.25, 0.25, length=SMOOTHNESS), colorkey=list( at=seq(-0.25, 0.25, length=SMOOTHNESS), col=COLOR_SCHEME_DIV))  + #expand = c(0.01, 0.01) +
  layer_(panel.2dsmoother(..., n = 200))
plot
png(paste(outdir, "/Terrain_prolifT_AICDdifference.png", sep=""), width = 2000, height = 1900, res = 300)
print(plot)
dev.off()
plot <- levelplot(liveT2 / Nfounder ~ CD40L_dose * Antigen_dose, data_difference, col.regions = COLOR_SCHEME_DIV, main = paste("Population fold-change at ", survivalT2, "hrs",sep=""),
                  panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
                  at=seq(-1.5, 1.5, length=SMOOTHNESS), colorkey=list( at=seq(-1.5, 1.5, length=SMOOTHNESS), col=COLOR_SCHEME_DIV) ) +
  layer_(panel.2dsmoother(..., n = 500, method = "loess"))
plot
png(paste(outdir, "/Terrain_liveT2_AICDdifference.png", sep=""), width = 2000, height = 1900, res = 300)
print(plot)
dev.off()






# plotting ratios between 2 plots
data_difference <- cbind(data_noDelay_withAICD[, 1:3], (data_noDelay_withAICD[, 4:6]/data_noDelay_noAICD[, 4:6]))
plot <- levelplot(log10(liveT1) ~ CD40L_dose * Antigen_dose, data_difference, col.regions = COLOR_SCHEME_DIV, main = paste("Cell survival rate at ", survivalT1, "hrs",sep=""),
                  panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
                  at=seq(-0.5, 0.5, length=SMOOTHNESS), colorkey=list( at=seq(-0.5, 0.5, length=SMOOTHNESS), col=COLOR_SCHEME_DIV))  +
  layer_(panel.2dsmoother(..., n = 200))
plot
png(paste(outdir, "/Terrain_live24_AICDratio.png", sep=""), width = 2000, height = 1900, res = 300)
print(plot)
dev.off()

plot <- levelplot(log10(prolifT/liveT1) ~ CD40L_dose * Antigen_dose, data_difference, col.regions = COLOR_SCHEME_DIV, main = paste("% Cumulative proliferative cell at ", proliferateT, "hrs",sep=""),
                  panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
                  at=seq(-0.5, 0.5, length=SMOOTHNESS), colorkey=list( at=seq(-0.5, 0.5, length=SMOOTHNESS), col=COLOR_SCHEME_DIV))  + #expand = c(0.01, 0.01) +
  layer_(panel.2dsmoother(..., n = 200))
plot
png(paste(outdir, "/Terrain_prolifT_AICDratio.png", sep=""), width = 2000, height = 1900, res = 300)
print(plot)
dev.off()
plot <- levelplot(log10(liveT2) ~ CD40L_dose * Antigen_dose, data_difference, col.regions = COLOR_SCHEME_DIV, main = paste("Population fold-change at ", survivalT2, "hrs",sep=""),
                  panel = panel.levelplot.points, cex = 1.2, xlab="CD40L dose", ylab="Antigen dose",
                  at=seq(-0.5, 0.5, length=SMOOTHNESS), colorkey=list( at=seq(-0.5, 0.5, length=SMOOTHNESS), col=COLOR_SCHEME_DIV) ) +
  layer_(panel.2dsmoother(..., n = 500, method = "loess"))
plot
png(paste(outdir, "/Terrain_liveT2_AICDratio.png", sep=""), width = 2000, height = 1900, res = 300)
print(plot)
dev.off()