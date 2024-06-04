# Compatible with output from Helen's Julia code
library(ggplot2)
library(dplyr)
countCells <-function(BcellLineage, mat, MAX_TIME) {
  cellCount <- nrow(BcellLineage)
  for (i in 1:cellCount) {
    current_gen <- BcellLineage$generation[i]
    if ((current_gen <= MAX_GEN) && (ceiling(BcellLineage$birthday[i]) <= MAX_TIME)) {
      for (j in ceiling(BcellLineage$birthday[i]):min(MAX_TIME, ceiling(BcellLineage$abs_fate_t[i]))) {
        mat[current_gen, j+1] = mat[current_gen, j+1] + 1
      }
    }
  }
  return(mat)
}

EXPDose <- "BCR-Clow"
# EXPDose <- "nostim"

in_dir <- "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/Bcell\ modeling\ project/RMSD\ evaluation/"
ExpLineageFile <- ""
seq_Dose = c("5P1", "5P3", "5P5", "5P8", 
            "9P1", "9P3", "9P5", "9P8", 
            "8P1", "8P3", "8P5", "8P8")
if (EXPDose == "high") {
  ExpLineageFile <- paste(in_dir, "CD40-high7.csv", sep="")
} else if (EXPDose == "med") {
  ExpLineageFile <- paste(in_dir, "CD40-med7.csv", sep="")
} else if (EXPDose == "low") {
  ExpLineageFile <- paste(in_dir, "CD40-low7.csv", sep="")
} else if (EXPDose == "nostim") {
  ExpLineageFile <- paste(in_dir, "CD40-nostim7.csv", sep="")
} else if (EXPDose == "BCR-Chigh") {
  ExpLineageFile <- paste(in_dir, "BCR-CD40-high7.csv", sep="")
} else if (EXPDose == "BCR-Clow") {
  ExpLineageFile <- paste(in_dir, "BCR-CD40-low7.csv", sep="")
} else if (EXPDose %in% seq_Dose) {
  ExpLineageFile <- paste(in_dir, EXPDose, ".csv", sep="")
}

ExpLineage <- read.table(file=ExpLineageFile, sep=",", header = FALSE)

eachFounder <- 125
Nfounder <- 0
MAX_GEN <- 7
MAX_TIME <- 96
binLength <- 1
cellCount <- nrow(BcellLineage)
mat = matrix(0, MAX_GEN, MAX_TIME+1)

# HomeDir = "/Users/helenhuang/Desktop/lineages_125_"
# HomeDir = "/Users/helenhuang/Desktop/Simulated\ output/lineages_125_"
HomeDir <- "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/BCR-CD40\ crosstalk\ project/withAICD/"
# HomeDir = "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/BCR-CD40\ crosstalk\ project/Sequential\ stim/"
# Modifier = "CD40F_LBCR"
# Modifier = "LCD40_HBCR_1hr_"
# Modifier = "LCD40_HBCR_8hr_"
Modifier = "Dose15_125cells_iter"
# Modifier = "nostimF_"
# Modifier = "CD40F_L" # after tuning cell fate modules (F=final)
# Modifier = "CD40_low" # before tuning cell fate modules
# Modifier = "nostim"

# Modifier = "lineages_125_HCD40_HBCR_8hr_" # 5P
# Modifier = "lineages_125_HCD40_LBCR_8hr_" # 9P
# Modifier = "lineages_125_LCD40_HBCR_8hr_" # 8P

# File1 = c(14, 15)
# File1 = c(11, 12, 13)
# File1 = c(4:6)
# File1 = c(1,2,4,6,7,8,55,56,57)
# File1 = c(1:12)
# File1 = c("101", "102", "103")
# File1 = c("104", "102", "103")
File1 = c(1:8)
# File1 = c(50:57)
for (m in 1:length(File1)) {
  FilePath = paste(HomeDir, Modifier, File1[m], ".txt",sep="")
  # FilePath = HomeDir
  if (file.exists(FilePath)) {
    BcellLineage <- read.table(file=FilePath, sep="\t", header = TRUE)
    mat = countCells(BcellLineage, mat, MAX_TIME)
    Nfounder = Nfounder+eachFounder
  }
}

# sumstats of model output
# timepoints <- mat[, c(24, 38, 50, 72, 96, 120)]
timepoints <- mat[, c(24, 38, 50, 72, 96)]
Ntimepoints <- 5

timepointsSums <- colSums(timepoints)
timepointsP <- t(t(timepoints) / timepointsSums)
popExpansion0 <- timepointsSums / Nfounder
popExpansion15 <- timepointsSums / timepointsSums[1]

RMSD.generation <- sqrt(sum((ExpLineage[4:10, 2:(1+Ntimepoints)] - timepointsP)^2))
RMSD.expansion15 <- sqrt(sum((ExpLineage[2, 2:(1+Ntimepoints)] - popExpansion15)^2 / (Ntimepoints-1) / max(ExpLineage[2, 2:(1+Ntimepoints)])))
RMSD.expansion0 <- sqrt(sum((ExpLineage[3, 2:(1+Ntimepoints)] - popExpansion0)^2 / Ntimepoints / max(ExpLineage[3, 2:(1+Ntimepoints)])))
RMSD.generation
RMSD.expansion15 + RMSD.expansion0
RMSD.total <- RMSD.generation + (RMSD.expansion15 + RMSD.expansion0)
RMSD.total


# generation vs. time area plot
generation <- rep(1:(MAX_GEN), each=MAX_TIME+1)
generation <- as.character(generation)
time <- rep(0:MAX_TIME,times=(MAX_GEN))
cells <- as.vector(t(mat)) / Nfounder
data <- data.frame(time, cells, generation)
data$generation <- factor(data$generation , levels=c("7", "6", "5", "4", "3", "2", "1") )
# stacked area chart
ggplot(data, aes(x=time, y=cells, fill=generation)) + geom_area() + xlab("Time (hour)") +
  scale_fill_manual(values = c("grey20", "grey30", "grey40", "grey50", "grey59", "grey70", "grey87")) + 
  theme_bw() + 
  scale_x_continuous(limits = c(0, MAX_TIME), expand = c(0.02, 0.02), breaks = seq(0, MAX_TIME, 24)) +
  scale_y_continuous(limits = c(0, 3.0), expand = c(0.02, 0.02)) +
  # scale_fill_grey(start = 0, end = .84) + theme_bw() + ylim(0, 420) + xlim(0, 120) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




