# Compatible with output from Helen's Julia code
library(ggplot2)
library(dplyr)
# ExpLineageFile2 <- "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/Bcell\ modeling\ project/Model\ fit\ evaluation/HX004\ CD40-low.csv"
ExpLineageFile2 <- "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/Bcell\ modeling\ project/Model\ fit\ evaluation/BCR-CD40-high7.csv"
Dose <- "high"
ExpLineageFile <- ""
if (Dose == "high") {
  ExpLineageFile <- "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/Bcell\ modeling\ project/Model\ fit\ evaluation/CD40-high7.csv"
} else if (Dose == "medium") {
  ExpLineageFile <- "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/Bcell\ modeling\ project/Model\ fit\ evaluation/CD40-med7.csv"
} else if (Dose == "low") {
  ExpLineageFile <- "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/Bcell\ modeling\ project/Model\ fit\ evaluation/CD40-low7.csv"
}
ExpLineage <- read.table(file=ExpLineageFile, sep=",", header = FALSE)
ExpLineage2 <- read.table(file=ExpLineageFile2, sep=",", header = FALSE)

# ExpLineage[4:11, 2:7] - ExpLineage2[4:11, 2:7]
# (ExpLineage[4:11, 2:7] - ExpLineage2[4:11, 2:7])^2
# sum((ExpLineage[4:11, 2:7] - ExpLineage2[4:11, 2:7])^2)
RMSD.generation <- sqrt(sum((ExpLineage[4:11, 2:7] - ExpLineage2[4:11, 2:7])^2))
RMSD.expansion15 <- sqrt(sum((ExpLineage[2, 2:7] - ExpLineage2[2, 2:7])^2 / 5 / max(ExpLineage[2, 2:7], ExpLineage2[2, 2:7])))
RMSD.expansion0 <- sqrt(sum((ExpLineage[3, 2:7] - ExpLineage2[3, 2:7])^2 / 6 / max(ExpLineage[3, 2:7], ExpLineage2[3, 2:7])))
RMSD.generation
RMSD.expansion15 + RMSD.expansion0
RMSD.total <- RMSD.generation + (RMSD.expansion15 + RMSD.expansion0) /2
RMSD.total