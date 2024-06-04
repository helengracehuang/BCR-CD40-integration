# Compatible with output from Helen's Julia code (BCR pulsing CD40 stimulation)
library(ggplot2)
library(dplyr)
countCells <-function(BcellLineage, cells) {
  cellCount <- nrow(BcellLineage)
  for (i in 1:cellCount) {
    current_gen <- BcellLineage$generation[i]
    if (current_gen <= MAX_GEN) {
      for (j in ceiling(BcellLineage$birthday[i]):ceiling(BcellLineage$abs_fate_t[i])) {
        cells[j] = cells[j] + 1
      }
    }
  }
  return(cells)
}

HomeDir = "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/BCR-CD40\ crosstalk\ project/Sequential\ stim/lineages_125_HCD40_HBCR_"
File1 = c("107", "108", "109") # version of simulation
BcellLineage1.1 <- read.table(file=paste(HomeDir, "1hr_", File1[1], ".txt",sep=""), sep="\t", header = TRUE)
BcellLineage1.2 <- read.table(file=paste(HomeDir, "1hr_", File1[2], ".txt",sep=""), sep="\t", header = TRUE)
BcellLineage1.3 <- read.table(file=paste(HomeDir, "1hr_", File1[3], ".txt",sep=""), sep="\t", header = TRUE)
File3 = c("107", "108", "109")
BcellLineage3.1 <- read.table(file=paste(HomeDir, "3hr_", File3[1], ".txt",sep=""), sep="\t", header = TRUE)
BcellLineage3.2 <- read.table(file=paste(HomeDir, "3hr_", File3[2], ".txt",sep=""), sep="\t", header = TRUE)
BcellLineage3.3 <- read.table(file=paste(HomeDir, "3hr_", File3[3], ".txt",sep=""), sep="\t", header = TRUE)
# HomeDir = "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/BCR-CD40\ crosstalk\ project/Sequential\ stim/lineages_125_LCD40_HBCR_"
File5 = c("101", "102", "103")
BcellLineage5.1 <- read.table(file=paste(HomeDir, "5hr_", File5[1], ".txt",sep=""), sep="\t", header = TRUE)
BcellLineage5.2 <- read.table(file=paste(HomeDir, "5hr_", File5[2], ".txt",sep=""), sep="\t", header = TRUE)
BcellLineage5.3 <- read.table(file=paste(HomeDir, "5hr_", File5[3], ".txt",sep=""), sep="\t", header = TRUE)
File8 = c("104", "102", "103")
BcellLineage8.1 <- read.table(file=paste(HomeDir, "8hr_", File8[1], ".txt",sep=""), sep="\t", header = TRUE)
BcellLineage8.2 <- read.table(file=paste(HomeDir, "8hr_", File8[2], ".txt",sep=""), sep="\t", header = TRUE)
BcellLineage8.3 <- read.table(file=paste(HomeDir, "8hr_", File8[3], ".txt",sep=""), sep="\t", header = TRUE)

BcellLineages <- list(BcellLineage1.1, BcellLineage1.2, BcellLineage1.3, 
                   BcellLineage3.1, BcellLineage3.2, BcellLineage3.3, 
                   BcellLineage5.1, BcellLineage5.2, BcellLineage5.3, 
                   BcellLineage8.1, BcellLineage8.2, BcellLineage8.3)
Nfounder <- 125
MAX_GEN <- 8
MAX_TIME <- 96
binLength <- 1
replicates <- 3
cells <- rep(0, MAX_TIME)

# for (m in 1:length(File1)) {
#   FilePath = paste(HomeDir, Modifier, File1[m], ".txt",sep="")
#   if (file.exists(FilePath)) {
#     BcellLineage <- read.table(file=FilePath, sep="\t", header = TRUE)
#     cells = countCells(BcellLineage, cells)
#     Nfounder = Nfounder+eachFounder
#   }
# }

AllCells <- list()
for (m in 1:length(BcellLineages)) {
  cellCount <- nrow(BcellLineages[[m]])
  cells <- rep(0, MAX_TIME)
  for (i in 1:cellCount) {
    if (ceiling(BcellLineages[[m]]$birthday[i]) <= MAX_TIME) {
      for (j in ceiling(BcellLineages[[m]]$birthday[i]):min(MAX_TIME, ceiling(BcellLineages[[m]]$abs_fate_t[i]))) {
        cells[j] = cells[j] + 1
      }
    }
  }
  AllCells[[length(AllCells) + 1]] <- cells
}

cells1.mean = (AllCells[[1]] + AllCells[[2]] + AllCells[[3]]) /replicates
cells1.sd = sqrt(((AllCells[[1]] - cells1.mean)^2 + (AllCells[[2]] - cells1.mean)^2 + (AllCells[[3]] - cells1.mean)^2) /replicates)
# cells1.margin <- qt(0.975,df=replicates-1)*cells1.sd/sqrt(replicates)
cells3.mean = (AllCells[[3+1]] + AllCells[[3+2]] + AllCells[[3+3]]) /replicates
cells3.sd = sqrt(((AllCells[[3+1]] - cells3.mean)^2 + (AllCells[[3+2]] - cells3.mean)^2 + (AllCells[[3+3]] - cells3.mean)^2) /replicates)

cells5.mean = (AllCells[[6+1]] + AllCells[[6+2]] + AllCells[[6+3]]) /replicates
cells5.sd = sqrt(((AllCells[[6+1]] - cells5.mean)^2 + (AllCells[[6+2]] - cells5.mean)^2 + (AllCells[[6+3]] - cells5.mean)^2) /replicates)

cells8.mean = (AllCells[[9+1]] + AllCells[[9+2]] + AllCells[[9+3]]) /replicates
cells8.sd = sqrt(((AllCells[[9+1]] - cells8.mean)^2 + (AllCells[[9+2]] - cells8.mean)^2 + (AllCells[[9+3]] - cells8.mean)^2) /replicates)

# cells10.mean = (AllCells[[12+1]] + AllCells[[12+2]] + AllCells[[12+3]]) /replicates
# cells10.sd = sqrt(((AllCells[[12+1]] - cells10.mean)^2 + (AllCells[[12+2]] - cells10.mean)^2 + (AllCells[[12+3]] - cells10.mean)^2) /replicates)

cells1.mean = cells1.mean / Nfounder
cells1.sd = cells1.sd / Nfounder
cells3.mean = cells3.mean / Nfounder
cells3.sd = cells3.sd / Nfounder
cells5.mean = cells5.mean / Nfounder
cells5.sd = cells5.sd / Nfounder
cells8.mean = cells8.mean / Nfounder
cells8.sd = cells8.sd / Nfounder
# cells10.mean = cells10.mean / Nfounder
# cells10.sd = cells10.sd / Nfounder

# generation vs. time area plot
time <- rep(1:MAX_TIME,times=4) # change to 5 if adding 10hr pulse
cells.mean <- c(as.vector(t(cells1.mean)), as.vector(t(cells3.mean)), as.vector(t(cells5.mean)), as.vector(t(cells8.mean)))
cells.sd <- c(as.vector(t(cells1.sd)), as.vector(t(cells3.sd)), as.vector(t(cells5.sd)), as.vector(t(cells8.sd)))
condition <- rep(c("1hr pulse", "3hr pulse", "5hr pulse", "8hr pulse"), each=MAX_TIME)
data <- data.frame(time, cells.mean, cells.sd, condition)
data$condition <- factor(data$condition , levels=c("1hr pulse", "3hr pulse", "5hr pulse", "8hr pulse") )
# stacked area chart
ggplot(data, aes(x=time, y=cells.mean, color=condition)) + geom_line(linewidth = 1) + xlab("Time (hour)") + ylab("Live cell count (fold change)") +
  # scale_color_manual(values=c("#C00000", "#EE7100", "#FFC000", "#DBDF00", "#8CC547")) + theme_bw() + 
  scale_color_manual(values=c("#C00000", "#EE7100", "#FFC000", "#DBDF00")) + theme_bw() + 
  geom_ribbon(aes(ymin=(cells.mean-cells.sd), ymax=(cells.mean+cells.sd)), linetype=2, alpha=0.1) +
  scale_x_continuous(limits = c(0, MAX_TIME), expand = c(0, 0), breaks = seq(0, MAX_TIME, 24)) +
  scale_y_continuous(limits = c(0, 3.0), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ggplot(data, aes(x=time, y=cells, color=condition)) + geom_line(size = 1) + xlab("Time (hour)") +
#   scale_color_manual(values=c("#C00000", "#C00000", "#C00000", "#C00000")) + theme_bw() + ylim(0, 400) + xlim(0, 120) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



