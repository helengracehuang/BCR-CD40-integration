library(ggplot2)
library(dplyr)
BCRhigh <- read.table(file="/Users/helenhuang/Desktop/HBCR_Td_5.txt", sep="\t", header = TRUE)
BCRmed <- read.table(file="/Users/helenhuang/Desktop/MBCR_Td_4.txt", sep="\t", header = TRUE)
BCRlow <- read.table(file="/Users/helenhuang/Desktop/LBCR_Td_4.txt", sep="\t", header = TRUE)
Nostim <- read.table(file="/Users/helenhuang/Desktop/nostim_Td_2.txt", sep="\t", header = TRUE)
CD40low <- read.table(file="/Users/helenhuang/Desktop/LCD40_Td_1.txt", sep="\t", header = TRUE)
CD40med <- read.table(file="/Users/helenhuang/Desktop/MCD40_Td_2.txt", sep="\t", header = TRUE)
CD40high <- read.table(file="/Users/helenhuang/Desktop/HCD40_Td_2.txt", sep="\t", header = TRUE)

############# Violin plot of Td distribution ###############
Dose <- rep(c("high", "med", "low", "nostim"), each=125)
mergedTd <- rbind(BCRhigh, BCRmed, BCRlow, Nostim)
# mergedTd <- rbind(CD40high, CD40med, CD40low, Nostim)

mergedTd$Dose <- Dose
mergedTd$Dose <- factor(mergedTd$Dose, levels=c("nostim", "low", "med", "high") )
mergedTd$Td[mergedTd$Td >= 72] <- NA
ggplot(mergedTd, aes(y=Td, x=Dose, fill=Dose)) + 
  geom_violin() + geom_boxplot(width=0.05, fill="white") +
  scale_y_continuous(limits = c(0, 72), expand = c(0.02, 0.02), breaks = seq(0, 120, 24)) +
  scale_fill_brewer(palette="Reds") + theme_classic() + xlab("anti-BCR dose") + ylab("Td (hrs)")
  # scale_fill_brewer(palette="Blues") + theme_classic() + xlab("anti-CD40 dose") + ylab("Td (hrs)")

############# Plot survival curve ###############
BcellSurvival <- list(BCRhigh, BCRmed, BCRlow, Nostim)
# BcellSurvival <- list(CD40high, CD40med, CD40low, Nostim)
Nfounder <- 125
MAX_GEN <- 8
MAX_TIME <- 72
binLength <- 1

AllCells <- list()
for (m in 1:length(BcellSurvival)) {
  cellCount <- nrow(BcellSurvival[[m]])
  cells <- rep(0, MAX_TIME)
  for (i in 1:cellCount) {
    for (j in 0:ceiling(BcellSurvival[[m]]$Td[i])) {
      cells[j+1] = cells[j+1] + 1
    }
  }
  AllCells[[length(AllCells) + 1]] <- cells
}

time <- rep(0:MAX_TIME,times=4)
# cells.num <- c(AllCells[[1]], AllCells[[2]], AllCells[[3]], AllCells[[4]])
cells.num <- c(AllCells[[1]][1:73], AllCells[[2]][1:73], AllCells[[3]][1:73], AllCells[[4]][1:73])
condition <- rep(c("high", "med", "low", "nostim"), each=MAX_TIME+1)
data <- data.frame(time, cells.num, condition)
data$condition <- factor(data$condition , levels=c("nostim", "low", "med", "high") )
# survival curve
ggplot(data, aes(x=time, y=cells.num, color=condition)) + geom_line(size = 1) + xlab("Time (hour)") + ylab("# survived cells") +
  scale_color_brewer(palette="Reds") + theme_bw() +
  # scale_color_brewer(palette="Blues") + theme_bw() +
  scale_x_continuous(limits = c(0, 72), expand = c(0, 0), breaks = seq(0, 72, 12)) +
  scale_y_continuous(limits = c(1, 125), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())






############# Plot 24hr bar plot ###############
BcellSurvival <- list(BCRhigh, BCRmed, BCRlow, Nostim, CD40high, CD40med, CD40low, Nostim)

AllCells <- list()
for (m in 1:length(BcellSurvival)) {
  cellCount <- nrow(BcellSurvival[[m]])
  cells <- rep(0, MAX_TIME)
  for (i in 1:cellCount) {
    for (j in 0:ceiling(BcellSurvival[[m]]$Td[i])) {
      cells[j+1] = cells[j+1] + 1
    }
  }
  AllCells[[length(AllCells) + 1]] <- cells
}

# cells.num <- c(AllCells[[1]], AllCells[[2]], AllCells[[3]], AllCells[[4]])
cells.num <- c(AllCells[[1]][24], AllCells[[2]][24], AllCells[[3]][24], AllCells[[4]][24],
               AllCells[[5]][24], AllCells[[6]][24], AllCells[[7]][24], AllCells[[8]][24]) / Nfounder * 100
condition <- rep(c("BCR-H", "BCR-M", "BCR-L", "nostim", "CD40-H", "CD40-M", "CD40-L", "nastim"), each=1)
data <- data.frame(cells.num, condition)
data$condition <- factor(data$condition , levels=c("nostim", "CD40-L", "CD40-M", "CD40-H", "nastim", "BCR-L", "BCR-M", "BCR-H") )
# survival curve
ggplot(data, aes(x = condition, y = cells.num, fill = condition)) +
  geom_bar(aes(), stat = "identity", position = "dodge") + scale_y_continuous(limits = c(0, 80)) +
  scale_fill_manual(values = c("#DEDEDE", "#BDD7E7", "#6AAFD6", "#2271B5", "#DEDEDE", "#FBAE91", "#FC6A4A", "#CB1C1D")) +
  labs(x = "Stimuli", y = "% Live cells at 24hrs") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
