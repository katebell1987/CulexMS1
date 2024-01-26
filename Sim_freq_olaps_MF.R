#Analysis of frequency and sanity checks of simulated data
#07142023 by M.Fritz

###### SETTING UP WORKSPACE ######

setwd("~/Desktop/Bell_window_sims/")

Sim_full_results <- read.table("./sim_window_files/All_sims.txt", header = F)
Sim_overlaps_results <- read.table("olaps.txt", header = F)
Num_pos_sims <- read.table("num_olaps.txt", header = F)

###### SANITY CHECK SIM RESULTS ######

#How well distributed are the simulations along the chromosomes?
#Avg.num should be the equivalent of the dist windows in empirical dataset
Num_per_Chrom <- tapply(Sim_full_results$V2, Sim_full_results$V1, length)#divide by 1000 >
Num_per_Chrom/1000 

#what is the rightmost (greatest) sim window start per chr?
Max_per_Chrom <- tapply(Sim_full_results$V2, Sim_full_results$V1, max)



#How well distributed are the overlaps along the chromosomes?
Num_per_Chrom <- tapply(Sim_overlaps_results$V2, Sim_overlaps_results$V1, length)#divide by 1000 to get avg/chrom
Num_per_Chrom/1000 #less than half on Chr1 compared to Chr2 & 3 b/c Chr1 has fewer genes

#what is the rightmost (greatest) sim overlap start per sim_window_files chr?
Max_per_Chrom <- tapply(Sim_overlaps_results$V2, Sim_overlaps_results$V1, max)


#line plot dist simulated values per chromosome
Chr1 <- subset(Sim_full_results, V1 == 1) 
Chr2 <- subset(Sim_full_results, V1 == 2)
Chr3 <- subset(Sim_full_results, V1 == 3)

png("Chr1_simWin_dens.png", units = "px", height = 500, width = 800)
hist(Chr1$V2, xlim = c(0, 132876167),  
     breaks = seq(from=0, to=133000000, by=500000), main = "", xlab = "Window Start Position")
dev.off()

png("Chr2_simWin_dens.png", units = "px", height = 500, width = 800)
hist(Chr2$V2, xlim = c(0, 225161840),  
     breaks = seq(from=0, to=226000000, by=500000), main = "", xlab = "Window Start Position")
dev.off()

png("Chr3_simWin_dens.png", units = "px", height = 500, width = 800)
hist(Chr3$V2, xlim = c(0, 199153858), 
     breaks = seq(from=0, to=200000000, by=500000), main = "", xlab = "Window Start Position")
dev.off()


###### FREQUENCY HISTOGRAM RESULTS ######

#table of overlaps
table(Num_pos_sims$V1)

png("sim_freq_plot.png", units = "px", height = 500, width = 800)
par(mar = c(5.1, 5.1, 4.1, 2.1))
hist(Num_pos_sims$V1, main = "", ylab = "Frequency", xlab = "Overlaps per simulation", cex.axis = 1.2, cex.lab = 1.4)

#add vertical lines at x=10 (cutoff) and x = mean overlaps per simulation
abline(v=10, col = "red", lwd = 1.2)
abline(v=mean(Num_pos_sims$V1), col = "blue", lwd = 1.2)

dev.off()


#hist dist olaps per chromosome - gives sense of olap spacing and sanity check for bedtools intersect
#ylim should be greater than for sims plots above, b/c each sim window may have more than one chemosensory gene overlap

Chr1_olaps <- subset(Sim_overlaps_results, V1 == 1)
Chr2_olaps <- subset(Sim_overlaps_results, V1 == 2)
Chr3_olaps <- subset(Sim_overlaps_results, V1 == 3)

png("Chr1_simOlaps_dens.png", units = "px", height = 500, width = 800)
hist(Chr1_olaps$V2, xlim = c(0, 132876167),
     breaks = seq(from=0, to=133000000, by=500000), main = "", xlab = "Window Start Position")
dev.off()

png("Chr2_simOlaps_dens.png", units = "px", height = 500, width = 800)
hist(Chr2_olaps$V2, xlim = c(0, 225161840),
     breaks = seq(from=0, to=226000000, by=500000), main = "", xlab = "Window Start Position")
dev.off()

png("Chr3_simOlaps_dens.png", units = "px", height = 500, width = 800)
hist(Chr3_olaps$V2, xlim = c(0, 199153858),
     breaks = seq(from=0, to=200000000, by=500000), main = "", xlab = "Window Start Position")
dev.off()

