###Script to identify over-representation of ORs and OBPs in genomic windows associated with Cpip host preference
###Written 07142023 by M. Fritz

 
###### SETTING UP WORKSPACE ######


#loading data
setwd("~/Desktop/Bell_window_sims/")
Data <- read.table("Wstats_full_windowSize.txt", header = T)

#getting significant windows
Wstats_Data <- subset(Data, Wstat >= 3.9687)

#getting means and sd per chromosome for rnorm
Means_per_Chrom <- tapply(Wstats_Data$WindowSize, Wstats_Data$CHR, mean)
Sd_per_Chrom <- tapply(Wstats_Data$WindowSize, Wstats_Data$CHR, sd)
Num_per_Chrom <- tapply(Wstats_Data$WindowSize, Wstats_Data$CHR, length)
Max_per_Chrom <- tapply(Wstats_Data$WindowSize, Wstats_Data$CHR, max)

#looking at histograms of window sizes by chromosome - checking for normality
#pdf(file = "hist_windSize_byChr.pdf")
#par(mfrow = c(3,1))
    #tapply(Wstats_Data$WindowSize, Wstats_Data$CHR,hist)
   
#dev.off()

###### GENERATING WINDOW SIZES AND STARTS PER CHR ######

#window size sim Chr1
#set.seed(1)
#WS1 <- abs(rnorm(Num_per_Chrom[1], Means_per_Chrom[1], Sd_per_Chrom[1]))

#sanity check to match distributions
#pdf(file = "Sim_windSize_Chr1.pdf")
#hist(WS1)
#dev.off()

#distributions look roughly the same, so building datasets

#Chr IDs
Chr1 <- rep(1, n = Num_per_Chrom[1])
Chr2 <- rep(2, n = Num_per_Chrom[2])
Chr3 <- rep(3, n = Num_per_Chrom[3]) 

#window size
WSi1 <- abs(rnorm(Num_per_Chrom[1], Means_per_Chrom[1], Sd_per_Chrom[1]))
WSi2 <- abs(rnorm(Num_per_Chrom[2], Means_per_Chrom[2], Sd_per_Chrom[2]))
WSi3 <- abs(rnorm(Num_per_Chrom[3], Means_per_Chrom[3], Sd_per_Chrom[3]))

#window starts
WSta1 <- sample.int((132876167-(1.5*Max_per_Chrom[1])), Num_per_Chrom[1])
WSta2 <- sample.int((225161840-(1.5*Max_per_Chrom[2])), Num_per_Chrom[2])
WSta3 <- sample.int((199153858-(1.5*Max_per_Chrom[3])), Num_per_Chrom[3])

#window stops
WSto1 <- round((WSta1 + WSi1), digits = 0)
WSto2 <- round((WSta2 + WSi2), digits = 0)
WSto3 <- round((WSta3 + WSi3), digits = 0)

#Gathering sim windows by Chr
Win1 <- data.frame(cbind(Chr1, WSta1, WSto1))
colnames(Win1) <- c('Chr','WindowStart','WindowStop')
Win2 <- data.frame(cbind(Chr2, WSta2, WSto2))
colnames(Win2) <- c('Chr','WindowStart','WindowStop')
Win3 <- data.frame(cbind(Chr3, WSta3, WSto3))
colnames(Win3) <- c('Chr','WindowStart','WindowStop')

#Full sim
Full_sim <- data.frame(rbind(Win1, Win2, Win3))
Uniq_FileName <- as.character(Full_sim[1,2])

write.table(Full_sim, file = paste0(Uniq_FileName, ".sim"), row.names = FALSE, col.names = FALSE, sep = "\t", eol = "\n")

