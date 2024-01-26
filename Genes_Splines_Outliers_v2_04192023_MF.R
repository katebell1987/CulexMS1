#Fritz Computing notes for overlapping genomic windows HostPref RDA manuscript 
#19Apr2023
#MF

###### PREPPING WORKSPACE #######
x <- c("dplyr", "GenWin", "ggpubr", "ggplot2") 
lapply(x, FUN = function(X) {do.call("library", list(X))}) #loading libraries

setwd("~/Downloads/Splines_ForKate/Splines_ForKate")

#loading file with outliers according to q-value
outliers <- read.csv("qvalue_outliers.csv", header = T)

#loading full dataset - pvalue for each locus.
pvals <- read.csv("full_pqvalues.csv", header = T)

#loading liftoff gff and manually curated chemosensory gene files to get mRNA coordinates
gff <- read.table("2023_03_12_sortMergedManualRefSeqAnnotation.gtf", sep="\t", header=FALSE, comment.char="#",
                  na.strings=".", stringsAsFactors=FALSE,
                  quote="", fill=FALSE)


####### DATA CHECK ######

#How are outliers distributed in the genome - which of the 3 chromosomes?
table(outliers$chromosome) #interesting that most are found on Chromosome 3.

#How many of these are clustered within a single tag a chromosome?
Diffs <- tapply(outliers$position, outliers$chromosome, diff)#gets diff between marker position along each Chr
Diffs #~7 appear < 100bp from one another on chromosomes 2 & 3, meaning on a single tag.   

#sanity check - making sure Diffs worked....
outliers$position[1] - outliers$position[2]#matches first obs for Diffs above

#sanity check - gff
gff %>% count(V3)#getting num transcripts = 28448


###### USE SMOOTHED SPLINES TO GENERATE WINDOWS ######

pvals$transP <- -log10(pvals$p.values) #transform pvals for plot

#subset chromosomes
CHR1 <- subset(pvals, Chromosome == "JADDTP010000001.1")
CHR2 <- subset(pvals, Chromosome == "JADDTP010000002.1")
CHR3 <- subset(pvals, Chromosome == "JADDTP010000003.1")

#generate splines
CHR1_SA <- splineAnalyze(CHR1$transP, CHR1$position, smoothness = 100, plotRaw = TRUE, plotWindows = TRUE, method = 4)
CHR2_SA <- splineAnalyze(CHR2$transP, CHR2$position, smoothness = 100, plotRaw = TRUE, plotWindows = TRUE, method = 4)
CHR3_SA <- splineAnalyze(CHR3$transP, CHR3$position, smoothness = 100, plotRaw = TRUE, plotWindows = TRUE, method = 4)

#Regenerate windowed dataset using 95% quantile as cutoff for Wstat (Beissinger et al. 2015 used 99%)
Wstats_CHR1 <- data.frame(CHR1_SA$windowData)
Wstats_CHR1$CHR <- rep("JADDTP010000001.1", times=nrow(Wstats_CHR1))

Wstats_CHR2 <- data.frame(CHR2_SA$windowData)
Wstats_CHR2$CHR <- rep("JADDTP010000002.1", times=nrow(Wstats_CHR2))

Wstats_CHR3 <- data.frame(CHR3_SA$windowData)
Wstats_CHR3$CHR <- rep("JADDTP01000003.1", times=nrow(Wstats_CHR3))

#generating full dataset to get genome-wide 95% quantile
Wstats_full <- rbind(Wstats_CHR1, Wstats_CHR2, Wstats_CHR3)
Lim_Wstat <- quantile(Wstats_full$Wstat, prob=c(.95)) #3.97

#top Wstats by Chromosome
Top_Wstats_CHR1 <- subset(Wstats_CHR1, Wstat >= Lim_Wstat)
Top_Wstats_CHR2 <- subset(Wstats_CHR2, Wstat >= Lim_Wstat)
Top_Wstats_CHR3 <- subset(Wstats_CHR3, Wstat >= Lim_Wstat)


###### CHECKING CORRELATIONS BETWEEN WINDOW SIZE, SNP COUNT AND WSTAT ######

#getting window size
Wstats_full$WindowSize <- Wstats_full$WindowStop - Wstats_full$WindowStart


#plotting corr of SNP count by window size with spearman 
SNPcountBySize <- data.frame(Wstats_full$WindowSize, Wstats_full$SNPcount)

colnames(SNPcountBySize)[1] ="WindowSize"
colnames(SNPcountBySize)[2] ="SNPcount"

head(SNPcountBySize)

ggplot(SNPcountBySize, aes( x=WindowSize, y=SNPcount ))+
  geom_point(pch=19) + stat_cor(method = "kendall", label.x = -5, label.y = 90) + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#plotting distributions of Wstats by SNP counts for Sig and NonSig windows 
SNPcountBySize$Wstat <- Wstats_full$Wstat

Sig <- subset(SNPcountBySize, Wstat >=3.97)
NonSig <- subset(SNPcountBySize, Wstat < 3.97)

df<-data.frame(X=c(Sig$SNPcount,NonSig$SNPcount),Grp=rep(c("Sig","NonSig"),times=c(nrow(Sig),nrow(NonSig))))

boxplot(X~Grp,data=df, notch = T, ylab = "SNP count", xlab="")


###### CHECKING WSTAT OVERLAP WITH QVALS ######

overlaps_CHR1 <- data.frame()

for (i in 1:nrow(Top_Wstats_CHR1)){
  int1 <- as.numeric(Top_Wstats_CHR1$WindowStart[[i]])
  int2 <- as.numeric(Top_Wstats_CHR1$WindowStop[[i]])
  print(int1)
  print(int2)
  ov <- subset(outliers[c(1:8),], position %in% int1:int2)
  overlaps_CHR1 <- rbind(overlaps_CHR1, ov)
} 

overlaps_CHR2 <- data.frame()

for (i in 1:nrow(Top_Wstats_CHR2)){
  int1 <- as.numeric(Top_Wstats_CHR2$WindowStart[[i]])
  int2 <- as.numeric(Top_Wstats_CHR2$WindowStop[[i]])
  print(int1)
  print(int2)
  ov <- subset(outliers[c(9:32),], position %in% int1:int2)
  overlaps_CHR2 <- rbind(overlaps_CHR2, ov)
} 

overlaps_CHR3 <- data.frame()

for (i in 1:nrow(Top_Wstats_CHR3)){
  int1 <- as.numeric(Top_Wstats_CHR3$WindowStart[[i]])
  int2 <- as.numeric(Top_Wstats_CHR3$WindowStop[[i]])
  print(int1)
  print(int2)
  ov <- subset(outliers[c(33:80),], position %in% int1:int2)
  overlaps_CHR3 <- rbind(overlaps_CHR3, ov)
} 

#proportion of outliers in Spline Windows
(sum(nrow(overlaps_CHR1),nrow(overlaps_CHR2),nrow(overlaps_CHR3)))/nrow(outliers) #0.81 woohoo!


###### PREPPING FOR BEDTOOLS ######

#Splines bedfile
Out_forJoin <- rbind(Top_Wstats_CHR1, Top_Wstats_CHR2, Top_Wstats_CHR3)
Out_forJoin <- Out_forJoin[,c(6,1,2)]
names(Out_forJoin) <- c('chrom', 'chromStart','chromEnd') #prepping for bedtools
Out_forJoin$chrom <- as.numeric(as.factor(Out_forJoin$chrom))

write.table(Out_forJoin, file = "RDA_SPLINES_Out_forJoin.bed", col.names = F, row.names = F, sep = "\t", eol = "\n")

#gff bedfile
gff_genes_withMito <- subset(gff, V3 == "transcript") #getting mRNA coordinates only
gff_genes <- subset(gff_genes_withMito, V1 != "NC_014574.1")
gff_forJoin <- gff_genes[,c(1,4,5)]
str(gff_forJoin)
names(gff_forJoin)<- c('chrom', 'chromStart','chromEnd') #prepping for bedtools
gff_forJoin$chrom <- as.numeric(as.factor(gff_forJoin$chrom))

write.table(gff_forJoin, file = "gff_forJoin.bed", col.names = F, row.names = F, sep = "\t", eol = "\n")


###### BEDTOOLS V.2.27.1 COMMANDS FOR OVERLAPS ######

#bedtools intersect -a RDA_SPLINES_Out_forJoin.bed -b gff_forJoin.bed -wo | wc -l #gives the number of overlaps as 1911
#bedtools intersect -a RDA_SPLINES_Out_forJoin.bed -b gff_forJoin.bed -wo > olaps.txt


###### OVERLAPS FROM BEDTOOLS OUTPUT ######

#reading in and prepping bedtools output
olaps <- read.table("olaps.txt", header = F)
olaps$gene_coords <- paste0(olaps$V4,"_",olaps$V5)


#getting mRNA information
gff_genes$gene_coords <- paste0(as.numeric(as.factor(gff_genes$V1)),"_",gff_genes$V4)
str(gff_genes$gene_coords)
gff_genes$gene_coords <- as.character(gff_genes$gene_coords)

Uniq_gene_coords <- data.frame(unique(olaps$gene_coords))

names(Uniq_gene_coords)[1] <- "gene_coords"
str(Uniq_gene_coords$gene_coords)
Uniq_gene_coords$gene_coords <- as.character(Uniq_gene_coords$gene_coords)

table(olaps$V1)

###just pulling out ORs and OBPs
all_gene_olaps <- merge(Uniq_gene_coords, gff_genes, by = "gene_coords")

#based on row names for ORs and OBPs, subsetting dataframe.  36 total genes
chemosens_gene_olaps <- all_gene_olaps[c(318,362,828,829,830,831,832,833,979,980,981,982,983,984,985,986,987,1017,1134,1136,1227,1228,1230,1231,1232,1233,1234,1235,1236,1237,1560,1568,1780,1788,1791,1803),]
write.table(chemosens_gene_olaps, file = "chemosens_gene_olaps.txt", sep = "\t", row.names = F)

#how many genes overlap per chromosome?
table(chemosens_gene_olaps$V1)

