#This script checks that the genotypes of human and chick seekers are very similar for Kate and Anna.
#01152021
#MF

###### PREP WORKSPACE ######
x <- c("gridExtra", "ggplot2", "reshape2") 
lapply(x, FUN = function(X) {do.call("library", list(X))})

setwd("~/Downloads/KateBell_Hered_MS_2022/RDA_FilesForMegan")


####### LOADING DATA #######
genos_qvals <- read.csv("genotypes_qvalue_outliers.csv", header = T)


####### WIDE TO LONG REFORMAT #######
long_dat <- reshape(genos_qvals, direction="long", 
        varying=c(colnames(genos_qvals[,c(7:86)])), 
        timevar = "locus",
        v.names="genoProb",
        idvar='ind')


#sanity check on reformat
ind1 <- subset(long_dat, ind =="cchh_046")
head(genos_qvals)
head(ind1)

head(genos_qvals[,c(81:86)])
tail(ind1)


####### SUBSET THE DATA FOR BOX PLOTS GENOPROBS #######
Rem_Nor_chix <- subset(long_dat, popID != "n" & Response == "c")
Rem_Nor_chix$locus <- as.factor(Rem_Nor_chix$locus)

Rem_Nor_swit <- subset(long_dat, popID != "n" & Response == "s")
Rem_Nor_swit$locus <- as.factor(Rem_Nor_swit$locus)

Rem_Nor_hum <- subset(long_dat, popID != "n" & Response == "h")
Rem_Nor_hum$locus <- as.factor(Rem_Nor_hum$locus)



####### GENERATE BOX PLOTS OF GENOPROBS BY HOST/CHOICE #######

png("GenoProbs_byHostandLocus_hum.png", width = 1800, height = 800, units = "px") #observations per locus
ggplot(Rem_Nor_hum, aes(x=locus, y=genoProb, fill = HumanHost)) + 
  geom_boxplot(notch = FALSE)+
  #geom_jitter(position=position_dodge(1))+
  labs(title="Human Responders",
       x ="Genotype Probabilities", y = "Top Variance Contributing Loci") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey30"),
        panel.border = element_rect(colour = "grey30", fill=NA, size=1.2), strip.text = element_text(size=15),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.x=element_text(angle=45, hjust=1, size = 11),
        axis.text.y=element_text(size = 11))

dev.off()

png("GenoProbs_byHostandLocus_swit.png", width = 1800, height = 800, units = "px")
ggplot(Rem_Nor_swit, aes(x=locus, y=genoProb, fill = HumanHost)) + 
  geom_boxplot(notch = FALSE)+
  #geom_jitter(position=position_dodge(1))+
  labs(title="Switchers",
       x ="Genotype Probabilities", y = "Top Variance Contributing Loci") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey30"),
        panel.border = element_rect(colour = "grey30", fill=NA, size=1.2), strip.text = element_text(size=15),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.x=element_text(angle=45, hjust=1, size = 11),
        axis.text.y=element_text(size = 11))

dev.off()

png("GenoProbs_byHostandLocus_chix.png", width = 1800, height = 800, units = "px")
ggplot(Rem_Nor_chix, aes(x=locus, y=genoProb, fill = HumanHost)) + 
  geom_boxplot(notch = FALSE)+
  #geom_jitter(position=position_dodge(1))+
  labs(title="Chick Responders",
       x ="Genotype Probabilities", y = "Top Variance Contributing Loci") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey30"),
        panel.border = element_rect(colour = "grey30", fill=NA, size=1.2), strip.text = element_text(size=15),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.x=element_text(angle=45, hjust=1, size = 11),
        axis.text.y=element_text(size = 11))

dev.off()


####### SUBSET THE DATA FOR MEANS PLOT #######
Rem_Nor_K <- subset(genos_qvals, popID != "n" & HumanHost =="K")
summary(Rem_Nor_K[,c(1:6)]) #sanity check

Rem_Nor_A <- subset(genos_qvals, popID != "n" & HumanHost =="A")
summary(Rem_Nor_A[,c(1:6)]) #sanity check


####### LOOP THRU LOCI - MEAN GENO PROBS PER HOST/CHOICE #######
Kate_means <- data.frame()

for (i in colnames(Rem_Nor_K[,c(7:max(ncol(Rem_Nor_K)))])){
  Kate_means <- rbind(Kate_means, tapply(Rem_Nor_K[[i]], Rem_Nor_K$Response, mean))
}

colnames(Kate_means) <- c("genoProb_chix_K", "genoProb_hum_K", "genoProb_swit_K")



Anna_means <- data.frame()

for (i in colnames(Rem_Nor_A[,c(7:max(ncol(Rem_Nor_A)))])){
  Anna_means <- rbind(Anna_means, tapply(Rem_Nor_A[[i]], Rem_Nor_A$Response, mean))
}

colnames(Anna_means) <- c("genoProb_chix_A", "genoProb_hum_A", "genoProb_swit_A")

#dataframe for corr test and plot of mean genotype probs
full_df <- cbind(Kate_means, Anna_means)


####### PLOTTING RELATIONSHIP GENO PROBS BY HOST/CHOICE #######
png("GenoProb_Corrs.png", width = 1800, height = 800, units = "px") #the means

p1 <- ggplot(full_df, aes(x=genoProb_hum_K, y=genoProb_hum_A)) + 
  labs(title="Human Responders",
       x ="Genotype Probabilities - Hum 1", y = "Genotype Probabilities - Hum 2") +
  geom_point(shape=16, size = 3, color="black") +
  geom_smooth(method=lm) + 
  theme_classic() + 
  theme(plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5),
        axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))

p2 <- ggplot(full_df, aes(x=genoProb_swit_K, y=genoProb_swit_A)) + 
  labs(title="Switchers",
       x ="Genotype Probabilities - Hum 1", y = "Genotype Probabilities - Hum 2") +
  geom_point(shape=16, size = 3, color="black") +
  geom_smooth(method=lm) + 
  theme_classic() + 
  theme(plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5),
        axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))

p3 <- ggplot(full_df, aes(x=genoProb_chix_K, y=genoProb_chix_A)) + 
  labs(title="Avian Responders",
       x ="Genotype Probabilities - Hum 1", y = "Genotype Probabilities - Hum 2") +
  geom_point(shape=16, size = 3, color="black") +
  geom_smooth(method=lm) + 
  theme_classic() + 
  theme(plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5),
        axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))

grid.arrange(p1, p2, p3, nrow = 1)

dev.off()


###### CORRELATION COEFFICIENTS PER HUMAN/CHOICE #######

####using pearson assumes bivariate normal dist.  Since data are bounded between 0-2, 
####probably not normally distributed.  Using Kendall's tau, which is non-parametric

#Geno prob corrs for Human Resp
cor.test(full_df$genoProb_hum_K, full_df$genoProb_hum_A, method = "kendall", conf.level = 0.95)

#Geno prob corrs for Switchers
cor.test(full_df$genoProb_swit_K, full_df$genoProb_swit_A, method = "kendall", conf.level = 0.95)

#Geno prob corrs for Avian Resp
cor.test(full_df$genoProb_chix_K, full_df$genoProb_chix_A, method = "kendall", conf.level = 0.95)

