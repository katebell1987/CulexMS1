################################################################################
#                     Code to run RDA analysis                                 #
#                     for Culex preference data                                #
################################################################################

#### set working directory
setwd("~/Dropbox/UMD_Research/KLB_AAB_GBS/NewAlignment/Kate/entropy/Results/")

### Load libraries ####
#install.packages("vegan")
library(vegan)
require(ggplot2)
require(grid)

#### Read in genotype probabilities for K=2 from entropy #####
g3c0<-read.csv("gprob3c0.txt") # chain 1
g3c1<-read.csv("gprob3c1.txt") # chain 2
combg3<-((g3c0[,-1]+g3c1[,-1])/2) # combine chains


#### Read in IDs with human host and host preference information ####
ids<-read.csv("~/Dropbox/UMD_Research/KLB_AAB_GBS/NewAlignment/Kate/indiv_ids_hostinfo.csv",header=T)
ids$popID<-as.factor(substring(ids$pop,1,1)) # make sure populations are correct

#### Run conditional RDA ####
conRda<-rda(combg3~ ids$Response +Condition(ids$popID)) # Run model

#### Run an anova ####
set.seed(966124) # set seed
anova(conRda,permutations = 9999,parallel = 6,by="margin") # Run anova, need to set seed and permutations 
# Permutation test for rda under NA model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 9999
# Model: rda(formula = combg3 ~ ids$Response + Condition(ids$popID))
# Df Variance      F Pr(>F)  
# ids$Response   2     95.0 1.1774  0.022 *
# Residual     110   4436.5                
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
set.seed(NULL)

#### Explore output from model ####
RsquareAdj(conRda)$r.squared ## Unadjusted R2
RsquareAdj(conRda)$adj.r.squared ## Adjusted R2
summary(eigenvals(conRda))

#### Make GG Plots of individuals for RDA1 vs RDA2 ####
scor <- scores(conRda, display=c("sp", "cn", "bp"), scaling=3,choices = c(1:2)) # get scores
numeric_centroids <- data.frame(scor$centroids) # get centroids for host preference
numeric_centroids$Name<-c("C","H","S") # add names to be displayed on plot
individual_scores <- data.frame(scores(conRda,choices = c(1:2))$sites) # get scores for each individual
individual_scores$individuals<-ids$ind # add id names
bg1<-c((rep("#FFD70040",48)),(rep("#00bfff40",44)),(rep("#228B2240",23))) # set colors
col1<-c((rep("gold",48)),(rep("deepskyblue",44)),(rep("forestgreen",23)))
symb1<-as.data.frame(matrix(nrow=115,ncol=1,NA)) # data frame for symbols
for (i in 1:length(ids$Response)){ # for loop to assign pch values for host choice type
  if (ids$Response[i]=="h") symb1[i,1]<-24 #triangle
  else if (ids$Response[i]=="c") symb1[i,1]<-21 #circle
  else symb1[i,1]<-22 #square
}

pdf("~/Dropbox/UMD_Research/MS1_prefGBS/Figures/RDA_1vs2_conditionalRDA.pdf",width=10, height=8) # make the plot, write to a pdf
(RDA_plot <- ggplot(numeric_centroids, aes(x = RDA1, y = RDA2))+labs(x="RDA1 1.162%",y="RDA2 0.934%",color="Population")+geom_vline(xintercept=c(0, 0), linetype = "dashed",col="grey")+geom_hline(yintercept=c(0,0), linetype = "dashed",col="grey")+
    geom_point(data = individual_scores, aes(x = RDA1, y = RDA2),
               size = 3,
               colour=col,
               pch=symb[,1],
               bg = bg) +
    geom_text(aes(label = Name),size=3)  +
    coord_cartesian(x = c(-15,15), y = c(-15,15))+ theme_bw()) + theme( panel.grid.major = element_blank(),
                                                                        panel.grid.minor = element_blank())
dev.off()

#### Write data to .RData file ####
save.image(file="~/Dropbox/UMD_Research/MS1_prefGBS/RData/ConRDA_CulexPref.RData")

############################# end of script ####################################