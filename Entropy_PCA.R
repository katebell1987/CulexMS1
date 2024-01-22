#### set working directory ####
setwd("~/Dropbox/UMD_Research/KLB_AAB_GBS/NewAlignment/Kate/entropy/Results/")

#### Read in Data ####
g3c0<-read.csv("gprob3c0.txt")
g3c1<-read.csv("gprob3c1.txt")

combg3<-((g3c0[,-1]+g3c1[,-1])/2)


#### Read in IDs with human host and host preference information ####
ids<-read.csv("../../indiv_ids_hostinfo.csv",header=T)
ids$popID<-as.factor(substring(ids$pop,1,1)) ## Populations

## Set Colors and Symbols
bg<-c((rep("#FFD70040",48)),(rep("#00bfff40",44)),(rep("#228B2240",23))) # set colors
col<-c((rep("gold",48)),(rep("deepskyblue",44)),(rep("forestgreen",23)))
symb<-as.data.frame(matrix(nrow=115,ncol=1,NA)) # data frame for symbols
for (i in 1:length(ids$Response)){ # for loop to assign pch values for host choice type
  if (ids$Response[i]=="h") symb[i,1]<-24 #triangle
  else if (ids$Response[i]=="c") symb[i,1]<-21 #circle
  else symb[i,1]<-22 #square
}

gprob<-cbind(ids,combg3)
gprob[1:5,1:10]
gprobMat<-gprob[,-c(1:6)]
gprobMat[1:5,1:5]

#### PCA ####
pcaout<-prcomp(gprobMat,center=T,scale=F)
summary(pcaout) ## Proportion of Variance  0.4591  0.1089  0.0208  0.01529

save.image("~/Dropbox/UMD_Research/MS1_prefGBS/RData/entropy_PCA.RData")
# Plot symbols in two columns, shifted to the left by 3 and 1 respectively

pdf("~/Dropbox/UMD_Research/MS1_prefGBS/Figures/PCA_gprobsK3_1vs2.pdf")
plot(pcaout$x[,1],pcaout$x[,2],type="n",xlab="PC1 45.91%",ylab="PC2 10.89%")
points(pcaout$x[gprob[,6] == "c",1],pcaout$x[gprob[,6] == "c",2],col=col[gprob[,6] == "c"],pch=symb[gprob[,6] == "c",],bg=bg[gprob[,6] == "c"],cex=1)
points(pcaout$x[gprob[,6] == "e",1],pcaout$x[gprob[,6] == "e",2],col=col[gprob[,6] == "e"],pch=symb[gprob[,6] == "e",],bg=bg[gprob[,6] == "e"],cex=1)
points(pcaout$x[gprob[,6] == "n",1],pcaout$x[gprob[,6] == "n",2],col=col[gprob[,6] == "n"],pch=symb[gprob[,6] == "n",],bg=bg[gprob[,6] == "n"],cex=1)
dev.off()
legend('bottomleft', legend=c("BG1 S","BG1 A","BG1 H","AG2 S","AG2 A","AG2 H","AG3 S","AG3 A","AG3 H"),pch=c(22,21,24,22,21,24,22,21,24),col=c("gold","gold","gold","skyblue","skyblue","skyblue","forestgreen","forestgreen","forestgreen"),bg=,cex=1,pt.cex=1.2,ncol=3)

pdf("PCA_gprobsK3_1vs3.pdf")
plot(pcaout$x[,1],pcaout$x[,3],type="n",xlab="PC1 45.91%",ylab="PC3 1.9%")
points(pcaout$x[gprob[,6] == "c",1],pcaout$x[gprob[,6] == "c",3],col=col[gprob[,6] == "c"],pch=symb[gprob[,6] == "c",],cex=1)
points(pcaout$x[gprob[,6] == "e",1],pcaout$x[gprob[,6] == "e",3],col=col[gprob[,6] == "e"],pch=symb[gprob[,6] == "e",],cex=1)
points(pcaout$x[gprob[,6] == "n",1],pcaout$x[gprob[,6] == "n",3],col=col[gprob[,6] == "n"],pch=symb[gprob[,6] == "n",],cex=1)
legend('bottomleft', legend=c("BG1 S","BG1 A","BG1 H","AG2 S","AG2 A","AG2 H","AG3 S","AG3 A","AG3 H"),pch=c(22,21,24,22,21,24,22,21,24),col=c("gold","gold","gold","skyblue","skyblue","skyblue","forestgreen","forestgreen","forestgreen"),cex=1,pt.cex=1.2,ncol=3)
dev.off()

pdf("PCA_gprobsK3_2vs3.pdf")
plot(pcaout$x[,2],pcaout$x[,3],type="n",xlab="PC2 10.89%",ylab="PC3 1.9%")
points(pcaout$x[gprob[,6] == "c",2],pcaout$x[gprob[,6] == "c",3],col=col[gprob[,6] == "c"],pch=symb[gprob[,6] == "c",],cex=1)
points(pcaout$x[gprob[,6] == "e",2],pcaout$x[gprob[,6] == "e",3],col=col[gprob[,6] == "e"],pch=symb[gprob[,6] == "e",],cex=1)
points(pcaout$x[gprob[,6] == "n",2],pcaout$x[gprob[,6] == "n",3],col=col[gprob[,6] == "n"],pch=symb[gprob[,6] == "n",],cex=1)
legend('bottomleft', legend=c("BG1 S","BG1 A","BG1 H","AG2 S","AG2 A","AG2 H","AG3 S","AG3 A","AG3 H"),pch=c(22,21,24,22,21,24,22,21,24),col=c("gold","gold","gold","skyblue","skyblue","skyblue","forestgreen","forestgreen","forestgreen"),cex=1,pt.cex=1.2,ncol=3)
dev.off()


