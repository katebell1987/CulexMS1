##### Set WD ####
setwd("~/Dropbox/UMD_Research/MS1_prefGBS/")

#### load data ####
load("RData/entropy_PCA.RData") ## entropy PCA
load("RData/ConRDA_CulexPref.RData") ## RDA
load("RData/barplot.RData") ## barplot

#### Set up layout ####
pdf("~/Dropbox/UMD_Research/MS1_prefGBS/Figures/PCA_RDA_barplot.pdf",height=10,width=15)
nf <- layout( matrix(c(1,2,3,3), nrow=2, byrow=TRUE) )
plot(pcaout$x[,1],pcaout$x[,2],type="n",xlab="",ylab="",main="(a) PCA of Genotype Probabilities",cex.main=2.25)
points(pcaout$x[gprob[,6] == "c",1],pcaout$x[gprob[,6] == "c",2],col=col[gprob[,6] == "c"],pch=symb[gprob[,6] == "c",],bg=bg[gprob[,6] == "c"],cex=2)
points(pcaout$x[gprob[,6] == "e",1],pcaout$x[gprob[,6] == "e",2],col=col[gprob[,6] == "e"],pch=symb[gprob[,6] == "e",],bg=bg[gprob[,6] == "e"],cex=2)
points(pcaout$x[gprob[,6] == "n",1],pcaout$x[gprob[,6] == "n",2],col=col[gprob[,6] == "n"],pch=symb[gprob[,6] == "n",],bg=bg[gprob[,6] == "n"],cex=2)
abline(h=0,lty="dashed",col="lightgrey")
abline(v=0,lty="dashed",col="lightgrey")
mtext("PC1 45.91%",side=1,cex=1.5,line=2.5)
mtext("PC2 10.89%",side=2,cex=1.5,line=2.5)
plot(individual_scores$RDA1,individual_scores$RDA2,xlim=c(-15,15),ylim=c(-15,15),pch=symb[,1],bg=bg,col=col,xlab="",ylab="",main="(b) RDA on Host Choice",cex.main=2.2,cex=2)
text(numeric_centroids$RDA1,numeric_centroids$RDA2,labels=numeric_centroids$Name)
abline(h=0,lty="dashed",col="lightgrey")
abline(v=0,lty="dashed",col="lightgrey")
mtext("RDA1 1.162%",side=1,cex=1.5,line=2.5)
mtext("RDA2 0.934%",side=2,cex=1.5,line=2.5)
barplot(t.means3, col=colors,border=NA, beside=F,space=space,las=1,main="(c) Admixture Proportions K = 3",cex.main=2.25)
mtext("BG1",side=1,at=25,cex=1.5,line=1)
mtext("AG2",side=1,at=70,cex=1.5,line=1)
mtext("AG3",side=1,at=105,cex=1.5,line=1)
mtext("Admixture Proportion q",side=2,line=2.5,cex=1.5)
dev.off()
