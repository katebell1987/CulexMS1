#### Set Workign directory ####
setwd("~/Dropbox/UMD_Research/KLB_AAB_GBS/NewAlignment/Kate/")

#### Read in positions ####
dat<-read.csv("all.loci.txt",header=F)
#dat<-read.csv("Kate25i21Positions.txt",header=F)

#### Fix names ####
IDs <- strsplit(as.character(dat[,1]),':') ## remove file ending
IDs2<-do.call(rbind, IDs)
Update<-as.data.frame(as.numeric(IDs2[,2]))

#### Plot of density of markers across chromosome ####
Chr1<-as.data.frame(Update$`as.numeric(IDs2[, 2])`[IDs2[,1]=="JADDTP010000001.1"])
Chr2<-as.data.frame(Update$`as.numeric(IDs2[, 2])`[IDs2[,1]=="JADDTP010000002.1"])
Chr3<-as.data.frame(Update$`as.numeric(IDs2[, 2])`[IDs2[,1]=="JADDTP010000003.1"])

pdf("DensityPlotSNPs.pdf")
par(mfrow=c(3,1))
plot(density(Chr1[,1]),main="Density of SNPs Chr1",xlab="Position")
plot(density(Chr2[,1]),main="Density of SNPs Chr2",xlab="Position")
plot(density(Chr3[,1]),main="Density of SNPs Chr3",xlab="Position")
dev.off()

#### Sliding window plots #####

winstart<-1
winend<-132876167
stepsize<-10000
windowsize<-500000
numberofsteps1<-(winend/windowsize)*(windowsize/stepsize)
wins<-winstart
wine<-winstart+windowsize
out1<-as.data.frame(matrix(NA,nrow=numberofsteps1,ncol=4))
for (i in 1:numberofsteps1){
        out1[i,4]<-length(which(Chr1>= wins & Chr1< wine))
        out1[i,3]<-(wins+wine)/2
        out1[i,2]<-wine
        out1[i,1]<-wins
        wins<-wins+stepsize
        wine<-wine+stepsize
}
colnames(out1)<-c("Start Position","End Position","Midpoint","Number of SNPs")


winstart<-1
winend<-225161840
stepsize<-10000
windowsize<-500000
numberofsteps2<-(winend/windowsize)*(windowsize/stepsize)
wins<-winstart
wine<-winstart+windowsize
out2<-as.data.frame(matrix(NA,nrow=numberofsteps2,ncol=4))
for (i in 1:numberofsteps2){
        out2[i,4]<-length(which(Chr2>= wins & Chr2< wine))
        out2[i,3]<-(wins+wine)/2
        out2[i,2]<-wine
        out2[i,1]<-wins
        wins<-wins+stepsize
        wine<-wine+stepsize
}

winstart<-1
winend<-201550677
stepsize<-10000
windowsize<-500000
numberofsteps3<-(winend/windowsize)*(windowsize/stepsize)
wins<-winstart
wine<-winstart+windowsize
out3<-as.data.frame(matrix(NA,nrow=numberofsteps3,ncol=4))
for (i in 1:numberofsteps3){
        out3[i,4]<-length(which(Chr3>= wins & Chr3< wine))
        out3[i,3]<-(wins+wine)/2
        out3[i,2]<-wine
        out3[i,1]<-wins
        wins<-wins+stepsize
        wine<-wine+stepsize
}

pdf("~/Dropbox/UMD_Research/MS1_prefGBS/Figures/SNPposition_SW.pdf",width=8,height=8)
par(mfrow=c(3,1),mai=c(0.65,0.75,0.1,0.1))
plot(out1[,4],xaxt = "n",xlab="Chromosome 1",ylab="",ylim=c(0,95),type="l",cex.lab=2,cex.axis=1.5)
loc<-seq(1,numberofsteps1,100)
lab<-seq(1,ceiling(numberofsteps1/100),1)
axis(1, at=loc, labels=lab)
plot(out2[,4],xaxt = "n",xlab="Chromosome 2",ylab="",ylim=c(0,95),type="l",cex.lab=2,cex.axis=1.5)
loc<-seq(1,numberofsteps2,100)
lab<-seq(1,ceiling(numberofsteps2/100),1)
axis(1, at=loc, labels=lab)
plot(out3[,4],xaxt = "n",xlab="Chromosome 3",ylab="",ylim=c(0,95),type="l",cex.lab=2,cex.axis=1.5)
loc<-seq(1,numberofsteps3,100)
lab<-seq(1,ceiling(numberofsteps3/100),1)
axis(1, at=loc, labels=lab,cex=2)
par(mfrow=c(1,1),mai=c(0.5,1,0.8,1))
mtext("Number of SNPs per 500KB Window",side=2,cex=1.75,line=3)
dev.off()

##### Histogram of plots per 500kb window ####
pdf("Hist_SnpsPer500kb.pdf")
par(mfrow=c(3,1))
hist(out1[,4],main="Histogram of SNPs per 500kb window",xlab="Number of SNPs Chr1")
text(62,1600,"Mean = 17.31\n Median=16\n Min=0\n Max=73\n Total= 4,634\n Chr Size=13.3MB")

hist(out2[,4],main="Histogram of SNPs per 500kb window",xlab="Number of SNPs Chr2")
text(80,2300,"Mean = 25.06\n Median=24\n Min=0\n Max=95\n Total= 11,301\n Chr Size=22.5MB")

hist(out3[,4],main="Histogram of SNPs per 500kb window",xlab="Number of SNPs Chr3")
text(70,2600,"Mean = 21.27\n Median=20\n Min=0\n Max=76\n Total= 8,594\n Chr Size=20.2MB")

dev.off()

##################################################################

pdf("SNPposition.pdf",width=15,height=10)
par(mfrow=c(3,1))
Chr1<-as.data.frame(Update$`as.numeric(IDs2[, 2])`[IDs2[,1]=="JADDTP010000001.1"])
max(Chr1)
window<-seq(100000,133000000,100000)
Chr1window<-as.data.frame(matrix(NA,nrow=length(window),ncol=2))
for(i in 1:length(window)){
        if(i == 1){Chr1window[i,2]<-length(which(Chr1< window[i]))
                Chr1window[i,1]<-window[i]}
        else {Chr1window[i,2]<-length(which(Chr1>=window[(i-1)] & Chr1< window[i]))
        Chr1window[i,1]<-window[i]}
}
sum(Chr1window$V2) ## Check the numbers match
yl<-max(Chr1window$V2)
Chr1window$V2[Chr1window$V2==0]<-NA
plot(Chr1window$V2,xaxt = "n",xlab="Chromosome 1 MB",ylab="Number of SNPs per 100kb Window",ylim=c(0,yl),type="l")
out<-seq(1,1330,10)
lab<-seq(1,133.8,1)
axis(1, at=out, labels=lab)

Chr2<-as.data.frame(Update$`as.numeric(IDs2[, 2])`[IDs2[,1]=="JADDTP010000002.1"])
max(Chr2)
window2<-seq(100000,226000000,100000)
Chr2window<-as.data.frame(matrix(NA,nrow=length(window2),ncol=2))
for(i in 1:length(window2)){
        if(i == 1){Chr2window[i,2]<-length(which(Chr2< window2[i]))
        Chr2window[i,1]<-window2[i]}
        else {Chr2window[i,2]<-length(which(Chr2>=window2[(i-1)] & Chr2< window2[i]))
        Chr2window[i,1]<-window2[i]}
}
yl<-max(Chr2window$V2)
sum(Chr2window$V2) ## Check the numbers match
Chr2window$V2[Chr2window$V2==0]<-NA
plot(Chr2window$V2,xaxt = "n",xlab="Chromosome 2 MB",ylab="Number of SNPs per 100kb Window",ylim=c(0,yl))
out<-seq(1,2260,10)
lab<-seq(1,226,1)
axis(1, at=out, labels=lab)

Chr3<-as.data.frame(Update$`as.numeric(IDs2[, 2])`[IDs2[,1]=="JADDTP010000003.1"])
max(Chr3)
window3<-seq(100000,202000000,100000)
Chr3window<-as.data.frame(matrix(NA,nrow=length(window3),ncol=2))
for(i in 1:length(window3)){
        if(i == 1){Chr3window[i,2]<-length(which(Chr3< window3[i]))
        Chr3window[i,1]<-window3[i]}
        else {Chr3window[i,2]<-length(which(Chr3>=window3[(i-1)] & Chr3< window3[i]))
        Chr3window[i,1]<-window3[i]}
}
yl<-max(Chr3window$V2)
sum(Chr3window$V2) ## Check the numbers match
Chr3window$V2[Chr3window$V2==0]<-NA
plot(Chr3window$V2,xaxt = "n",xlab="Chromosome 3 MB",ylab="Number of SNPs per 100kb Window",ylim=c(0,yl))
out<-seq(1,2020,10)
lab<-seq(1,202,1)
axis(1, at=out, labels=lab)

dev.off()







###########################################################################################################
dat<-read.csv("Kate20i21Positions.txt",header=F)
IDs <- strsplit(as.character(dat[,1]),':') ## remove file ending
IDs2<-do.call(rbind, IDs)
Update<-as.data.frame(as.numeric(IDs2[,2]))


Chr1<-as.data.frame(Update$`as.numeric(IDs2[, 2])`[IDs2[,1]=="JADDTP010000001.1"])
Chr1$Chr<-rep(0.7,length(Chr1[,1]))
colnames(Chr2)<-c("Position","Chr")
Chr2<-as.data.frame(Update$`as.numeric(IDs2[, 2])`[IDs2[,1]=="JADDTP010000002.1"])
Chr2$Chr<-rep(1.9,length(Chr2[,1]))
colnames(Chr2)<-c("Position","Chr")
Chr3<-as.data.frame(Update$`as.numeric(IDs2[, 2])`[IDs2[,1]=="JADDTP010000003.1"])
Chr3$Chr<-rep(3.1,length(Chr3[,1]))
colnames(Chr3)<-c("Position","Chr")
datplot<-rbind(Chr1,Chr2,Chr3)
size<-as.data.frame(c(132876167,225161840,201550677))
size$Chr<-c("1","2","3")

pdf("test.pdf",width=20,height=5)
bp <- barplot(size[,1], border=NA, col="grey80",horiz = TRUE,xlab="Position")
## location of the vertical segments
with(datplot,
     segments(
       datplot$Position,
       datplot$Chr-0.5,
       datplot$Position,
       datplot$Chr+0.5,
       col="darkgrey",
       lwd=(1e-10000), 
       lend=1
     )
)
#text(x=(-1),y=c(0.7,1.9,3.1),labels = "Chr1")
dev.off()


pdf("Plot_filtered25i21.pdf",width=20,height=5)
bp <- barplot(size[,1], border=NA, col="grey80",horiz = TRUE,xlab="Position")
## location of the vertical segments
with(datplot,
     segments(
       datplot$Position,
       datplot$Chr-0.5,
       datplot$Position,
       datplot$Chr+0.5,
       col="darkgrey",
       lwd=(1e-10000), 
       lend=1
     )
)
#text(x=(-1),y=c(0.7,1.9,3.1),labels = "Chr1")
dev.off()

Chr1<-as.data.frame(dat$V2[dat[,1]=="JADDTP010000001.1"])
Chr1$Chr<-rep(0.7,length(Chr1[,1]))
colnames(Chr1)<-c("Position","Chr")
Chr2<-as.data.frame(dat$V2[dat[,1]=="JADDTP010000002.1"])
Chr2$Chr<-rep(1.9,length(Chr2[,1]))
colnames(Chr2)<-c("Position","Chr")
Chr3<-as.data.frame(dat$V2[dat[,1]=="JADDTP010000003.1"])
Chr3$Chr<-rep(3.1,length(Chr3[,1]))
colnames(Chr3)<-c("Position","Chr")
datplot<-rbind(Chr1,Chr2,Chr3)
size<-as.data.frame(c(132876167,225161840,201550677))
size$Chr<-c("1","2","3")

####


