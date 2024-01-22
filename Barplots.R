### Set WD
setwd("~/Dropbox/UMD_Research/KLB_AAB_GBS/NewAlignment/Kate/entropy/Results/")
### Structure bar plots of q
names<-read.csv("../../indiv_ids.txt",header=TRUE)
nind<-115

###################### Set up Info for Plot ####
space<-c(rep(0,48),1,rep(0,43),1,rep(0,22))

legend_names<-c("BG1:","AG2","AG3")

####################### Read in Data #########################

############################# k2 #############################
q2c0<-read.csv("q2c0.txt",header=T)
q2c1<-read.csv("q2c1.txt",header=T)

## Combine chains and get means ##
q2meanCombined<-cbind(q2c0$mean,q2c1$mean) 
q2means<-apply(q2meanCombined,1,FUN=mean) # 1 = rows

## k2 ##
q2meanMat<-cbind(q2means[1:nind],q2means[(1 + nind):(nind * 2)])
rowSums<-apply(q2meanMat,1,FUN=sum)
tempMat<-cbind(q2meanMat,rowSums)
q2meanFixed<-matrix(NA,nrow = dim(q2meanMat)[1],ncol = dim(q2meanMat)[2])
for (i in 1:dim(tempMat)[1]) {
  for (j in 1:2){
    q2meanFixed[i,j]<-tempMat[i,j]/tempMat[i,3]
  }
}
q2meanMat[1:5,1:2]
q2meanFixed[1:5,1:2]
# check:
apply(q2meanFixed,1,FUN=sum)

cluster21<-cbind(names,q2meanFixed[,1])
cluster22<-cbind(names,q2meanFixed[,2])
means2<-data.frame(cbind(cluster21[,4],cluster22[,4]))
t.means2<-t(means2)


################################ k3 #############################
q3c0<-read.csv("q3c0.txt",header=T)
q3c1<-read.csv("q3c1.txt",header=T)

## Combine chains and get means ##
q3meanCombined<-cbind(q3c0$mean,q3c1$mean) 
q3means<-apply(q3meanCombined,1,FUN=mean) # 1 = rows

## k3 ##
q3meanMat<-cbind(q3means[1:nind],q3means[(1 + nind):(nind * 2)],q3means[(1+nind*2):(nind*3)])
rowSums<-apply(q3meanMat,1,FUN=sum)
tempMat<-cbind(q3meanMat,rowSums)
q3meanFixed<-matrix(NA,nrow = dim(q3meanMat)[1],ncol = dim(q3meanMat)[2])
for (i in 1:dim(tempMat)[1]) {
  for (j in 1:3){
    q3meanFixed[i,j]<-tempMat[i,j]/tempMat[i,4]
  }
}
q3meanMat[1:5,1:3]
q3meanFixed[1:5,1:3]
# check:
apply(q3meanFixed,1,FUN=sum)

cluster31<-cbind(names,q3meanFixed[,1])
cluster32<-cbind(names,q3meanFixed[,2])
cluster33<-cbind(names,q3meanFixed[,3])
means3<-data.frame(cbind(cluster31[,4],cluster32[,4],cluster33[,4]))
t.means3<-t(means3)

################################### k4 ###########################
q4c0<-read.csv("q4c0.txt",header=T)
q4c1<-read.csv("q4c1.txt",header=T)

## Combine chains and get means ##
q4meanCombined<-cbind(q4c0$mean,q4c1$mean) 
q4means<-apply(q4meanCombined,1,FUN=mean) # 1 = rows

## k4 ##
q4meanMat<-cbind(q4means[1:nind],q4means[(1 + nind):(nind * 2)],q4means[(1+nind*2):(nind*3)],q4means[(1+nind*3):(nind*4)])
rowSums<-apply(q4meanMat,1,FUN=sum)
tempMat<-cbind(q4meanMat,rowSums)
q4meanFixed<-matrix(NA,nrow = dim(q4meanMat)[1],ncol = dim(q4meanMat)[2])
for (i in 1:dim(tempMat)[1]) {
  for (j in 1:4){
    q4meanFixed[i,j]<-tempMat[i,j]/tempMat[i,5]
  }
}
q4meanMat[1:5,]
q4meanFixed[1:5,]
# check:
apply(q4meanFixed,1,FUN=sum)

cluster41<-cbind(names,q4meanFixed[,1])
cluster42<-cbind(names,q4meanFixed[,2])
cluster43<-cbind(names,q4meanFixed[,3])
cluster44<-cbind(names,q4meanFixed[,4])
means4<-data.frame(cbind(cluster41[,4],cluster42[,4],cluster43[,4],cluster44[,4]))
t.means4<-t(means4)

############################### k5 ############################
q5c0<-read.csv("q5c0.txt",header=T)
q5c1<-read.csv("q5c1.txt",header=T)

## Combine chains and get means ##
q5meanCombined<-cbind(q5c0$mean,q5c1$mean) 
q5means<-apply(q5meanCombined,1,FUN=mean) # 1 = rows

## k5 ##
q5meanMat<-cbind(q5means[1:nind],q5means[(1 + nind):(nind * 2)],q5means[(1+nind*2):(nind*3)],q5means[(1+nind*3):(nind*4)],q5means[(1+nind*4):(nind*5)])
rowSums<-apply(q5meanMat,1,FUN=sum)
tempMat<-cbind(q5meanMat,rowSums)
q5meanFixed<-matrix(NA,nrow = dim(q5meanMat)[1],ncol = dim(q5meanMat)[2])
for (i in 1:dim(tempMat)[1]) {
  for (j in 1:5){
    q5meanFixed[i,j]<-tempMat[i,j]/tempMat[i,6]
  }
}
q5meanMat[1:5,]
q5meanFixed[1:5,]
# check:
apply(q5meanFixed,1,FUN=sum)

cluster51<-cbind(names,q5meanFixed[,1])
cluster52<-cbind(names,q5meanFixed[,2])
cluster53<-cbind(names,q5meanFixed[,3])
cluster54<-cbind(names,q5meanFixed[,4])
cluster55<-cbind(names,q5meanFixed[,5])
means5<-data.frame(cbind(cluster51[,4],cluster52[,4],cluster53[,4],cluster54[,4],cluster55[,4]))
t.means5<-t(means5)

################################# k6 ########################
q6c0<-read.csv("q6c0.txt",header=T)
q6c1<-read.csv("q6c1.txt",header=T)

## Combine chains and get means ##
q6meanCombined<-cbind(q6c0$mean,q6c1$mean) 
q6means<-apply(q6meanCombined,1,FUN=mean) # 1 = rows

## k6 ##
q6meanMat<-cbind(q6means[1:nind],q6means[(1 + nind):(nind * 2)],q6means[(1+nind*2):(nind*3)],q6means[(1+nind*3):(nind*4)],q6means[(1+nind*4):(nind*5)],q6means[(1+nind*5):(nind*6)])
rowSums<-apply(q6meanMat,1,FUN=sum)
tempMat<-cbind(q6meanMat,rowSums)
q6meanFixed<-matrix(NA,nrow = dim(q6meanMat)[1],ncol = dim(q6meanMat)[2])
for (i in 1:dim(tempMat)[1]) {
  for (j in 1:6){
    q6meanFixed[i,j]<-tempMat[i,j]/tempMat[i,7]
  }
}
q6meanMat[1:5,]
q6meanFixed[1:5,]
# check:
apply(q6meanFixed,1,FUN=sum)

cluster61<-cbind(names,q6meanFixed[,1])
cluster62<-cbind(names,q6meanFixed[,2])
cluster63<-cbind(names,q6meanFixed[,3])
cluster64<-cbind(names,q6meanFixed[,4])
cluster65<-cbind(names,q6meanFixed[,5])
cluster66<-cbind(names,q6meanFixed[,6])
means6<-data.frame(cbind(cluster61[,4],cluster62[,4],cluster63[,4],cluster64[,4],cluster65[,4],cluster66[,4]))
t.means6<-t(means6)

save.image("~/Dropbox/UMD_Research/MS1_prefGBS/RData/barplot.RData")

####################### Plot 2 through 5 ####################
pdf("K2_K6_30xi21.pdf",height=6,width=9)
grid <- matrix(c(1:6),nrow=6,ncol=1,byrow=TRUE)
layout(grid)
par(mar=c(0,0,0,0))

## Plot k2 ##
par(fig=c(0.075,0.995,0.77,0.95),new=T) ## 3rd = bottom, 4th = top 0.02 apart
par(mar=c(0,1.5,0.5,1.5))
#colors= c( "#313695","#a50026") ##9970ab"
colors= c( "gold","forestgreen") ##9970ab"
barplot(t.means2, col=colors,border=NA, beside=F,space=space,las=1)

## Plot k3 ##
par(fig=c(0.075,0.995,0.59,0.76),new=T)
par(mar=c(0,1.5,0.5,1.5))
#colors= c( "#a50026","#fdae61","#313695")
colors= c("forestgreen","skyblue","gold") ## Red = northfield, blue = eva, orange = cal1
barplot(t.means3, col=colors,border=NA, beside=F,space=space,las=1)


## Plot k4 ##
par(fig=c(0.075,0.995,0.41,0.58),new=T)
par(mar=c(0,1.5,0.5,1.5))
#colors= c("#313695","#a50026","#fdae61","#fee090")
colors= c("gold","forestgreen","deepskyblue","skyblue")
barplot(t.means4, col=colors,border=NA, beside=F,space=space,las=1)

## Plot k5 ##
par(fig=c(0.075,0.995,0.23,0.40),new=T)
par(mar=c(0,1.5,0.5,1.5))
#colors= c("#a50026","#313695","#4575b4","#fdae61","#fee090")
colors= c("forestgreen","gold","goldenrod","deepskyblue","skyblue")
barplot(t.means5, col=colors,border=NA, beside=F,space=space,las=1)

## Plot k6 ##
par(fig=c(0.075,0.995,0.05,0.22),new=T)
par(mar=c(0,1.5,0.5,1.5))
#colors= c( "#313695","#fdae61","#a50026","#fee090","#4575b4","#ffffbf")
colors= c( "gold","deepskyblue","forestgreen","skyblue","goldenrod","deepskyblue3")
barplot(t.means6, col=colors,border=NA, beside=F,space=space,las=1)

# Labels
par(fig=c(0,1,0,1),new=TRUE)
par(mar=c(0,0,0,0));par(oma=c(0,0,0,0));par(xpd=TRUE)
plot(0:1,0:1,bty="n",xlab="",type="n",xaxt="n",yaxt="n")
#legend("topright",legend=legend_names,cex=0.5,inset=c(0,0),xpd=T)
text(0,0.5,"Admixture Proportions",cex=1.5,srt=90,font=3)
#text(0.5,0.99,"",cex=1.5,font=4)
text(0.02,0.89,"K=2",cex=1.2,srt=90,font=1)
text(0.02,0.69,"K=3",cex=1.2,srt=90,font=1)
text(0.02,0.49,"K=4",cex=1.2,srt=90,font=1)
text(0.02,0.3,"K=5",cex=1.2,srt=90,font=1)
text(0.02,0.1,"K=6",cex=1.2,srt=90,font=1)
text(0.3,0,"BG1",cex=1.5)
text(0.65,0,"AG2",cex=1.5)
text(0.92,0,"AG3",cex=1.5)
dev.off()



## Order Pops
pops<-unique(names[,2])
t.names<-t(names)
qnames<-rbind(t.names,t.means2)
list_pops<-list()
for (j in 1:length(pops)){
  mylist<-list()
  k<-1
  for(i in 1:115){
    
    if (qnames[2,i]==pops[j]){
      mylist[[k]]<-qnames[,i]
      k<-k+1}
  }
  list_pops[[j]]<-do.call(rbind,mylist)
}

# Check Numbers
out<-matrix(nrow=30,ncol=2)
for (i in 1:length(pops)){
  out[i,]<-dim(list_pops[[i]])
}

# Cal1
Cal<-rbind(list_pops[[1]],list_pops[[2]],list_pops[[3]],list_pops[[4]],list_pops[[5]],list_pops[[6]],list_pops[[7]],list_pops[[8]],list_pops[[9]],list_pops[[10]],list_pops[[11]])
unique(Cal[,2])

#Eva
Eva<-rbind(list_pops[[12]],list_pops[[13]],list_pops[[14]],list_pops[[15]],list_pops[[16]],list_pops[[17]],list_pops[[18]],list_pops[[19]],list_pops[[20]],list_pops[[21]],list_pops[[22]])
unique(Eva[,2])            

#anna<-cpe,dop,fcr,lks,ybg
North<-rbind(list_pops[[23]],list_pops[[24]],list_pops[[25]],list_pops[[26]],list_pops[[27]],list_pops[[28]],list_pops[[29]],list_pops[[30]])
unique(North[,2])

ordered_pops<-rbind(Cal,Eva,North)
t.ordered_pops2<-t(ordered_pops[,-c(1:2)])
unique(ordered_pops[,2])