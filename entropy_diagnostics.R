# comment out the k's you don't need
# automated diagnostics collection and t.medians examination / fix sum_to_1 begins around line 860
# prints .csv's for: ESSs for combined chains, GR stats, sum_to_1, sum_to_1 fix
# saves fixed t.medians as list in .RData object
# set max k here first:
k=6 # set max k
library(coda)
# k=2
q2c0<-read.csv("q2c0_mcmc.txt",header=F);q2c0mcmc<-as.mcmc(q2c0[,-1])
nind <- dim(q2c0)[1]/2
q2c1<-read.csv("q2c1_mcmc.txt",header=F);q2c1mcmc<-as.mcmc(q2c1[,-1])
gelmanStats_k2<-numeric()
for (i in 1:nind){
  gelmanStats_k2[i]<-matrix(unlist(gelman.diag(mcmc.list(as.mcmc(q2c0mcmc[i,]),as.mcmc(q2c1mcmc[i,])))))[1]
}
tq2c0<-t(q2c0mcmc);ESS_k2c0<-numeric()
for (i in 1:nind){
  ESS_k2c0[i]<-effectiveSize(tq2c0[,i]) 
}
tq2c1<-t(q2c1mcmc);ESS_k2c1<-numeric()
for (i in 1:nind){
  ESS_k2c1[i]<-effectiveSize(tq2c1[,i]) 
}

# k=3
q3c0<-read.csv("q3c0_mcmc.txt",header=F);q3c0mcmc<-as.mcmc(q3c0[,-1])
q3c1<-read.csv("q3c1_mcmc.txt",header=F);q3c1mcmc<-as.mcmc(q3c1[,-1])
gelmanStats_k3<-numeric() 
for (i in 1:nind){
  gelmanStats_k3[i]<-matrix(unlist(gelman.diag(mcmc.list(as.mcmc(q3c0mcmc[i,]),as.mcmc(q3c1mcmc[i,])))))[1]
}
tq3c0<-t(q3c0mcmc);ESS_k3c0<-numeric()
for (i in 1:nind){
  ESS_k3c0[i]<-effectiveSize(tq3c0[,i]) 
}
tq3c1<-t(q3c1mcmc);ESS_k3c1<-numeric()
for (i in 1:nind){
  ESS_k3c1[i]<-effectiveSize(tq3c1[,i]) 
}
# k=4
q4c0<-read.csv("q4c0_mcmc.txt",header=F);q4c0mcmc<-as.mcmc(q4c0[,-1]) 
q4c1<-read.csv("q4c1_mcmc.txt",header=F);q4c1mcmc<-as.mcmc(q4c1[,-1])
gelmanStats_k4<-numeric() 
for (i in 1:nind){
  gelmanStats_k4[i]<-matrix(unlist(gelman.diag(mcmc.list(as.mcmc(q4c0mcmc[i,]),as.mcmc(q4c1mcmc[i,])))))[1]
}
tq4c0<-t(q4c0mcmc);ESS_k4c0<-numeric()
for (i in 1:nind){
  ESS_k4c0[i]<-effectiveSize(tq4c0[,i]) 
}
tq4c1<-t(q4c1mcmc);ESS_k4c1<-numeric()
for (i in 1:nind){
  ESS_k4c1[i]<-effectiveSize(tq4c1[,i]) 
}
# k=5
q5c0<-read.csv("q5c0_mcmc.txt",header=F);q5c0mcmc<-as.mcmc(q5c0[,-1]) 
q5c1<-read.csv("q5c1_mcmc.txt",header=F);q5c1mcmc<-as.mcmc(q5c1[,-1])
gelmanStats_k5<-numeric() 
for (i in 1:nind){
  gelmanStats_k5[i]<-matrix(unlist(gelman.diag(mcmc.list(as.mcmc(q5c0mcmc[i,]),as.mcmc(q5c1mcmc[i,])))))[1]
}
tq5c0<-t(q5c0mcmc);ESS_k5c0<-numeric()
for (i in 1:nind){
  ESS_k5c0[i]<-effectiveSize(tq5c0[,i]) 
}
tq5c1<-t(q5c1mcmc);ESS_k5c1<-numeric()
for (i in 1:nind){
  ESS_k5c1[i]<-effectiveSize(tq5c1[,i]) 
}
# k=6
q6c0<-read.csv("q6c0_mcmc.txt",header=F);q6c0mcmc<-as.mcmc(q6c0[,-1]) 
q6c1<-read.csv("q6c1_mcmc.txt",header=F);q6c1mcmc<-as.mcmc(q6c1[,-1])
gelmanStats_k6<-numeric() 
for (i in 1:nind){
  gelmanStats_k6[i]<-matrix(unlist(gelman.diag(mcmc.list(as.mcmc(q6c0mcmc[i,]),as.mcmc(q6c1mcmc[i,])))))[1]
}
tq6c0<-t(q6c0mcmc);ESS_k6c0<-numeric()
for (i in 1:nind){
  ESS_k6c0[i]<-effectiveSize(tq6c0[,i]) 
}
tq6c1<-t(q6c1mcmc);ESS_k6c1<-numeric()
for (i in 1:nind){
  ESS_k6c1[i]<-effectiveSize(tq6c1[,i]) 
}
######automations
# ESS's - combined chains
comboESSz<-matrix(NA,nrow = 6,ncol = 6)
for(i in 2:k){
  combo<-(eval(parse(text=paste0("ESS_k",i,"c0",sep=""))))+(eval(parse(text=paste0("ESS_k",i,"c1",sep=""))))
  comboESSz[i,]<-round(summary(combo),2)
}
comboESSz<-comboESSz[-1,]
colnames(comboESSz)<-c("min","1st","median","mean","3rd","max")
write.csv(comboESSz,"ESS_combinedChains.csv",row.names = F,quote = F)
# GR
gelmanList<-list()
for(i in 2:k){
  gelmanList[[i]]<-eval(parse(text=paste0("gelmanStats_k",i,sep="")))
}
gelmanSumZ<-matrix(NA,nrow = length(gelmanList),ncol = 6)
for(i in 1:length(gelmanList)){
  gelmanSumZ[i,]<-round(as.numeric(summary(gelmanList[[i]])),6)
}
gelmanSumZ<-gelmanSumZ[-1,]
colnames(gelmanSumZ)<-c("min","1st","median","mean","3rd","max")
write.csv(gelmanSumZ,"gelmanStats.csv",row.names = F,quote = F)
save.image("diagnostics.RData")


