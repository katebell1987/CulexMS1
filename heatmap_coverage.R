#### Set Working Directory ####
setwd("~/Dropbox/UMD_Research/KLB_AAB_GBS/NewAlignment/Kate/")

#### Read in Data ####
dat <- read.csv("all.loci.txtdoublefiltered_variants.mpgl_coverage.csv", header = F) # CSV text file from coverage_calc_bam perl script. Basically a count of the number of reads per individual, per site 
dat[1:5,1:5]
dim(dat)  #115 24,530
# rows = ind.s
# col = counts of reads for each variant and first column is ind IDs
str(dat)
out<-lapply(dat[,-c(1)],as.numeric)
out1<-as.data.frame(out)

host<-read.csv("indiv_ids_hostinfo.csv",header=T)

mean_loci.cover_chick<-as.data.frame(apply(dat[host$Response=="c",-1],2,mean))
mean_loci.cover_human<-as.data.frame(apply(dat[host$Response=="h",-1],2,mean))
comb<-cbind(mean_loci.cover_chick,mean_loci.cover_human)
comb$difference<-comb$`apply(dat[host$Response == "c", -1], 2, mean)`-comb$`apply(dat[host$Response == "h", -1], 2, mean)`

###### Get outliers #####
pval<-read.csv("entropy/Results/pvalue_outliers.csv",header=T)
library(tidyr)
pval_fix<-pval %>%
  separate(entropy_name, 
           into = c("text", "num"), 
           sep = "(?<=[A-Za-z])(?=[0-9])"
  )

just_num<-as.numeric(pval_fix$num)
just_num_fix<-just_num+1 ## add one so locus matches the row number

mean_loci.cover_all<-apply(dat[,-1],2,mean)   #2 = cols - so mean coverage for loci 
all_mean<-mean(mean_loci.cover_all) # distribution is right skewed so let's use quantiles not sd. 
lb_all<-quantile(mean_loci.cover_all, probs = .025)
ub_all<-quantile(mean_loci.cover_all, probs = .975)

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.00    6.00    9.00   10.73   14.00   33.00 
##### Lets make a plot with histogram of differences along y axis and then average coverage along z axis ####
pdf("~/Dropbox/UMD_Research/MS1_prefGBS/Figures/hist_coverageByGroup.pdf",height=8,width=14)
par(mar=c(0,0,0,0))
plot(0:100,0:100,bty="n",ylab="",type="n",yaxt="n",xaxt="n")
par(fig=c(0.05,0.95,0.15,0.95),new=T) #0.1,0.95,0.1,0.95 # x1 x2 y1 y2
#par(mar=c(0,0,0,0)) #bottom, left, top, right #1.5,0,1.5,3
hist(comb$difference,breaks=100,main="")
axis(4, seq(0, 2500, length.out=11), seq(0,50,by=5))
points(x=comb$difference[just_num_fix],y=((comb$`apply(dat[host$Response == "c", -1], 2, mean)`[just_num_fix]*(2500/50))),pch=21,bg="#fb9a9980")
points(x=comb$difference[just_num_fix],y=((comb$`apply(dat[host$Response == "h", -1], 2, mean)`[just_num_fix]*(2500/50))),pch=24,bg="#fb9a9980")
abline(h=(all_mean*(2500/50)),col="#74add1")
abline(h=lb_all*(2500/50),lty='longdash',col="#74add1")
abline(h=ub_all*(2500/50),lty='longdash',col="#74add1")
abline(v=mean(comb$difference),col="#74add1")
abline(v=quantile(comb$difference, probs = .975),col="#74add1",lty="longdash")
abline(v=quantile(comb$difference, probs = .025),col="#74add1",lty="longdash")
mtext("Frequency",side=2,line=2.25)
mtext("Mean Number of Reads Per SNP",side=4,line=2.25)
mtext("Difference in Mean Number of Reads Per SNP\n Between Avian and Human Responders",side=1,line=3.35)
dev.off()


out<-cbind(pval,comb[just_num_fix,])
colnames(out)<-c("chromosome","position","entropy_name","RDA1","RDA2","pvalue","Chick_mean","Human_mean","Difference")

write.table(out,"~/Dropbox/UMD_Research/MS1_prefGBS/RData/locus_id_diffcov.csv",row.names =FALSE,col.names = TRUE,quote=F,sep=",")


###### Get outliers for qvalue  #####
qval<-read.csv("entropy/Results/qvalue_outliers.csv",header=T)
library(tidyr)
qval_fix<-qval %>%
  separate(entropy_name, 
           into = c("text", "num"), 
           sep = "(?<=[A-Za-z])(?=[0-9])"
  )

just_numq<-as.numeric(qval_fix$num)
just_num_fixq<-just_numq+1 ## add one so locus matches the row number

##### Lets make a plot with histogram of differences along y axis and then average coverage along z axis ####
pdf("~/Dropbox/UMD_Research/MS1_prefGBS/Figures/hist_coverageByGroup_qvalues_updated.pdf",height=8,width=14)
par(mar=c(0,0,0,0))
plot(0:100,0:100,bty="n",ylab="",type="n",yaxt="n",xaxt="n")
par(fig=c(0.15,0.85,0.2,0.9),new=T) #0.1,0.95,0.1,0.95 # x1 x2 y1 y2
#par(mar=c(0,0,0,0)) #bottom, left, top, right #1.5,0,1.5,3
hist(comb$difference,breaks=100,main="",cex.axis=2)
axis(4, seq(0, 2500, length.out=11), seq(0,50,by=5),cex.axis=2)
points(x=comb$difference[just_num_fixq],y=((comb$`apply(dat[host$Response == "c", -1], 2, mean)`[just_num_fixq]*(2500/50))),pch=21,bg="#fb9a9980")
points(x=comb$difference[just_num_fixq],y=((comb$`apply(dat[host$Response == "h", -1], 2, mean)`[just_num_fixq]*(2500/50))),pch=24,bg="#fb9a9980")
abline(h=(all_mean*(2500/50)),col="#74add1")
abline(h=lb_all*(2500/50),lty='longdash',col="#74add1")
abline(h=ub_all*(2500/50),lty='longdash',col="#74add1")
abline(v=mean(comb$difference),col="#74add1")
abline(v=quantile(comb$difference, probs = .975),col="#74add1",lty="longdash")
abline(v=quantile(comb$difference, probs = .025),col="#74add1",lty="longdash")
mtext("Frequency",side=2,line=4,cex=2)
mtext("Mean Number of Reads Per SNP",side=4,line=4,cex=2)
mtext("Difference in Mean Number of Reads Per SNP\n Between Avian and Human Responders",side=1,line=5,cex=2)
dev.off()


out<-cbind(qval,comb[just_num_fixq,])
colnames(out)<-c("chromosome","position","entropy_name","RDA1","RDA2","qvalue","Chick_mean","Human_mean","Difference")

write.table(out,"~/Dropbox/UMD_Research/MS1_prefGBS/RData/locus_id_diffcov_qvalue.csv",row.names =FALSE,col.names = TRUE,quote=F,sep=",")


