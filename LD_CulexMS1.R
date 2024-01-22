######## R code to make summaries of LD calculations ##########
# Bootstrapping function --------------------------------------------------
boot.fn <- function(x, N=5000) {
  Int.1 <- replicate(N, mean(sample(x, size= length(x), replace=T)))
  Int.CI <- quantile(Int.1, probs=c(0.025,0.975))
  Int.CI
}

setwd("/Volumes/G-DRIVE mobile USB/LD/")

all_possible_snps<-as.data.frame(matrix(nrow=15,ncol=7))

# First make summaries for all SNPs possible per population, and per behavioral response
setwd("~/Dropbox/UMD_Research/MS1_prefGBS/RData/LD/")
cal_chr1<-read.csv("cal1_chr1.geno.ld",header=T,sep="\t")
cal_chr1_fix<-cal_chr1[cal_chr1[,5]!="NaN",]
all_possible_snps[1,1]<-"BG1"
all_possible_snps[1,2]<-"Chr1"
all_possible_snps[1,3]<-4634
all_possible_snps[1,4]<-length(unique(c(cal_chr1_fix$POS1,cal_chr1_fix$POS2)))
all_possible_snps[1,5]<-mean(cal_chr1_fix$R.2)
out<-boot.fn(cal_chr1_fix$R.2)
all_possible_snps[1,6]<-out[1]
all_possible_snps[1,7]<-out[2]

eva_chr1<-read.csv("eva_chr1.geno.ld",header=T,sep="\t")
eva_chr1_fix<-eva_chr1[eva_chr1$R.2!="NaN",]
all_possible_snps[2,1]<-"AG2"
all_possible_snps[2,2]<-"Chr1"
all_possible_snps[2,3]<-4634
all_possible_snps[2,4]<-length(unique(c(eva_chr1_fix$POS1,eva_chr1_fix$POS2)))
all_possible_snps[2,5]<-mean(eva_chr1_fix$R.2)
out<-boot.fn(eva_chr1_fix$R.2)
all_possible_snps[2,6]<-out[1]
all_possible_snps[2,7]<-out[2]

nor_chr1<-read.csv("nor_chr1.geno.ld",header=T,sep="\t")
nor_chr1_fix<-nor_chr1[nor_chr1$R.2!="NaN",]
all_possible_snps[3,1]<-"AG3"
all_possible_snps[3,2]<-"Chr1"
all_possible_snps[3,3]<-4634
all_possible_snps[3,4]<-length(unique(c(nor_chr1_fix$POS1,nor_chr1_fix$POS2)))
all_possible_snps[3,5]<-mean(nor_chr1_fix$R.2)
out<-boot.fn(nor_chr1_fix$R.2)
all_possible_snps[3,6]<-out[1]
all_possible_snps[3,7]<-out[2]

chick_chr1<-read.csv("chick_chr1.geno.ld",header=T,sep="\t")
chick_chr1_fix<-chick_chr1[chick_chr1$R.2!="NaN",]
all_possible_snps[4,1]<-"Avian Responders"
all_possible_snps[4,2]<-"Chr1"
all_possible_snps[4,3]<-4634
all_possible_snps[4,4]<-length(unique(c(chick_chr1_fix$POS1,chick_chr1_fix$POS2)))
all_possible_snps[4,5]<-mean(chick_chr1_fix$R.2)
out<-boot.fn(chick_chr1_fix$R.2)
all_possible_snps[4,6]<-out[1]
all_possible_snps[4,7]<-out[2]

human_chr1<-read.csv("human_chr1.geno.ld",header=T,sep="\t")
human_chr1_fix<-human_chr1[human_chr1$R.2!="NaN",]
all_possible_snps[5,1]<-"Human Responders"
all_possible_snps[5,2]<-"Chr1"
all_possible_snps[5,3]<-4634
all_possible_snps[5,4]<-length(unique(c(human_chr1_fix$POS1,human_chr1_fix$POS2)))
all_possible_snps[5,5]<-mean(human_chr1_fix$R.2)
out<-boot.fn(human_chr1_fix$R.2)
all_possible_snps[5,6]<-out[1]
all_possible_snps[5,7]<-out[2]

cal_chr2<-read.csv("cal1_chr2.geno.ld",header=T,sep="\t")
cal_chr2_fix<-cal_chr2[cal_chr2$R.2!="NaN",]
all_possible_snps[6,1]<-"BG1"
all_possible_snps[6,2]<-"chr2"
all_possible_snps[6,3]<-11301
all_possible_snps[6,4]<-length(unique(c(cal_chr2_fix$POS1,cal_chr2_fix$POS2)))
all_possible_snps[6,5]<-mean(cal_chr2_fix$R.2)
out<-boot.fn(cal_chr2_fix$R.2)
all_possible_snps[6,6]<-out[1]
all_possible_snps[6,7]<-out[2]

eva_chr2<-read.csv("eva_chr2.geno.ld",header=T,sep="\t")
eva_chr2_fix<-eva_chr2[eva_chr2$R.2!="NaN",]
all_possible_snps[7,1]<-"AG2"
all_possible_snps[7,2]<-"chr2"
all_possible_snps[7,3]<-11301
all_possible_snps[7,4]<-length(unique(c(eva_chr2_fix$POS1,eva_chr2_fix$POS2)))
all_possible_snps[7,5]<-mean(eva_chr2_fix$R.2)
out<-boot.fn(eva_chr2_fix$R.2)
all_possible_snps[7,6]<-out[1]
all_possible_snps[7,7]<-out[2]

nor_chr2<-read.csv("nor_chr2.geno.ld",header=T,sep="\t")
nor_chr2_fix<-nor_chr2[nor_chr2$R.2!="NaN",]
all_possible_snps[8,1]<-"AG3"
all_possible_snps[8,2]<-"chr2"
all_possible_snps[8,3]<-11301
all_possible_snps[8,4]<-length(unique(c(nor_chr2_fix$POS1,nor_chr2_fix$POS2)))
all_possible_snps[8,5]<-mean(nor_chr2_fix$R.2)
out<-boot.fn(nor_chr2_fix$R.2)
all_possible_snps[8,6]<-out[1]
all_possible_snps[8,7]<-out[2]

chick_chr2<-read.csv("chick_chr2.geno.ld",header=T,sep="\t")
chick_chr2_fix<-chick_chr2[chick_chr2$R.2!="NaN",]
all_possible_snps[9,1]<-"Avian Responders"
all_possible_snps[9,2]<-"chr2"
all_possible_snps[9,3]<-11301
all_possible_snps[9,4]<-length(unique(c(chick_chr2_fix$POS1,chick_chr2_fix$POS2)))
all_possible_snps[9,5]<-mean(chick_chr2_fix$R.2)
out<-boot.fn(chick_chr2_fix$R.2)
all_possible_snps[9,6]<-out[1]
all_possible_snps[9,7]<-out[2]

human_chr2<-read.csv("human_chr2.geno.ld",header=T,sep="\t")
human_chr2_fix<-human_chr2[human_chr2$R.2!="NaN",]
all_possible_snps[10,1]<-"Human Responders"
all_possible_snps[10,2]<-"chr2"
all_possible_snps[10,3]<-11301
all_possible_snps[10,4]<-length(unique(c(human_chr2_fix$POS1,human_chr2_fix$POS2)))
all_possible_snps[10,5]<-mean(human_chr2_fix$R.2)
out<-boot.fn(human_chr2_fix$R.2)
all_possible_snps[10,6]<-out[1]
all_possible_snps[10,7]<-out[2]

cal_chr3<-read.csv("cal1_chr3.geno.ld",header=T,sep="\t")
cal_chr3_fix<-cal_chr3[cal_chr3$R.2!="NaN",]
all_possible_snps[11,1]<-"BG1"
all_possible_snps[11,2]<-"chr3"
all_possible_snps[11,3]<-8594
all_possible_snps[11,4]<-length(unique(c(cal_chr3_fix$POS1,cal_chr3_fix$POS2)))
all_possible_snps[11,5]<-mean(cal_chr3_fix$R.2)
out<-boot.fn(cal_chr3_fix$R.2)
all_possible_snps[11,6]<-out[1]
all_possible_snps[11,7]<-out[2]

eva_chr3<-read.csv("eva_chr3.geno.ld",header=T,sep="\t")
eva_chr3_fix<-eva_chr3[eva_chr3$R.2!="NaN",]
all_possible_snps[12,1]<-"AG2"
all_possible_snps[12,2]<-"chr3"
all_possible_snps[12,3]<-8594
all_possible_snps[12,4]<-length(unique(c(eva_chr3_fix$POS1,eva_chr3_fix$POS2)))
all_possible_snps[12,5]<-mean(eva_chr3_fix$R.2)
out<-boot.fn(eva_chr3_fix$R.2)
all_possible_snps[12,6]<-out[1]
all_possible_snps[12,7]<-out[2]

nor_chr3<-read.csv("nor_chr3.geno.ld 2",header=T,sep="\t")
nor_chr3_fix<-nor_chr3[nor_chr3$R.2!="NaN",]
all_possible_snps[13,1]<-"AG3"
all_possible_snps[13,2]<-"chr3"
all_possible_snps[13,3]<-8594
all_possible_snps[13,4]<-length(unique(c(nor_chr3_fix$POS1,nor_chr3_fix$POS2)))
all_possible_snps[13,5]<-mean(nor_chr3_fix$R.2)
out<-boot.fn(nor_chr3_fix$R.2)
all_possible_snps[13,6]<-out[1]
all_possible_snps[13,7]<-out[2]

chick_chr3<-read.csv("chick_chr3.geno.ld",header=T,sep="\t")
chick_chr3_fix<-chick_chr3[chick_chr3$R.2!="NaN",]
all_possible_snps[14,1]<-"Avian Responders"
all_possible_snps[14,2]<-"chr3"
all_possible_snps[14,3]<-8594
all_possible_snps[14,4]<-length(unique(c(chick_chr3_fix$POS1,chick_chr3_fix$POS2)))
all_possible_snps[14,5]<-mean(chick_chr3_fix$R.2)
out<-boot.fn(chick_chr3_fix$R.2)
all_possible_snps[14,6]<-out[1]
all_possible_snps[14,7]<-out[2]

human_chr3<-read.csv("human_chr3.geno.ld",header=T,sep="\t")
human_chr3_fix<-human_chr3[human_chr3$R.2!="NaN",]
all_possible_snps[15,1]<-"Human Responders"
all_possible_snps[15,2]<-"chr3"
all_possible_snps[15,3]<-8594
all_possible_snps[15,4]<-length(unique(c(human_chr3_fix$POS1,human_chr3_fix$POS2)))
all_possible_snps[15,5]<-mean(human_chr3_fix$R.2)
out<-boot.fn(human_chr3_fix$R.2)
all_possible_snps[15,6]<-out[1]
all_possible_snps[15,7]<-out[2]
write.table(all_possible_snps,file="summary_LD_allsnps.csv",row.names=F,col.names=T,quote=F,sep=",")

############ Just Overlaps SNPs ###############

cal_chr1<-read.csv("cal1_chr1.geno.ld",header=T,sep="\t")
cal_chr1_fix<-cal_chr1[cal_chr1[,5]!="NaN",]
eva_chr1<-read.csv("eva_chr1.geno.ld",header=T,sep="\t")
eva_chr1_fix<-eva_chr1[eva_chr1$R.2!="NaN",]
nor_chr1<-read.csv("nor_chr1.geno.ld",header=T,sep="\t")
nor_chr1_fix<-nor_chr1[nor_chr1$R.2!="NaN",]
chick_chr1<-read.csv("chick_chr1.geno.ld",header=T,sep="\t")
chick_chr1_fix<-chick_chr1[chick_chr1$R.2!="NaN",]
human_chr1<-read.csv("human_chr1.geno.ld",header=T,sep="\t")
human_chr1_fix<-human_chr1[human_chr1$R.2!="NaN",]
cal_eva_1<-merge(cal_chr1_fix,eva_chr1_fix,by = c("POS1","POS2"))
length(unique(cal_eva_1$POS1))
length(unique(cal_eva_1$POS2))
colnames(cal_eva_1)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva")
cal_eva_nor_1<-merge(cal_eva_1,nor_chr1_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_1$POS1))
length(unique(cal_eva_nor_1$POS2))
colnames(cal_eva_nor_1)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor")
cal_eva_nor_chick_1<-merge(cal_eva_nor_1,chick_chr1_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_chick_1$POS1))
length(unique(cal_eva_nor_chick_1$POS2))
colnames(cal_eva_nor_chick_1)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor","chr_chick","nind_chick","r2_chick")
cal_eva_nor_chick_human_1<-merge(cal_eva_nor_chick_1,human_chr1_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_chick_human_1$POS1))
length(unique(cal_eva_nor_chick_human_1$POS2))
colnames(cal_eva_nor_chick_human_1)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor","chr_chick","nind_chick","r2_chick","chr_human","nind_human","r2_human")

overlap_snps<-as.data.frame(matrix(nrow=15,ncol=6))
overlap_snps[1,1]<-"BG1"
overlap_snps[1,2]<-"chr1"
overlap_snps[1,3]<-length(unique(c(cal_eva_nor_chick_human_1$POS1,cal_eva_nor_chick_human_1$POS2)))
overlap_snps[1,4]<-mean(cal_eva_nor_chick_human_1$r2_cal)
out<-boot.fn(cal_eva_nor_chick_human_1$r2_cal)
overlap_snps[1,5]<-out[1]
overlap_snps[1,6]<-out[2]
overlap_snps[2,1]<-"AG2"
overlap_snps[2,2]<-"chr1"
overlap_snps[2,3]<-length(unique(c(cal_eva_nor_chick_human_1$POS1,cal_eva_nor_chick_human_1$POS2)))
overlap_snps[2,4]<-mean(cal_eva_nor_chick_human_1$r2_eva)
out<-boot.fn(cal_eva_nor_chick_human_1$r2_eva)
overlap_snps[2,5]<-out[1]
overlap_snps[2,6]<-out[2]
overlap_snps[3,1]<-"AG3"
overlap_snps[3,2]<-"chr1"
overlap_snps[3,3]<-length(unique(c(cal_eva_nor_chick_human_1$POS1,cal_eva_nor_chick_human_1$POS2)))
overlap_snps[3,4]<-mean(cal_eva_nor_chick_human_1$r2_nor)
out<-boot.fn(cal_eva_nor_chick_human_1$r2_nor)
overlap_snps[3,5]<-out[1]
overlap_snps[3,6]<-out[2]
overlap_snps[4,1]<-"Avian Responders"
overlap_snps[4,2]<-"chr1"
overlap_snps[4,3]<-length(unique(c(cal_eva_nor_chick_human_1$POS1,cal_eva_nor_chick_human_1$POS2)))
overlap_snps[4,4]<-mean(cal_eva_nor_chick_human_1$r2_chick)
out<-boot.fn(cal_eva_nor_chick_human_1$r2_chick)
overlap_snps[4,5]<-out[1]
overlap_snps[4,6]<-out[2]
overlap_snps[5,1]<-"Human Responders"
overlap_snps[5,2]<-"chr1"
overlap_snps[5,3]<-length(unique(c(cal_eva_nor_chick_human_1$POS1,cal_eva_nor_chick_human_1$POS2)))
overlap_snps[5,4]<-mean(cal_eva_nor_chick_human_1$r2_human)
out<-boot.fn(cal_eva_nor_chick_human_1$r2_human)
overlap_snps[5,5]<-out[1]
overlap_snps[5,6]<-out[2]

rm(cal_chr1,cal_chr1_fix,cal_eva_1,cal_eva_nor_1,cal_eva_nor_chick_human_1,cal_eva_nor_chick_1,nor_chr1,nor_chr1_fix,eva_chr1,eva_chr1_fix,human_chr1,human_chr1_fix,chick_chr1,chick_chr1_fix)
cal_chr2<-read.csv("cal1_chr2.geno.ld",header=T,sep="\t")
cal_chr2_fix<-cal_chr2[cal_chr2[,5]!="NaN",]
eva_chr2<-read.csv("eva_chr2.geno.ld",header=T,sep="\t")
eva_chr2_fix<-eva_chr2[eva_chr2$R.2!="NaN",]
nor_chr2<-read.csv("nor_chr2.geno.ld",header=T,sep="\t")
nor_chr2_fix<-nor_chr2[nor_chr2$R.2!="NaN",]
chick_chr2<-read.csv("chick_chr2.geno.ld",header=T,sep="\t")
chick_chr2_fix<-chick_chr2[chick_chr2$R.2!="NaN",]
human_chr2<-read.csv("human_chr2.geno.ld",header=T,sep="\t")
human_chr2_fix<-human_chr2[human_chr2$R.2!="NaN",]
rm(cal_chr2,eva_chr2,nor_chr2,chick_chr2,human_chr2)
cal_eva_2<-merge(cal_chr2_fix,eva_chr2_fix,by = c("POS1","POS2"))
length(unique(cal_eva_2$POS1))
length(unique(cal_eva_2$POS2))
rm(cal_chr2_fix,eva_chr2_fix)
colnames(cal_eva_2)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva")
cal_eva_nor_2<-merge(cal_eva_2,nor_chr2_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_2$POS1))
length(unique(cal_eva_nor_2$POS2))
rm(nor_chr2_fix)
colnames(cal_eva_nor_2)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor")
cal_eva_nor_chick_2<-merge(cal_eva_nor_2,chick_chr2_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_chick_2$POS1))
length(unique(cal_eva_nor_chick_2$POS2))
colnames(cal_eva_nor_chick_2)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor","chr_chick","nind_chick","r2_chick")
rm(chick_chr2_fix)
cal_eva_nor_chick_human_2<-merge(cal_eva_nor_chick_2,human_chr2_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_chick_human_2$POS1))
length(unique(cal_eva_nor_chick_human_2$POS2))
colnames(cal_eva_nor_chick_human_2)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor","chr_chick","nind_chick","r2_chick","chr_human","nind_human","r2_human")

overlap_snps[6,1]<-"BG1"
overlap_snps[6,2]<-"chr2"
overlap_snps[6,3]<-length(unique(c(cal_eva_nor_chick_human_2$POS1,cal_eva_nor_chick_human_2$POS2)))
overlap_snps[6,4]<-mean(cal_eva_nor_chick_human_2$r2_cal)
out<-boot.fn(cal_eva_nor_chick_human_2$r2_cal)
overlap_snps[6,5]<-out[1]
overlap_snps[6,6]<-out[2]
overlap_snps[7,1]<-"AG2"
overlap_snps[7,2]<-"chr2"
overlap_snps[7,3]<-length(unique(c(cal_eva_nor_chick_human_2$POS1,cal_eva_nor_chick_human_2$POS2)))
overlap_snps[7,4]<-mean(cal_eva_nor_chick_human_2$r2_eva)
out<-boot.fn(cal_eva_nor_chick_human_2$r2_eva)
overlap_snps[7,5]<-out[1]
overlap_snps[7,6]<-out[2]
overlap_snps[8,1]<-"AG3"
overlap_snps[8,2]<-"chr2"
overlap_snps[8,3]<-length(unique(c(cal_eva_nor_chick_human_2$POS1,cal_eva_nor_chick_human_2$POS2)))
overlap_snps[8,4]<-mean(cal_eva_nor_chick_human_2$r2_nor)
out<-boot.fn(cal_eva_nor_chick_human_2$r2_nor)
overlap_snps[8,5]<-out[1]
overlap_snps[8,6]<-out[2]
overlap_snps[9,1]<-"Avian Responders"
overlap_snps[9,2]<-"chr2"
overlap_snps[9,3]<-length(unique(c(cal_eva_nor_chick_human_2$POS1,cal_eva_nor_chick_human_2$POS2)))
overlap_snps[9,4]<-mean(cal_eva_nor_chick_human_2$r2_chick)
out<-boot.fn(cal_eva_nor_chick_human_2$r2_chick)
overlap_snps[9,5]<-out[1]
overlap_snps[9,6]<-out[2]
overlap_snps[10,1]<-"Human Responders"
overlap_snps[10,2]<-"chr2"
overlap_snps[10,3]<-length(unique(c(cal_eva_nor_chick_human_2$POS1,cal_eva_nor_chick_human_2$POS2)))
overlap_snps[10,4]<-mean(cal_eva_nor_chick_human_2$r2_human)
out<-boot.fn(cal_eva_nor_chick_human_2$r2_human)
overlap_snps[10,5]<-out[1]
overlap_snps[10,6]<-out[2]

rm(cal_eva_2,cal_eva_nor_2,cal_eva_nor_chick_2,cal_eva_nor_chick_human_2)
cal_chr3<-read.csv("cal1_chr3.geno.ld",header=T,sep="\t")
cal_chr3_fix<-cal_chr3[cal_chr3[,5]!="NaN",]
eva_chr3<-read.csv("eva_chr3.geno.ld",header=T,sep="\t")
eva_chr3_fix<-eva_chr3[eva_chr3$R.2!="NaN",]
nor_chr3<-read.csv("nor_chr3.geno.ld 2",header=T,sep="\t")
nor_chr3_fix<-nor_chr3[nor_chr3$R.2!="NaN",]
chick_chr3<-read.csv("chick_chr3.geno.ld",header=T,sep="\t")
chick_chr3_fix<-chick_chr3[chick_chr3$R.2!="NaN",]
human_chr3<-read.csv("human_chr3.geno.ld",header=T,sep="\t")
human_chr3_fix<-human_chr3[human_chr3$R.2!="NaN",]
rm(cal_chr3,eva_chr3,nor_chr3,chick_chr3,human_chr3)
cal_eva_3<-merge(cal_chr3_fix,eva_chr3_fix,by = c("POS1","POS2"))
length(unique(cal_eva_3$POS1))
length(unique(cal_eva_3$POS2))
rm(cal_chr3_fix,eva_chr3_fix)
colnames(cal_eva_3)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva")
cal_eva_nor_3<-merge(cal_eva_3,nor_chr3_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_3$POS1))
length(unique(cal_eva_nor_3$POS2))
rm(nor_chr3_fix)
colnames(cal_eva_nor_3)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor")
cal_eva_nor_chick_3<-merge(cal_eva_nor_3,chick_chr3_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_chick_3$POS1))
length(unique(cal_eva_nor_chick_3$POS2))
colnames(cal_eva_nor_chick_3)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor","chr_chick","nind_chick","r2_chick")
rm(chick_chr3_fix)
cal_eva_nor_chick_human_3<-merge(cal_eva_nor_chick_3,human_chr3_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_chick_human_3$POS1))
length(unique(cal_eva_nor_chick_human_3$POS2))
colnames(cal_eva_nor_chick_human_3)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor","chr_chick","nind_chick","r2_chick","chr_human","nind_human","r2_human")

overlap_snps[11,1]<-"BG1"
overlap_snps[11,2]<-"chr3"
overlap_snps[11,3]<-length(unique(c(cal_eva_nor_chick_human_3$POS1,cal_eva_nor_chick_human_3$POS2)))
overlap_snps[11,4]<-mean(cal_eva_nor_chick_human_3$r2_cal)
out<-boot.fn(cal_eva_nor_chick_human_3$r2_cal)
overlap_snps[11,5]<-out[1]
overlap_snps[11,6]<-out[2]
overlap_snps[12,1]<-"AG2"
overlap_snps[12,2]<-"chr3"
overlap_snps[12,3]<-length(unique(c(cal_eva_nor_chick_human_3$POS1,cal_eva_nor_chick_human_3$POS2)))
overlap_snps[12,4]<-mean(cal_eva_nor_chick_human_3$r2_eva)
out<-boot.fn(cal_eva_nor_chick_human_3$r2_eva)
overlap_snps[12,5]<-out[1]
overlap_snps[12,6]<-out[2]
overlap_snps[13,1]<-"AG3"
overlap_snps[13,2]<-"chr3"
overlap_snps[13,3]<-length(unique(c(cal_eva_nor_chick_human_3$POS1,cal_eva_nor_chick_human_3$POS2)))
overlap_snps[13,4]<-mean(cal_eva_nor_chick_human_3$r2_nor)
out<-boot.fn(cal_eva_nor_chick_human_3$r2_nor)
overlap_snps[13,5]<-out[1]
overlap_snps[13,6]<-out[2]
overlap_snps[14,1]<-"Avian Responders"
overlap_snps[14,2]<-"chr3"
overlap_snps[14,3]<-length(unique(c(cal_eva_nor_chick_human_3$POS1,cal_eva_nor_chick_human_3$POS2)))
overlap_snps[14,4]<-mean(cal_eva_nor_chick_human_3$r2_chick)
out<-boot.fn(cal_eva_nor_chick_human_3$r2_chick)
overlap_snps[14,5]<-out[1]
overlap_snps[14,6]<-out[2]
overlap_snps[15,1]<-"Human Responders"
overlap_snps[15,2]<-"chr3"
overlap_snps[15,3]<-length(unique(c(cal_eva_nor_chick_human_3$POS1,cal_eva_nor_chick_human_3$POS2)))
overlap_snps[15,4]<-mean(cal_eva_nor_chick_human_3$r2_human)
out<-boot.fn(cal_eva_nor_chick_human_3$r2_human)
overlap_snps[15,5]<-out[1]
overlap_snps[15,6]<-out[2]

write.table(overlap_snps,file="summary_LD_overlapsnps.csv",row.names=F,col.names = F,quote=F,sep=" ,")

#### all possible top snps ######
top_snps<-as.data.frame(matrix(nrow=15,ncol=7))
cal_topchr1<-read.csv("cal1_top_chr1.geno.ld",header=T,sep="\t")
cal_topchr1_fix<-cal_topchr1[cal_topchr1[,5]!="NaN",]
top_snps[1,1]<-"BG1"
top_snps[1,2]<-"topchr1"
top_snps[1,3]<-length(unique(c(cal_topchr1$POS1,cal_topchr1$POS2)))
top_snps[1,4]<-length(unique(c(cal_topchr1_fix$POS1,cal_topchr1_fix$POS2)))
top_snps[1,5]<-mean(cal_topchr1_fix$R.2)
out<-boot.fn(cal_topchr1_fix$R.2)
top_snps[1,6]<-out[1]
top_snps[1,7]<-out[2]

eva_topchr1<-read.csv("eva_top_chr1.geno.ld",header=T,sep="\t")
eva_topchr1_fix<-eva_topchr1[eva_topchr1[,5]!="NaN",]
top_snps[2,1]<-"AG2"
top_snps[2,2]<-"topchr1"
top_snps[2,3]<-length(unique(c(eva_topchr1$POS1,eva_topchr1$POS2)))
top_snps[2,4]<-length(unique(c(eva_topchr1_fix$POS1,eva_topchr1_fix$POS2)))
top_snps[2,5]<-mean(eva_topchr1_fix$R.2)
out<-boot.fn(eva_topchr1_fix$R.2)
top_snps[2,6]<-out[1]
top_snps[2,7]<-out[2]

nor_topchr1<-read.csv("nor_top_chr1.geno.ld",header=T,sep="\t")
nor_topchr1_fix<-nor_topchr1[nor_topchr1[,5]!="NaN",]
top_snps[3,1]<-"AG3"
top_snps[3,2]<-"topchr1"
top_snps[3,3]<-length(unique(c(nor_topchr1$POS1,nor_topchr1$POS2)))
top_snps[3,4]<-length(unique(c(nor_topchr1_fix$POS1,nor_topchr1_fix$POS2)))
top_snps[3,5]<-mean(nor_topchr1_fix$R.2)
out<-boot.fn(nor_topchr1_fix$R.2)
top_snps[3,6]<-out[1]
top_snps[3,7]<-out[2]

chick_topchr1<-read.csv("chick_top_chr1.geno.ld",header=T,sep="\t")
chick_topchr1_fix<-chick_topchr1[chick_topchr1[,5]!="NaN",]
top_snps[4,1]<-"Avian Responders"
top_snps[4,2]<-"topchr1"
top_snps[4,3]<-length(unique(c(chick_topchr1$POS1,chick_topchr1$POS2)))
top_snps[4,4]<-length(unique(c(chick_topchr1_fix$POS1,chick_topchr1_fix$POS2)))
top_snps[4,5]<-mean(chick_topchr1_fix$R.2)
out<-boot.fn(chick_topchr1_fix$R.2)
top_snps[4,6]<-out[1]
top_snps[4,7]<-out[2]

human_topchr1<-read.csv("human_top_chr1.geno.ld",header=T,sep="\t")
human_topchr1_fix<-human_topchr1[human_topchr1[,5]!="NaN",]
top_snps[5,1]<-"Human Responders"
top_snps[5,2]<-"topchr1"
top_snps[5,3]<-length(unique(c(human_topchr1$POS1,human_topchr1$POS2)))
top_snps[5,4]<-length(unique(c(human_topchr1_fix$POS1,human_topchr1_fix$POS2)))
top_snps[5,5]<-mean(human_topchr1_fix$R.2)
out<-boot.fn(human_topchr1_fix$R.2)
top_snps[5,6]<-out[1]
top_snps[5,7]<-out[2]

cal_topchr2<-read.csv("cal1_top_chr2.geno.ld",header=T,sep="\t")
cal_topchr2_fix<-cal_topchr2[cal_topchr2[,5]!="NaN",]
top_snps[6,1]<-"BG1"
top_snps[6,2]<-"topchr2"
top_snps[6,3]<-length(unique(c(cal_topchr2$POS1,cal_topchr2$POS2)))
top_snps[6,4]<-length(unique(c(cal_topchr2_fix$POS1,cal_topchr2_fix$POS2)))
top_snps[6,5]<-mean(cal_topchr2_fix$R.2)
out<-boot.fn(cal_topchr2_fix$R.2)
top_snps[6,6]<-out[1]
top_snps[6,7]<-out[2]

eva_topchr2<-read.csv("eva_top_chr2.geno.ld",header=T,sep="\t")
eva_topchr2_fix<-eva_topchr2[eva_topchr2[,5]!="NaN",]
top_snps[7,1]<-"AG2"
top_snps[7,2]<-"topchr2"
top_snps[7,3]<-length(unique(c(eva_topchr2$POS1,eva_topchr2$POS2)))
top_snps[7,4]<-length(unique(c(eva_topchr2_fix$POS1,eva_topchr2_fix$POS2)))
top_snps[7,5]<-mean(eva_topchr2_fix$R.2)
out<-boot.fn(eva_topchr2_fix$R.2)
top_snps[7,6]<-out[1]
top_snps[7,7]<-out[2]

nor_topchr2<-read.csv("nor_top_chr2.geno.ld",header=T,sep="\t")
nor_topchr2_fix<-nor_topchr2[nor_topchr2[,5]!="NaN",]
top_snps[8,1]<-"AG3"
top_snps[8,2]<-"topchr2"
top_snps[8,3]<-length(unique(c(nor_topchr2$POS1,nor_topchr2$POS2)))
top_snps[8,4]<-length(unique(c(nor_topchr2_fix$POS1,nor_topchr2_fix$POS2)))
top_snps[8,5]<-mean(nor_topchr2_fix$R.2)
out<-boot.fn(nor_topchr2_fix$R.2)
top_snps[8,6]<-out[1]
top_snps[8,7]<-out[2]

chick_topchr2<-read.csv("chick_top_chr2.geno.ld",header=T,sep="\t")
chick_topchr2_fix<-chick_topchr2[chick_topchr2[,5]!="NaN",]
top_snps[9,1]<-"Avian Responders"
top_snps[9,2]<-"topchr2"
top_snps[9,3]<-length(unique(c(chick_topchr2$POS1,chick_topchr2$POS2)))
top_snps[9,4]<-length(unique(c(chick_topchr2_fix$POS1,chick_topchr2_fix$POS2)))
top_snps[9,5]<-mean(chick_topchr2_fix$R.2)
out<-boot.fn(chick_topchr2_fix$R.2)
top_snps[9,6]<-out[1]
top_snps[9,7]<-out[2]

human_topchr2<-read.csv("human_top_chr2.geno.ld",header=T,sep="\t")
human_topchr2_fix<-human_topchr2[human_topchr2[,5]!="NaN",]
top_snps[10,1]<-"Human Responders"
top_snps[10,2]<-"topchr2"
top_snps[10,3]<-length(unique(c(human_topchr2$POS1,human_topchr2$POS2)))
top_snps[10,4]<-length(unique(c(human_topchr2_fix$POS1,human_topchr2_fix$POS2)))
top_snps[10,5]<-mean(human_topchr2_fix$R.2)
out<-boot.fn(human_topchr2_fix$R.2)
top_snps[10,6]<-out[1]
top_snps[10,7]<-out[2]

cal_topchr3<-read.csv("cal1_top_chr3.geno.ld",header=T,sep="\t")
cal_topchr3_fix<-cal_topchr3[cal_topchr3[,5]!="NaN",]
top_snps[11,1]<-"BG1"
top_snps[11,2]<-"topchr3"
top_snps[11,3]<-length(unique(c(cal_topchr3$POS1,cal_topchr3$POS2)))
top_snps[11,4]<-length(unique(c(cal_topchr3_fix$POS1,cal_topchr3_fix$POS2)))
top_snps[11,5]<-mean(cal_topchr3_fix$R.2)
out<-boot.fn(cal_topchr3_fix$R.2)
top_snps[11,6]<-out[1]
top_snps[11,7]<-out[2]

eva_topchr3<-read.csv("eva_top_chr3.geno.ld",header=T,sep="\t")
eva_topchr3_fix<-eva_topchr3[eva_topchr3[,5]!="NaN",]
top_snps[12,1]<-"AG2"
top_snps[12,2]<-"topchr3"
top_snps[12,3]<-length(unique(c(eva_topchr3$POS1,eva_topchr3$POS2)))
top_snps[12,4]<-length(unique(c(eva_topchr3_fix$POS1,eva_topchr3_fix$POS2)))
top_snps[12,5]<-mean(eva_topchr3_fix$R.2)
out<-boot.fn(eva_topchr3_fix$R.2)
top_snps[12,6]<-out[1]
top_snps[12,7]<-out[2]

nor_topchr3<-read.csv("nor_top_chr3.geno.ld",header=T,sep="\t")
nor_topchr3_fix<-nor_topchr3[nor_topchr3[,5]!="NaN",]
top_snps[13,1]<-"AG3"
top_snps[13,2]<-"topchr3"
top_snps[13,3]<-length(unique(c(nor_topchr3$POS1,nor_topchr3$POS2)))
top_snps[13,4]<-length(unique(c(nor_topchr3_fix$POS1,nor_topchr3_fix$POS2)))
top_snps[13,5]<-mean(nor_topchr3_fix$R.2)
out<-boot.fn(nor_topchr3_fix$R.2)
top_snps[13,6]<-out[1]
top_snps[13,7]<-out[2]

chick_topchr3<-read.csv("chick_top_chr3.geno.ld",header=T,sep="\t")
chick_topchr3_fix<-chick_topchr3[chick_topchr3[,5]!="NaN",]
top_snps[14,1]<-"Avian Responders"
top_snps[14,2]<-"topchr3"
top_snps[14,3]<-length(unique(c(chick_topchr3$POS1,chick_topchr3$POS2)))
top_snps[14,4]<-length(unique(c(chick_topchr3_fix$POS1,chick_topchr3_fix$POS2)))
top_snps[14,5]<-mean(chick_topchr3_fix$R.2)
out<-boot.fn(chick_topchr3_fix$R.2)
top_snps[14,6]<-out[1]
top_snps[14,7]<-out[2]

human_topchr3<-read.csv("human_top_chr3.geno.ld",header=T,sep="\t")
human_topchr3_fix<-human_topchr3[human_topchr3[,5]!="NaN",]
top_snps[15,1]<-"Human Responders"
top_snps[15,2]<-"topchr3"
top_snps[15,3]<-length(unique(c(human_topchr3$POS1,human_topchr3$POS2)))
top_snps[15,4]<-length(unique(c(human_topchr3_fix$POS1,human_topchr3_fix$POS2)))
top_snps[15,5]<-mean(human_topchr3_fix$R.2)
out<-boot.fn(human_topchr3_fix$R.2)
top_snps[15,6]<-out[1]
top_snps[15,7]<-out[2]

write.table(top_snps,file="summary_LD_topsnps.csv",row.names=F,col.names = F,quote=F,sep=" ,")

############ Overlap top SNPs ###############
overlap_top_snps<-as.data.frame(matrix(nrow=15,ncol=6))

cal_top_chr1<-read.csv("cal1_top_chr1.geno.ld",header=T,sep="\t")
cal_top_chr1_fix<-cal_top_chr1[cal_top_chr1[,5]!="NaN",]
eva_top_chr1<-read.csv("eva_top_chr1.geno.ld",header=T,sep="\t")
eva_top_chr1_fix<-eva_top_chr1[eva_top_chr1$R.2!="NaN",]
nor_top_chr1<-read.csv("nor_top_chr1.geno.ld",header=T,sep="\t")
nor_top_chr1_fix<-nor_top_chr1[nor_top_chr1$R.2!="NaN",]
chick_top_chr1<-read.csv("chick_top_chr1.geno.ld",header=T,sep="\t")
chick_top_chr1_fix<-chick_top_chr1[chick_top_chr1$R.2!="NaN",]
human_top_chr1<-read.csv("human_top_chr1.geno.ld",header=T,sep="\t")
human_top_chr1_fix<-human_top_chr1[human_top_chr1$R.2!="NaN",]
cal_eva_top_1<-merge(cal_top_chr1_fix,eva_top_chr1_fix,by = c("POS1","POS2"))
length(unique(cal_eva_top_1$POS1))
length(unique(cal_eva_top_1$POS2))
colnames(cal_eva_top_1)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva")
cal_eva_nor_top_1<-merge(cal_eva_top_1,nor_top_chr1_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_top_1$POS1))
length(unique(cal_eva_nor_top_1$POS2))
colnames(cal_eva_nor_top_1)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor")
cal_eva_nor_chick_top_1<-merge(cal_eva_nor_top_1,chick_top_chr1_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_chick_top_1$POS1))
length(unique(cal_eva_nor_chick_top_1$POS2))
colnames(cal_eva_nor_chick_top_1)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor","chr_chick","nind_chick","r2_chick")
cal_eva_nor_chick_human_top_1<-merge(cal_eva_nor_chick_top_1,human_top_chr1_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_chick_human_top_1$POS1))
length(unique(cal_eva_nor_chick_human_top_1$POS2))
colnames(cal_eva_nor_chick_human_top_1)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor","chr_chick","nind_chick","r2_chick","chr_human","nind_human","r2_human")

overlap_top_snps<-as.data.frame(matrix(nrow=15,ncol=6))
overlap_top_snps[1,1]<-"BG1"
overlap_top_snps[1,2]<-"chr1"
overlap_top_snps[1,3]<-length(unique(c(cal_eva_nor_chick_human_top_1$POS1,cal_eva_nor_chick_human_top_1$POS2)))
overlap_top_snps[1,4]<-mean(cal_eva_nor_chick_human_top_1$r2_cal)
out<-boot.fn(cal_eva_nor_chick_human_top_1$r2_cal)
overlap_top_snps[1,5]<-out[1]
overlap_top_snps[1,6]<-out[2]
overlap_top_snps[2,1]<-"AG2"
overlap_top_snps[2,2]<-"chr1"
overlap_top_snps[2,3]<-length(unique(c(cal_eva_nor_chick_human_top_1$POS1,cal_eva_nor_chick_human_top_1$POS2)))
overlap_top_snps[2,4]<-mean(cal_eva_nor_chick_human_top_1$r2_eva)
out<-boot.fn(cal_eva_nor_chick_human_top_1$r2_eva)
overlap_top_snps[2,5]<-out[1]
overlap_top_snps[2,6]<-out[2]
overlap_top_snps[3,1]<-"AG3"
overlap_top_snps[3,2]<-"chr1"
overlap_top_snps[3,3]<-length(unique(c(cal_eva_nor_chick_human_top_1$POS1,cal_eva_nor_chick_human_top_1$POS2)))
overlap_top_snps[3,4]<-mean(cal_eva_nor_chick_human_top_1$r2_nor)
out<-boot.fn(cal_eva_nor_chick_human_top_1$r2_nor)
overlap_top_snps[3,5]<-out[1]
overlap_top_snps[3,6]<-out[2]
overlap_top_snps[4,1]<-"Avian Responders"
overlap_top_snps[4,2]<-"chr1"
overlap_top_snps[4,3]<-length(unique(c(cal_eva_nor_chick_human_top_1$POS1,cal_eva_nor_chick_human_top_1$POS2)))
overlap_top_snps[4,4]<-mean(cal_eva_nor_chick_human_top_1$r2_chick)
out<-boot.fn(cal_eva_nor_chick_human_top_1$r2_chick)
overlap_top_snps[4,5]<-out[1]
overlap_top_snps[4,6]<-out[2]
overlap_top_snps[5,1]<-"Human Responders"
overlap_top_snps[5,2]<-"chr1"
overlap_top_snps[5,3]<-length(unique(c(cal_eva_nor_chick_human_top_1$POS1,cal_eva_nor_chick_human_top_1$POS2)))
overlap_top_snps[5,4]<-mean(cal_eva_nor_chick_human_top_1$r2_human)
out<-boot.fn(cal_eva_nor_chick_human_top_1$r2_human)
overlap_top_snps[5,5]<-out[1]
overlap_top_snps[5,6]<-out[2]

rm(cal_top_chr1,cal_top_chr1_fix,cal_eva_top_1,cal_eva_nor_top_1,cal_eva_nor_chick_human_top_1,cal_eva_nor_chick_top_1,nor_top_chr1,nor_top_chr1_fix,eva_top_chr1,eva_top_chr1_fix,human_top_chr1,human_top_chr1_fix,chick_top_chr1,chick_top_chr1_fix)
cal_top_chr2<-read.csv("cal1_top_chr2.geno.ld",header=T,sep="\t")
cal_top_chr2_fix<-cal_top_chr2[cal_top_chr2[,5]!="NaN",]
eva_top_chr2<-read.csv("eva_top_chr2.geno.ld",header=T,sep="\t")
eva_top_chr2_fix<-eva_top_chr2[eva_top_chr2$R.2!="NaN",]
nor_top_chr2<-read.csv("nor_top_chr2.geno.ld",header=T,sep="\t")
nor_top_chr2_fix<-nor_top_chr2[nor_top_chr2$R.2!="NaN",]
chick_top_chr2<-read.csv("chick_top_chr2.geno.ld",header=T,sep="\t")
chick_top_chr2_fix<-chick_top_chr2[chick_top_chr2$R.2!="NaN",]
human_top_chr2<-read.csv("human_top_chr2.geno.ld",header=T,sep="\t")
human_top_chr2_fix<-human_top_chr2[human_top_chr2$R.2!="NaN",]
rm(cal_top_chr2,eva_top_chr2,nor_top_chr2,chick_top_chr2,human_top_chr2)
cal_eva_top_2<-merge(cal_top_chr2_fix,eva_top_chr2_fix,by = c("POS1","POS2"))
length(unique(cal_eva_top_2$POS1))
length(unique(cal_eva_top_2$POS2))
rm(cal_top_chr2_fix,eva_top_chr2_fix)
colnames(cal_eva_top_2)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva")
cal_eva_nor_top_2<-merge(cal_eva_top_2,nor_top_chr2_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_top_2$POS1))
length(unique(cal_eva_nor_top_2$POS2))
rm(nor_top_chr2_fix)
colnames(cal_eva_nor_top_2)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor")
cal_eva_nor_chick_top_2<-merge(cal_eva_nor_top_2,chick_top_chr2_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_chick_top_2$POS1))
length(unique(cal_eva_nor_chick_top_2$POS2))
colnames(cal_eva_nor_chick_top_2)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor","chr_chick","nind_chick","r2_chick")
rm(chick_top_chr2_fix)
cal_eva_nor_chick_human_top_2<-merge(cal_eva_nor_chick_top_2,human_top_chr2_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_chick_human_top_2$POS1))
length(unique(cal_eva_nor_chick_human_top_2$POS2))
colnames(cal_eva_nor_chick_human_top_2)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor","chr_chick","nind_chick","r2_chick","chr_human","nind_human","r2_human")

overlap_top_snps[6,1]<-"BG1"
overlap_top_snps[6,2]<-"chr2"
overlap_top_snps[6,3]<-length(unique(c(cal_eva_nor_chick_human_top_2$POS1,cal_eva_nor_chick_human_top_2$POS2)))
overlap_top_snps[6,4]<-mean(cal_eva_nor_chick_human_top_2$r2_cal)
out<-boot.fn(cal_eva_nor_chick_human_top_2$r2_cal)
overlap_top_snps[6,5]<-out[1]
overlap_top_snps[6,6]<-out[2]
overlap_top_snps[7,1]<-"AG2"
overlap_top_snps[7,2]<-"chr2"
overlap_top_snps[7,3]<-length(unique(c(cal_eva_nor_chick_human_top_2$POS1,cal_eva_nor_chick_human_top_2$POS2)))
overlap_top_snps[7,4]<-mean(cal_eva_nor_chick_human_top_2$r2_eva)
out<-boot.fn(cal_eva_nor_chick_human_top_2$r2_eva)
overlap_top_snps[7,5]<-out[1]
overlap_top_snps[7,6]<-out[2]
overlap_top_snps[8,1]<-"AG3"
overlap_top_snps[8,2]<-"chr2"
overlap_top_snps[8,3]<-length(unique(c(cal_eva_nor_chick_human_top_2$POS1,cal_eva_nor_chick_human_top_2$POS2)))
overlap_top_snps[8,4]<-mean(cal_eva_nor_chick_human_top_2$r2_nor)
out<-boot.fn(cal_eva_nor_chick_human_top_2$r2_nor)
overlap_top_snps[8,5]<-out[1]
overlap_top_snps[8,6]<-out[2]
overlap_top_snps[9,1]<-"Avian Responders"
overlap_top_snps[9,2]<-"chr2"
overlap_top_snps[9,3]<-length(unique(c(cal_eva_nor_chick_human_top_2$POS1,cal_eva_nor_chick_human_top_2$POS2)))
overlap_top_snps[9,4]<-mean(cal_eva_nor_chick_human_top_2$r2_chick)
out<-boot.fn(cal_eva_nor_chick_human_top_2$r2_chick)
overlap_top_snps[9,5]<-out[1]
overlap_top_snps[9,6]<-out[2]
overlap_top_snps[10,1]<-"Human Responders"
overlap_top_snps[10,2]<-"chr2"
overlap_top_snps[10,3]<-length(unique(c(cal_eva_nor_chick_human_top_2$POS1,cal_eva_nor_chick_human_top_2$POS2)))
overlap_top_snps[10,4]<-mean(cal_eva_nor_chick_human_top_2$r2_human)
out<-boot.fn(cal_eva_nor_chick_human_top_2$r2_human)
overlap_top_snps[10,5]<-out[1]
overlap_top_snps[10,6]<-out[2]

rm(cal_eva_top_2,cal_eva_nor_top_2,cal_eva_nor_chick_top_2,cal_eva_nor_chick_human_top_2)
cal_top_chr3<-read.csv("cal1_top_chr3.geno.ld",header=T,sep="\t")
cal_top_chr3_fix<-cal_top_chr3[cal_top_chr3[,5]!="NaN",]
eva_top_chr3<-read.csv("eva_top_chr3.geno.ld",header=T,sep="\t")
eva_top_chr3_fix<-eva_top_chr3[eva_top_chr3$R.2!="NaN",]
nor_top_chr3<-read.csv("nor_top_chr3.geno.ld",header=T,sep="\t")
nor_top_chr3_fix<-nor_top_chr3[nor_top_chr3$R.2!="NaN",]
chick_top_chr3<-read.csv("chick_top_chr3.geno.ld",header=T,sep="\t")
chick_top_chr3_fix<-chick_top_chr3[chick_top_chr3$R.2!="NaN",]
human_top_chr3<-read.csv("human_top_chr3.geno.ld",header=T,sep="\t")
human_top_chr3_fix<-human_top_chr3[human_top_chr3$R.2!="NaN",]
rm(cal_top_chr3,eva_top_chr3,nor_top_chr3,chick_top_chr3,human_top_chr3)
cal_eva_top_3<-merge(cal_top_chr3_fix,eva_top_chr3_fix,by = c("POS1","POS2"))
length(unique(cal_eva_top_3$POS1))
length(unique(cal_eva_top_3$POS2))
rm(cal_top_chr3_fix,eva_top_chr3_fix)
colnames(cal_eva_top_3)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva")
cal_eva_nor_top_3<-merge(cal_eva_top_3,nor_top_chr3_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_top_3$POS1))
length(unique(cal_eva_nor_top_3$POS2))
rm(nor_top_chr3_fix)
colnames(cal_eva_nor_top_3)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor")
cal_eva_nor_chick_top_3<-merge(cal_eva_nor_top_3,chick_top_chr3_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_chick_top_3$POS1))
length(unique(cal_eva_nor_chick_top_3$POS2))
colnames(cal_eva_nor_chick_top_3)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor","chr_chick","nind_chick","r2_chick")
rm(chick_top_chr3_fix)
cal_eva_nor_chick_human_top_3<-merge(cal_eva_nor_chick_top_3,human_top_chr3_fix,by=c("POS1","POS2"))
length(unique(cal_eva_nor_chick_human_top_3$POS1))
length(unique(cal_eva_nor_chick_human_top_3$POS2))
colnames(cal_eva_nor_chick_human_top_3)<-c("POS1","POS2","chr_cal","nind_cal","r2_cal","chr_eva","nind_eva","r2_eva","chr_nor","nind_nor","r2_nor","chr_chick","nind_chick","r2_chick","chr_human","nind_human","r2_human")

overlap_top_snps[11,1]<-"BG1"
overlap_top_snps[11,2]<-"chr3"
overlap_top_snps[11,3]<-length(unique(c(cal_eva_nor_chick_human_top_3$POS1,cal_eva_nor_chick_human_top_3$POS2)))
overlap_top_snps[11,4]<-mean(cal_eva_nor_chick_human_top_3$r2_cal)
out<-boot.fn(cal_eva_nor_chick_human_top_3$r2_cal)
overlap_top_snps[11,5]<-out[1]
overlap_top_snps[11,6]<-out[2]
overlap_top_snps[12,1]<-"AG2"
overlap_top_snps[12,2]<-"chr3"
overlap_top_snps[12,3]<-length(unique(c(cal_eva_nor_chick_human_top_3$POS1,cal_eva_nor_chick_human_top_3$POS2)))
overlap_top_snps[12,4]<-mean(cal_eva_nor_chick_human_top_3$r2_eva)
out<-boot.fn(cal_eva_nor_chick_human_top_3$r2_eva)
overlap_top_snps[12,5]<-out[1]
overlap_top_snps[12,6]<-out[2]
overlap_top_snps[13,1]<-"AG3"
overlap_top_snps[13,2]<-"chr3"
overlap_top_snps[13,3]<-length(unique(c(cal_eva_nor_chick_human_top_3$POS1,cal_eva_nor_chick_human_top_3$POS2)))
overlap_top_snps[13,4]<-mean(cal_eva_nor_chick_human_top_3$r2_nor)
out<-boot.fn(cal_eva_nor_chick_human_top_3$r2_nor)
overlap_top_snps[13,5]<-out[1]
overlap_top_snps[13,6]<-out[2]
overlap_top_snps[14,1]<-"Avian Responders"
overlap_top_snps[14,2]<-"chr3"
overlap_top_snps[14,3]<-length(unique(c(cal_eva_nor_chick_human_top_3$POS1,cal_eva_nor_chick_human_top_3$POS2)))
overlap_top_snps[14,4]<-mean(cal_eva_nor_chick_human_top_3$r2_chick)
out<-boot.fn(cal_eva_nor_chick_human_top_3$r2_chick)
overlap_top_snps[14,5]<-out[1]
overlap_top_snps[14,6]<-out[2]
overlap_top_snps[15,1]<-"Human Responders"
overlap_top_snps[15,2]<-"chr3"
overlap_top_snps[15,3]<-length(unique(c(cal_eva_nor_chick_human_top_3$POS1,cal_eva_nor_chick_human_top_3$POS2)))
overlap_top_snps[15,4]<-mean(cal_eva_nor_chick_human_top_3$r2_human)
out<-boot.fn(cal_eva_nor_chick_human_top_3$r2_human)
overlap_top_snps[15,5]<-out[1]
overlap_top_snps[15,6]<-out[2]

write.table(overlap_top_snps,file="summary_LD_overlaptopsnps.csv",row.names=F,col.names = F,quote=F,sep=" ,")

mean(dat$cal1_r2)
out<-boot.fn(dat$cal1_r2)
out
