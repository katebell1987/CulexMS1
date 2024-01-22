################################################################################
#              Code to explore top SNPs from RDA analysis                      #
#                     of Culex preference data                                 #
################################################################################

#### Set working directory, load libraries and functions ####
setwd("~/Dropbox/UMD_Research/KLB_AAB_GBS/NewAlignment/Kate/entropy/Results/")
#install.packages("robust")
library(robust)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("qvalue") 
library(qvalue)
library(vegan)
library(ggplot2)


# Function to conduct a RDA based genome scan from https://github.com/Capblancq/RDA-landscape-genomics/blob/main/src/rdadapt.R #
rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

# function to select SNPs that are more than x SD form the mean, from https://popgen.nescent.org/2018-03-27_RDA_GEA.html #
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

#### Load RDA data ####
load("ConRDA_CulexPref.RData")

#### Read in data for chromosome and positions of SNPs in RDA ####
dat<-read.csv("../../all.loci.txt",header=FALSE) # read in data with chromosome names and positions of each SNP
IDs <- strsplit(as.character(dat[,1]),':') # split chromosome from position
IDs2<-do.call(rbind, IDs) # format

#### Get RDA scores and locus names from entropy ####
locus_scores <- scores(conRda, choices=c(1:2), display="species", scaling="none") # get RDA scores
locus_ent<-as.data.frame(colnames(combg3)) # get SNP names used in entropy

#### Estimate p values and q values (FDR) ####
rda_pq<-rdadapt(conRda,2) # calculate p values and q values (FDR) using function from Capblancq 2021
out<-as.data.frame(cbind(locus_ent,IDs2,rda_pq))
colnames(out)<-c("Entropy_Markername","Chromosome","position","p.values","q.values")
write.table(out,file="full_pqvalues.csv",col.names=TRUE,row.names=FALSE,sep=",",quote=FALSE)

#### Explore SNPs that show the strongest association with preference based on p values less than 0.001 ####
top_p.values<-rda_pq$p.values[which(rda_pq$p.values<0.001)] # what are the SNPs with top p values
RDA_scores_p<-locus_scores[which(rda_pq$p.values<0.001),] # get RDA scores for top SNPs
genome_names_p<-IDs2[which(rda_pq$p.values<0.001),] # get chromosome and positions for top SNPs
ent_names_p<-locus_ent[which(rda_pq$p.values<0.001),] # get names used in entropy for top SNPs
top_pSNPS<-as.data.frame(cbind(genome_names_p,ent_names_p,RDA_scores_p,top_p.values)) # combine into a data frame
colnames(top_pSNPS)<-c("chromosome","position","entropy_name","RDA1","RDA2","pvalue") # fix column names
write.table(top_pSNPS,file="pvalue_outliers.csv",col.names = TRUE,row.names = FALSE,sep=",",quote = FALSE)
geno_p<-combg3[,which(colnames(combg3)%in%ent_names_p)] # get genotype probabilities for top SNPs
geno_p_ind<-data.frame(cbind(ids,geno_p)) # add information about IDs of individuals
write.table(geno_p_ind,file="genotypes_pvalue_outliers.csv",col.names = TRUE,row.names=FALSE,sep=",",quote=FALSE)

#### Explore SNPs that show the strongest association with preference based on q values (FDR) less than 0.05 ####
top_q.values<-rda_pq$q.values[which(rda_pq$q.values<0.05)] # what are the SNPs with top q values
RDA_scores_q<-locus_scores[which(rda_pq$q.values<0.05),] # get RDA scores for top SNPs
genome_names_q<-IDs2[which(rda_pq$q.values<0.05),] # get chromosome and positions for top SNPs
ent_names_q<-locus_ent[which(rda_pq$q.values<0.05),] # get names used in entropy for top SNPs
top_qSNPS<-as.data.frame(cbind(genome_names_q,ent_names_q,RDA_scores_q,top_q.values)) # combine into a data frame
colnames(top_qSNPS)<-c("chromosome","position","entropy_name","RDA1","RDA2","qvalue") # fix column names
write.table(top_qSNPS,file="qvalue_outliers.csv",col.names = TRUE,row.names = FALSE,sep=",",quote = FALSE)
geno_q<-combg3[,which(colnames(combg3)%in%ent_names_q)] # get genotype probabilities for top SNPs
geno_q_ind<-data.frame(cbind(ids,geno_q)) # add information about IDs of individuals
write.table(geno_q_ind,file="genotypes_qvalue_outliers.csv",col.names = TRUE,row.names=FALSE,sep=",",quote=FALSE)

#### Explore SNPs that are more than 3 SD from mean for each RDA axis ####
cand1<- data.frame(outliers(locus_scores[,1],3)) # get scores for SNPs on RDA 1 more than 3SD from mean
ent_names_3sd_1<-as.data.frame(rownames(cand1)) # get entropy names for top SNPs
RDA_scores_3sd_1<-locus_scores[which(locus_ent[,1]%in%ent_names_3sd_1[,1]),2] # get RDA scores for top SNPs for RDA 2
genome_names_3sd_1<-IDs2[which(locus_ent[,1]%in%ent_names_3sd_1[,1]),] # get chromosome and positions for top SNPs
top_3sd_RDA1<-as.data.frame(cbind(genome_names_3sd_1,ent_names_3sd_1,cand1,RDA_scores_3sd_1)) # combine into a data frame
colnames(top_3sd_RDA1)<-c("chromosome","position","entropy_name","RDA1","RDA2") # fix column names
write.table(top_3sd_RDA1,file="3sd_RDA1_outliers.csv",col.names = TRUE,row.names = FALSE,sep=",",quote = FALSE)
geno_3sd_1<-combg3[,which(colnames(combg3)%in%ent_names_3sd_1[,1])] # get genotype probabilities for top SNPs
geno_3sd_1_ind<-data.frame(cbind(ids,geno_3sd_1)) # add information about IDs of individuals
write.table(geno_3sd_1_ind,file="genotypes_3sd_RDA1_outliers.csv",col.names = TRUE,row.names=FALSE,sep=",",quote=FALSE)

cand2<- data.frame(outliers(locus_scores[,2],3)) # get scores for SNPs on RDA 1 more than 3SD from mean
ent_names_3sd_2<-as.data.frame(rownames(cand2)) # get entropy names for top SNPs
RDA_scores_3sd_2<-locus_scores[which(locus_ent[,1]%in%ent_names_3sd_2[,1]),1] # get RDA scores for top SNPs for RDA 2
genome_names_3sd_2<-IDs2[which(locus_ent[,1]%in%ent_names_3sd_2[,1]),] # get chromosome and positions for top SNPs
top_3sd_RDA2<-as.data.frame(cbind(genome_names_3sd_2,ent_names_3sd_2,RDA_scores_3sd_2,cand2)) # combine into a data frame
colnames(top_3sd_RDA2)<-c("chromosome","position","entropy_name","RDA1","RDA2") # fix column names
write.table(top_3sd_RDA2,file="3sd_RDA2_outliers.csv",col.names = TRUE,row.names = FALSE,sep=",",quote = FALSE)
geno_3sd_2<-combg3[,which(colnames(combg3)%in%ent_names_3sd_2[,1])] # get genotype probabilities for top SNPs
geno_3sd_2_ind<-data.frame(cbind(ids,geno_3sd_2)) # add information about IDs of individuals
write.table(geno_3sd_2_ind,file="genotypes_3sd_RDA2_outliers.csv",col.names = TRUE,row.names=FALSE,sep=",",quote=FALSE)

#### Make plots that show p and q values by chromosome position and on RDA ####
outlier_info_plot<- IDs2[,1] # chromosome info for plot
outlier_info_plot[which(rda_pq$q.values<0.05)]<- "FDR outliers" # identify those SNPs with FDR (q) below 0.05
outlier_info_plot<- factor(outlier_info_plot,levels=c("JADDTP010000001.1","JADDTP010000002.1","JADDTP010000003.1","FDR outliers"))
TAB_manhatan<- data.frame(pos = 1:length(IDs2[,1]), pvalues = rda_pq$p.values, Outliers = outlier_info_plot) # make data frame

# Code from https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R #
OBP_OR_info<-read.csv("~/Dropbox/UMD_Research/MS1_prefGBS/RData/TableX_Wstat_chemosens_gene_overlaps-2.csv",header=TRUE)
OBP_OR_info$Chr[OBP_OR_info$Chr=="NC_051861.1"]<-"JADDTP010000001.1" 
OBP_OR_info$Chr[OBP_OR_info$Chr=="NC_051862.1"]<-"JADDTP010000002.1" 
OBP_OR_info$Chr[OBP_OR_info$Chr=="NC_051863.1"]<-"JADDTP010000003.1" 
chr<-c(IDs2[,1],OBP_OR_info$Chr) # chromosome info #24530:24564 = OBP location
pos<-c(IDs2[,2],OBP_OR_info$Start)# positions on chromosome
pos<-as.numeric(pos) 
if(!is.ordered(chr)) { # make sure chr is an ordered factor
  chr <- ordered(chr)
} else {
  chr <- chr[,drop=T]
}
if (any(pos>1e6)) pos<-pos/1e6 # make sure positions are in kbp
posmin <- tapply(pos,chr, min) # calculate absolute genomic position from relative chromosomal positions
posmax <- tapply(pos,chr, max)
posshift <- head(c(0,cumsum(posmax)),-1)
names(posshift) <- levels(chr)
genpos <- pos + posshift[chr]
max_p_outlier<-max(TAB_manhatan$pvalues[TAB_manhatan$Outliers=="FDR outliers"]) # get max p value for q outliers

ypos<-c(7,7,7,7.25,7.5,7.75,8,8.25,7,7.25,7.5,7.75,8,8.25,8.5,8.75,9,7,7,7.25,7,7.25,7,7.25,7.5,7.75,8,8.25,8.5,8.75,7,7.25,7,7.25,7.5,7.75)
ypos_line<-c(6.85,6.85,6.85,6.85,6.85,6.85,6.85,6,6,6,6,6,6,6,6,6,6,6.85,6.5,6.5,6.5,6.5,6,6,6,6,6,6,6,6,6.85,6.85,6.85,6.85,6.85,6.85)
ystart<-rep(0,36)

pdf("~/Dropbox/UMD_Research/MS1_prefGBS/Figures/manhattan_conRDA_pValue_qoutlier_UPDATE.pdf",height=6,width=12)
ggplot(data = TAB_manhatan)+ annotate("segment",x=genpos[24530:24565], y=ystart, xend=genpos[24530:24565], yend=ypos_line,size=0.05,col="grey50") +
  geom_point(aes(x=genpos[1:24529], y=-log10(pvalues),col = Outliers),size=0.85) +
  scale_color_manual(values = c("gray90", "grey84", "grey90","#f8514f")) +
  xlab("Chromosome") + ylab("-log10(p.values)") + scale_x_continuous(breaks=((posmax+posmin)/2+posshift),labels = c("1","2","3"))+
  geom_hline(yintercept=-log10(max_p_outlier), linetype="dashed", color = gray(.80), size=0.6) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 15) + 
  theme(legend.title=element_blank(),legend.position="none", legend.background = element_blank(), panel.grid = element_blank(), legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), legend.text=element_text(size=rel(1.2)), strip.text = element_text(size=11))+ylim(0,10) + annotate("text",x=c(genpos[24530:24537],genpos[24538:24546]-25,genpos[24547],genpos[24548:24549]+15,genpos[24550:24551]+20,genpos[24552:24559]+30,genpos[24560:24565]),y=ypos,label=OBP_OR_info$Gene,size=3) + annotate("segment",x=genpos[24538:24546],xend=genpos[24538:24546]-25,y=6,yend=6,size=0.05)+ annotate("segment",x=genpos[24538:24546]-25,xend=genpos[24538:24546]-25,y=6,yend=6.85,size=0.05)+ annotate("segment",x=genpos[24552:24559],xend=genpos[24552:24559]+30,y=6,yend=6,size=0.05)+annotate("segment",x=genpos[24552:24559]+30,xend=genpos[24552:24559]+30,y=6,yend=6.85,size=0.05)+annotate("segment",x=genpos[24550:24551],xend=genpos[24550:24551]+20,y=6.5,yend=6.5,size=0.05)+annotate("segment",x=genpos[24550:24551]+20,xend=genpos[24550:24551]+20,y=6.5,yend=6.85,size=0.05)+annotate("segment",x=genpos[24548:24549],xend=genpos[24548:24549]+15,y=6.5,yend=6.5,size=0.05)+annotate("segment",x=genpos[24548:24549]+15,xend=genpos[24548:24549]+15,y=6.5,yend=6.85,size=0.05)
dev.off()

## Formatting table for ggplot
locus_scores <- scores(conRda, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%top_pSNPS$entropy_name] <- "P value outliers"
TAB_loci$type[TAB_loci$names%in%top_qSNPS$entropy_name] <- "Q value outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "P value outliers", "Q value outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(conRda, choices=c(1,2), display="bp")) # pull the biplot scores
##TAB_var <- as.data.frame(scores(conRda, choices=c(1,2), display="cn",scaling="none")) # pull the biplot scores
pdf("~/Dropbox/UMD_Research/MS1_prefGBS/Figures/pvalues_conRDA.pdf",width=10,height=8)
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#fb9a99", "#e31a1c")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = c("human","switchers")), size = 2.5) +
  xlab("RDA 1 (1.16%)") + ylab("RDA 2 (0.93%)") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11) +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

dev.off()

pdf("~/Dropbox/UMD_Research/MS1_prefGBS/Figures/pvalues_conRDA_centroids.pdf",width=10,height=9)
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#fb9a99", "#e31a1c")) + geom_text(data = TAB_var, aes(x=RDA1, y=RDA2, label = c("Avian","Human","Switcher")), size = 2.5) +
  xlab("RDA 1 (1.16%)") + ylab("RDA 2 (0.93%)") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11) +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
dev.off()

pdf("~/Dropbox/UMD_Research/MS1_prefGBS/Figures/pvalues_conRDA_nolabel.pdf",width=10,height=9)
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#fb9a99", "#e31a1c"))+
  xlab("RDA 1 (1.16%)") + ylab("RDA 2 (0.93%)") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11) +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
dev.off()

