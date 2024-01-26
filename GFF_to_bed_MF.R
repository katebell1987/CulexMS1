###This is my script to prep the gff files for bedtools
###07142023 M. Fritz

###### PREPPING WORKSPACE ######
setwd("~/Desktop/Bell_window_sims/")

gff <- read.table("2023_03_12_sortMergedManualRefSeqAnnotation.gtf", sep="\t", header=FALSE, comment.char="#",
                  na.strings=".", stringsAsFactors=FALSE,
                  quote="", fill=FALSE)


###### GFF TO BEDFILE ######
gff_genes_withMito <- subset(gff, V3 == "transcript") #getting mRNA coordinates only
#gff_genes_withMito$V1 <- as.factor(gff_genes_withMito$V1)

#getting Ors and OBPs
gff_genes <- subset(gff_genes_withMito, V1 == "NC_051861.1" | V1 == "NC_051862.1" | V1 == "NC_051863.1")
summary(gff_genes$V1)#gives tot genes per chromosome

ids <- gsub("transcript_id", "", sapply(strsplit(gff_genes$V9, ";"), .subset,1)) #just getting out transcript_ids
gff_genes$genes <- ids

Ors <- gff_genes[grep("Or", gff_genes$genes, perl=TRUE), ]#selecting only Ors
OBPs <- gff_genes[grep("OBP", gff_genes$genes, perl=TRUE), ]#selecting only OBPs

Full_chemosensory <- rbind(Ors, OBPs)
Full_chemosensory <- Full_chemosensory[
  with(Full_chemosensory, order(V1, V4)),]

#getting bedfiles
gff_forJoin_full <- gff_genes[,c(1,4,5)] #for full dataset
str(gff_forJoin_full)
names(gff_forJoin_full)<- c('chrom', 'chromStart','chromEnd') #prepping for bedtools
gff_forJoin_full$chrom <- as.numeric(as.factor(gff_forJoin_full$chrom))

write.table(gff_forJoin_full, file = "gff_forJoin_full.bed", col.names = FALSE, row.names = FALSE, sep = "\t", eol = "\n")


gff_forJoin_chemo <- Full_chemosensory [,c(1,4,5)] #Only Ors and OBPs
str(gff_forJoin_chemo)
names(gff_forJoin_chemo)<- c('chrom', 'chromStart','chromEnd') #prepping for bedtools
gff_forJoin_chemo$chrom <- as.numeric(as.factor(gff_forJoin_chemo$chrom))

write.table(gff_forJoin_chemo, file = "gff_forJoin_chemo.bed", col.names = FALSE, row.names = FALSE, sep = "\t", eol = "\n")

