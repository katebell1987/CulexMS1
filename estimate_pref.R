################################################################################
##            Script to estimate multi-day preference for                    ##
##                       for Human 1 and 2                                    ##
################################################################################

setwd("~/Dropbox/UMD_Research/MS1_prefGBS/")

#### Install and load packages ####
#install.packages("bayespref")
#install.packages("MCMCpack")
#install.packages("ggpubr")
library(MCMCpack)
library(bayespref)
library(ggpubr)
library(ggplot2)

#### Read in human2 Data ####
data2<-read.csv("~/Dropbox/UMD_Research/Mapping_Exp/PrefTested.csv",header=T)
dat_mod2<-data2[,c(5:6,8)] # eva - 1, cal = 2, north = 3

#### Read in human1 data ####
data1<-read.csv("~/Dropbox/UMD_Research/Preference/RawData/Anna_fullDAT.csv")
data1<-as.data.frame(data1)

#### Organize human1 data ####
cal1<-subset(data1,strain=="cal1") # subset just to keep populations of interest
eva1<-subset(data1,strain=="eva")
nor1<-subset(data1,strain=="nor")   

cal1$mosq.id<-as.character(cal1$mosq.id) # Fix factor levels for subsets
IDsC<-as.character(unique(cal1$mosq.id))
levels(cal1$mosq.id)<-IDsC

eva1$mosq.id<-as.character(eva1$mosq.id)
IDsE<-as.character(unique(eva1$mosq.id))
levels(eva1$mosq.id)<-IDsE

nor1$mosq.id<-as.character(nor1$mosq.id)
IDsN<-as.character(unique(nor1$mosq.id))
levels(nor1$mosq.id)<-IDsN

response_cal1<-table(cal1$mosq.id,cal1$host) # Get frequency of response to host per individual
response_cal1<-as.data.frame.matrix(response_cal1)

response_eva1<-table(eva1$mosq.id,eva1$host)
response_eva1<-as.data.frame.matrix(response_eva1)

response_nor1<-table(nor1$mosq.id,nor1$host)
response_nor1<-as.data.frame.matrix(response_nor1)

popID1<-c(rep(2,length(response_cal1$h)),rep(1,length(response_eva1$h)),rep(3,length(response_nor1$h))) # add population ID

dat1<-rbind(response_cal1,response_eva1,response_nor1) # Bind together then add ID
dat_id1<-cbind(popID1,dat1)

# Remove individuals that didn't respond
fin_dat1<-subset(dat_id1,c>=1 | h >=1) ## cal 1 = 2, eva = 1, northfield = 3

########################### Cal ################################################
#### Combine data sets for Cal for human 1 and 2 ####
cal1_1<-fin_dat1[fin_dat1$popID1==2,1:3]
cal1_2<-dat_mod2[dat_mod2$PopID==2,]
c_cal<-c(cal1_1$c,cal1_2$C)
h_cal<-c(cal1_1$h,cal1_2$H)
HostID_cal<-c(rep(1,length(cal1_1$popID)),rep(2,length(cal1_2$PopID)))
comb_resp_cal<-cbind(HostID_cal,c_cal,h_cal) 

##### Run constrained model ####
calConstrained<-bayesPref(pData = comb_resp_cal, mcmcL = 50000, pops = TRUE, constrainP = c(1,1),dicburn = 2000) # dic =  -56.17882 

#### Check mixing ####
plot(calConstrained[[1]]$PopPref[1,2001:30000],xlab="MCMC step", ylab="PopPref")
effectiveSize(calConstrained[[1]]$PopPref[1,2001:30000]) #77.958

#### Get stats ####
cal_poppref<-credibleIntervals(prefres=c(calConstrained[[1]]),burn=2000,interval= 0.95)
cal_poppref$PopPref ### Population level preference for C [1] and H [2] rows, columns are the 2.5%, median, and 97.5% credible intervals. Constrained model

############################## Evanston ########################################
#### Combine data sets for Evanston for human 1 and 2 ####
Eva_1<-fin_dat1[fin_dat1$popID1==1,1:3]
Eva_2<-dat_mod2[dat_mod2$PopID==1,]
eva_c<-c(Eva_1$c,Eva_2$C)
eva_h<-c(Eva_1$h,Eva_2$H)
HostID_eva<-c(rep(1,length(Eva_1$popID1)),rep(2,length(Eva_2$PopID)))
comb_resp_Eva<-cbind(HostID_eva,eva_c,eva_h)

#### Run the model ####
evaConstrained<-bayesPref(pData = comb_resp_Eva, mcmcL = 50000, pops = TRUE, dicburn = 2000, constrainP = c(1,1)) ### dic =  -162.6227 

#### Check mixing ####
plot(evaConstrained[[1]]$PopPref[1,2001:50000],xlab="MCMC step", ylab="PopPref")
effectiveSize(evaConstrained[[1]]$PopPref[1,1001:50000])

#### Get stats ####
eva_poppref<-credibleIntervals(prefres=c(evaConstrained[[1]]),burn=2000,interval= 0.95)
eva_poppref$PopPref ### Population level preference for C [1] and H [2] rows, columns are the 2.5%, median, and 97.5% credible intervals. Constrained model

############################### North ##########################################
#### Combine data sets for Northfield for human 1 and 2 ####
Nor_1<-fin_dat1[fin_dat1$popID1==3,1:3]
Nor_2<-dat_mod2[dat_mod2$PopID==3,]
nor_c<-c(Nor_1$c,Nor_2$C)
nor_h<-c(Nor_1$h,Nor_2$H)
HostID_nor<-c(rep(1,length(Nor_1$popID1)),rep(2,length(Nor_2$PopID)))
comb_resp_Nor<-cbind(HostID_nor,nor_c,nor_h)

#### Run the model ####
norConstrained<-bayesPref(pData = comb_resp_Nor, mcmcL = 50000, pops = TRUE, dicburn = 2000, constrainP = c(1,1)) #dic =  -405.2342

#### Check mixing ####
plot(norConstrained[[1]]$PopPref[1,2001:50000],xlab="MCMC step", ylab="PopPref")
effectiveSize(norConstrained[[1]]$PopPref[1,2001:50000])

#### Get stats ####
north_poppref<-credibleIntervals(prefres=norConstrained[[1]],burn=2000,interval= 0.95)
north_poppref$PopPref


##### Make multi panel plots #####
burn=2000
colors_ind<-c("#bdbdbd40","#00000040")
colors_pop<-c("#bdbdbd","#000000")

pdf("~/Dropbox/UMD_Research/MS1_prefGBS/Figures/poppref_MDtesting.pdf",width=8,height=10)
par(mfrow=c(3,1))

##### cal ####
mcmc_cal<- dim(calConstrained[[1]]$PopPref)[2]
nCat_cal<- dim(calConstrained[[1]]$PopPref)[1]
nInd_cal<- dim(calConstrained[[1]]$IndPref)[1]

plot(0:1, 0:1, type = "n", ylim = c(0, 20), xlim = c(0,1), xlab = "Preference",ylab = "Probability Density", cex.lab = 1.2,cex.axis = 1.1,main="Below-ground 1 (BG1)")
for (i in 1:nInd_cal) {
  for (j in 1:nCat_cal) {
    lines(density(calConstrained[[1]]$IndPref[i, j, (burn + 1):mcmc_cal]), 
          col = colors_ind[j], lty = 2)
  }}
for (j in 1:nCat_cal) {
  lines(density(calConstrained[[1]]$PopPref[j, (burn + 1):mcmc_cal],adjust = 2),lty = 1, lwd =3,col=colors_pop[j])
}
legend(0.44,20, c("Avian","Human"), col = colors_pop, lty = c(rep(1, nCat_cal), rep(2, nCat_cal)), lwd = c(rep(2, nCat_cal)))


#### eva ####
mcmc_eva<- dim(evaConstrained[[1]]$PopPref)[2]
nCat_eva<- dim(evaConstrained[[1]]$PopPref)[1]
nInd_eva<- dim(evaConstrained[[1]]$IndPref)[1]

plot(0:1, 0:1, type = "n", ylim = c(0, 20), xlim = c(0,1), xlab = "Preference",ylab = "Probability Density", cex.lab = 1.2,cex.axis = 1.1,main="Above-ground 2 (AG2)")
for (i in 1:nInd_eva) {
  for (j in 1:nCat_eva) {
    lines(density(evaConstrained[[1]]$IndPref[i, j, (burn + 1):mcmc_eva]), 
          col = colors[j], lty = 2)
  }}
for (j in 1:nCat_eva) {
  lines(density(evaConstrained[[1]]$PopPref[j, (burn + 1):mcmc_eva],adjust = 2),lty = 1, lwd =3,col=colors[j]) 
}
legend(0.44,20, c("Avian","Human"), col = colors, lty = c(rep(1, nCat_eva), rep(2, nCat_eva)), lwd = c(rep(2, nCat_eva)))

#### north ####
mcmc_nor <- dim(norConstrained[[1]]$PopPref)[2]
nCat_nor <- dim(norConstrained[[1]]$PopPref)[1]
nInd_nor<- dim(norConstrained[[1]]$IndPref)[1]

plot(0:1, 0:1, type = "n", ylim = c(0, 20), xlim = c(0,1), xlab = "Preference",ylab = "Probability Density", cex.lab = 1.2,cex.axis = 1.1,main="Above-ground 3 (AG3)")
for (i in 1:nInd_nor) {
  for (j in 1:nCat_nor) {
    lines(density(norConstrained[[1]]$IndPref[i, j, (burn + 1):mcmc_nor]), 
          col = colors_ind[j], lty = 2)
  }}
for (j in 1:nCat_nor) {
  lines(density(norConstrained[[1]]$PopPref[j, (burn+1):mcmc],adjust = 2),lty = 1, lwd =3,col=colors_pop[j]) # Cal1
}
legend(0.44,20, c("Avian","Human"), col = colors_pop, lty = c(rep(1, nCat_nor), rep(2, nCat_nor)), lwd = c(rep(2, nCat_nor)))

dev.off()

save.image(file="~/Dropbox/UMD_Research/Preference/RawData/pop_pref_MDtesting.RData")

load("~/Dropbox/UMD_Research/Preference/RawData/pop_pref_MDtesting.RData")
pdf("~/Dropbox/UMD_Research/MS1_prefGBS/Figures/poppref_MDtesting_UPDATED.pdf",width=8,height=10)
par(mfrow=c(3,1))

##### cal ####
plot(0:1, 0:1, type = "n", ylim = c(0, 20), xlim = c(0,1), xlab = "Preference",ylab = "Probability Density", cex.lab = 1.2,cex.axis = 1.1,main="Below-ground 1 (BG1)")
# for (i in 1:nInd_cal) {
#   for (j in 1:nCat_cal) {
#     lines(density(calConstrained[[1]]$IndPref[i, j, (burn + 1):mcmc_cal]), 
#           col = colors_ind[j], lty = 2)
#   }}
for (j in 1:nCat_cal) {
  lines(density(calConstrained[[1]]$PopPref[j, (burn + 1):mcmc_cal],adjust = 2),lty = 1, lwd =3,col=colors_pop[j])
}
legend(0.44,20, c("Avian","Human"), col = colors_pop, lty = c(rep(1, nCat_cal), rep(2, nCat_cal)), lwd = c(rep(2, nCat_cal)))


#### eva ####
plot(0:1, 0:1, type = "n", ylim = c(0, 20), xlim = c(0,1), xlab = "Preference",ylab = "Probability Density", cex.lab = 1.2,cex.axis = 1.1,main="Above-ground 2 (AG2)")
# for (i in 1:nInd_eva) {
#   for (j in 1:nCat_eva) {
#     lines(density(evaConstrained[[1]]$IndPref[i, j, (burn + 1):mcmc_eva]), 
#           col = colors[j], lty = 2)
#   }}
for (j in 1:nCat_eva) {
  lines(density(evaConstrained[[1]]$PopPref[j, (burn + 1):mcmc_eva],adjust = 2),lty = 1, lwd =3,col=colors[j]) 
}
legend(0.44,20, c("Avian","Human"), col = colors, lty = c(rep(1, nCat_eva), rep(2, nCat_eva)), lwd = c(rep(2, nCat_eva)))

#### north ####
plot(0:1, 0:1, type = "n", ylim = c(0, 20), xlim = c(0,1), xlab = "Preference",ylab = "Probability Density", cex.lab = 1.2,cex.axis = 1.1,main="Above-ground 3 (AG3)")
# for (i in 1:nInd_nor) {
#   for (j in 1:nCat_nor) {
#     lines(density(norConstrained[[1]]$IndPref[i, j, (burn + 1):mcmc_nor]), 
#           col = colors_ind[j], lty = 2)
#   }}
for (j in 1:nCat_nor) {
  lines(density(norConstrained[[1]]$PopPref[j, (burn+1):mcmc],adjust = 2),lty = 1, lwd =3,col=colors_pop[j]) # Cal1
}
legend(0.44,20, c("Avian","Human"), col = colors_pop, lty = c(rep(1, nCat_nor), rep(2, nCat_nor)), lwd = c(rep(2, nCat_nor)))
dev.off()
