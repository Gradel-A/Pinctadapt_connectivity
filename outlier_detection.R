library(vcfR)
library(dartR)
library(adegenet)
library(tidyverse)
library(radiator)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(StAMPP)
library(reshape)
library(gdata)
library(PopGenome)
library(ape)
library(igraph)
library(FSA)
library(pcadapt)
library(OutFLANK)
library(VennDiagram)

#load the vcf data in the genligth format with good population information
load("/Users/antoinegradel/Downloads/puce_filtered_hwe_genligth.Rdata")

#set seed:
set.seed(123)

#get a genind
obj.gi <- gl2gi(obj.glx)


####pcadapt all dataset ####
#Create correct format for pcadapt, outflanks and bayescan:
tidy<-genind2df(obj.gi)
dim(tidy)

#Extract position and chromosomes
geno <- tidy
position <- obj.glx@loc.names # Positions in bp

#Get the populations
population <- pop(obj.gi)

#Create a matrix to fit all the info for PCADAPT
ts.pcadapt <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))

ts.pcadapt[geno == "AA"] <- 0
ts.pcadapt[geno == c("AC", "CA")] <- 1
ts.pcadapt[geno == "CC"] <- 2

#Remove first column that corresponded to population in geno
ts.pcadapt<-ts.pcadapt[,2:ncol(ts.pcadapt)]

#load the data
input.pcadapt <- read.pcadapt(t(ts.pcadapt), type = "pcadapt")

save(input.pcadapt, file = "/Users/antoinegradel/Downloads/input.pcadapt.RData")

#doing the statistics, first using 10 PC axes
x <- pcadapt(input = input.pcadapt, K = 10)

pdf(file= "pcadapt_screeplot.pdf",width=11.7,heigh=8.3)
plot(x, option = "screeplot")
dev.off() #screeplot indicates 4 axes are enough to capture the variation

#doing the analyses with regards to the screeplot: K=4 axes
x <- pcadapt(input = input.pcadapt, K = 4)

#Look at the groupings according to geographical location (color)
coly<- c("purple","red","red","red","gold","gold",
         "light blue","light blue","light blue","light blue","light blue", "light blue", rep("green", 12), "orange")
obj.gi@pop<-factor(obj.gi@pop, levels = c("INDO","NHV-P", "NHV-S","UAH","SCI", "MOP",
                                          "MARS", "GMBW","PEI","MAT","RIK","MOR",
                                          "RAR","ANA","TAH","TAK","AHE","MAN","TKP","KAT","ARA","TEA","KAU","MOT",
                                          "RVV"))

pdf(file= "pcadapt_pcaplot.pdf",width=11.7,heigh=8.3)
plot(x, option = "scores", pop = obj.gi@pop,col=coly)
dev.off()

#represent outliers and save the manhattan plot
pdf(file= "pcadapt_manhattan.pdf",width=11.7,heigh=8.3)
plot(x , option = "manhattan")
dev.off()

#using q-value to detect the outliers
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
outliers.pcadapt <- which(qval < alpha)
length(outliers.pcadapt) #3040 avec 0.05; 1365 avec 0.001

#Look at distribution of p-values:
plot(x, option = "qqplot")


#Look at loadings on PC axes: here we don't see anything but it's because it's not placed according to chromosomes
par(mfrow = c(1, 4))
for (i in 1:4)
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

#Write the outliers:
PC.adapt<-as.character(outliers.pcadapt)
write.table(PC.adapt,"outliers.pcadapt.txt",quote=F,sep=",")


#### outflank all dataset ####
#Replace NA values with "9":
ts.outflanks<-ts.pcadapt
ts.outflanks[is.na(ts.outflanks)] <- 9
colnames(ts.outflanks)<-seq(1:ncol(ts.outflanks))
pop<-obj.gi@pop

#Calculate baseline Fst
my_fst <- MakeDiploidFSTMat(ts.outflanks, locusNames = colnames(ts.outflanks), popNames = pop)

hist(my_fst$FST,breaks=50)

#determine whether some SNPs are statistical outliers
OF <- OutFLANK(my_fst,LeftTrimFraction=0.05,RightTrimFraction=0.05,NumberOfSamples=25,qthreshold=0.05)


OutFLANKResultsPlotter(OF,withOutliers=T,
                       NoCorr=T,binwidth=0.005,
                       Zoom=F,RightZoomFraction=0.05,titletext=NULL)

#which SNPs are statistical outliers?
P1 <- pOutlierFinderChiSqNoCorr(my_fst,Fstbar=OF$FSTNoCorrbar,
                                dfInferred=OF$dfInferred,qthreshold=0.01)
outliers.outflanks <- P1$OutlierFlag==TRUE #which of the SNPs are outliers?
table(outliers.outflanks) #162 outliers


plot(P1$LocusName,P1$FST,xlab="Position",ylab="FST",col=rgb(0,0,0,alpha=0.1))
points(P1$LocusName[outliers.outflanks],P1$FST[outliers.outflanks],col="magenta")

#Write outliers:
OutFlanks<-P1$LocusName[outliers.outflanks][1:155] ##attention 162 a adapté
write.table(OutFlanks,"outliers.outflanks.txt",quote=F,sep=",")


#### doing the same appraoch with only Society, Tuamotu, Austral and Gambier ####
#### sampling localities ####
obj.glx.pol <- gl.drop.pop(obj.glx, c("INDO","NHV-P", "NHV-S","UAH"))
obj.gi.pol <- gl2gi(obj.glx.pol)

####pcadapt polynesia dataset ####
#Create correct format for pcadapt, outflanks and bayescan:
tidy<-genind2df(obj.gi.pol)
dim(tidy)

#Extract position and chromosomes
geno <- tidy
position <- obj.glx@loc.names # Positions in bp

#Get the populations
population <- pop(obj.gi)

#Create a matrix to fit all the info for PCADAPT
ts.pcadapt <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))

ts.pcadapt[geno == "AA"] <- 0
ts.pcadapt[geno == c("AC", "CA")] <- 1
ts.pcadapt[geno == "CC"] <- 2

#Remove first column that corresponded to population in geno
ts.pcadapt<-ts.pcadapt[,2:ncol(ts.pcadapt)]

#load the data
input.pcadapt <- read.pcadapt(t(ts.pcadapt), type = "pcadapt")

save(input.pcadapt, file = "/Users/antoinegradel/Downloads/input.pcadapt.pol.RData")

#doing the statistics, first using 10 PC axes
x <- pcadapt(input = input.pcadapt, K = 10)

pdf(file= "pcadapt_screeplot.pdf",width=11.7,heigh=8.3)
plot(x, option = "screeplot")
dev.off() #screeplot indicates 4 axes are enough to capture the variation

#doing the analyses with regards to the screeplot: K=2 axes
x <- pcadapt(input = input.pcadapt, K = 2)

#Look at the groupings according to geographical location (color)
coly<- c("gold","gold",
         "light blue","light blue","light blue","light blue","light blue", "light blue", rep("green", 12), "orange")
obj.gi@pop<-factor(obj.gi@pop, levels = c("SCI", "MOP",
                                          "MARS", "GMBW","PEI","MAT","RIK","MOR",
                                          "RAR","ANA","TAH","TAK","AHE","MAN","TKP","KAT","ARA","TEA","KAU","MOT",
                                          "RVV"))

pdf(file= "pcadapt_pcaplot.pdf",width=11.7,heigh=8.3)
plot(x, option = "scores", pop = obj.gi@pop,col=coly)
dev.off()

#represent outliers and save the manhattan plot
pdf(file= "pcadapt_manhattan.pdf",width=11.7,heigh=8.3)
plot(x , option = "manhattan")
dev.off()

#using q-value to detect the outliers
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
outliers.pcadapt <- which(qval < alpha)
length(outliers.pcadapt) #3040 avec 0.05; 1365 avec 0.001

#Look at distribution of p-values:
plot(x, option = "qqplot")


#Look at loadings on PC axes: here we don't see anything but it's because it's not placed according to chromosomes
par(mfrow = c(1, 4))
for (i in 1:4)
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

#Write the outliers:
PC.adapt<-as.character(outliers.pcadapt)
write.table(PC.adapt,"outliers.pcadapt.txt",quote=F,sep=",")


#### outflank polynesia dataset ####
#Replace NA values with "9":
ts.outflanks<-ts.pcadapt
ts.outflanks[is.na(ts.outflanks)] <- 9
colnames(ts.outflanks)<-seq(1:ncol(ts.outflanks))
pop<-obj.gi@pop

#Calculate baseline Fst
my_fst <- MakeDiploidFSTMat(ts.outflanks, locusNames = colnames(ts.outflanks), popNames = pop)

hist(my_fst$FST,breaks=50)

#determine whether some SNPs are statistical outliers
OF <- OutFLANK(my_fst,LeftTrimFraction=0.05,RightTrimFraction=0.05,NumberOfSamples=25,qthreshold=0.05)


OutFLANKResultsPlotter(OF,withOutliers=T,
                       NoCorr=T,binwidth=0.005,
                       Zoom=F,RightZoomFraction=0.05,titletext=NULL)

#which SNPs are statistical outliers?
P1 <- pOutlierFinderChiSqNoCorr(my_fst,Fstbar=OF$FSTNoCorrbar,
                                dfInferred=OF$dfInferred,qthreshold=0.01)
outliers.outflanks <- P1$OutlierFlag==TRUE #which of the SNPs are outliers?
table(outliers.outflanks) #162 outliers


plot(P1$LocusName,P1$FST,xlab="Position",ylab="FST",col=rgb(0,0,0,alpha=0.1))
points(P1$LocusName[outliers.outflanks],P1$FST[outliers.outflanks],col="magenta")

#Write outliers:
OutFlanks<-P1$LocusName[outliers.outflanks][1:155] ##attention 162 a adapté
write.table(OutFlanks,"outliers.outflanks.pol.txt",quote=F,sep=",")












