library(dartR)
library(adegenet)

load("/puce_filtered_hwe_genligth.Rdata")

outliers.summary <- read.table("/outliers.venn.csv", sep = "\t", header = TRUE)


atleast1 <- unique(c(outliers.summary$at_least2, outliers.summary$pcadapt, outliers.summary$outflank, outliers.summary$bayescan))
obj.glx.neutral.1 <- gl.drop.loc(obj.glx, loc.list = na.omit(obj.glx@loc.names[atleast1]))
obj.glx.outliers.1 <- gl.keep.loc(obj.glx, loc.list = na.omit(obj.glx@loc.names[atleast1]))

obj.glx
obj.glx.neutral.1
obj.glx.outliers.1

atleast2 <- outliers.summary$at_least2
obj.glx.neutral.2 <- gl.drop.loc(obj.glx, loc.list = na.omit(obj.glx@loc.names[atleast2]))
obj.glx.outliers.2 <- gl.keep.loc(obj.glx, loc.list = na.omit(obj.glx@loc.names[atleast2]))

obj.glx.outliers.2

alltech <- outliers.summary$all
obj.glx.neutral.all <- gl.drop.loc(obj.glx, loc.list = na.omit(obj.glx@loc.names[alltech]))
obj.glx.outliers.all <- gl.keep.loc(obj.glx, loc.list = na.omit(obj.glx@loc.names[alltech]))

obj.glx.outliers.2

out.all.pcadapt <- c(na.omit(obj.glx@loc.names[outliers.summary$all]),
                     na.omit(obj.glx@loc.names[outliers.summary$bayescan_pcadapt]),
                     na.omit(obj.glx@loc.names[outliers.summary$pcadapt_outflank]),
                     na.omit(obj.glx@loc.names[outliers.summary$pcadapt]))

out.all.outflanck <- c(na.omit(obj.glx@loc.names[outliers.summary$all]),
                         na.omit(obj.glx@loc.names[outliers.summary$pcadapt_outflank]),
                         na.omit(obj.glx@loc.names[outliers.summary$outflank]))
                       
out.all.bayescan <- c(na.omit(obj.glx@loc.names[outliers.summary$all]),
                      na.omit(obj.glx@loc.names[outliers.summary$bayescan_pcadapt]),
                      na.omit(obj.glx@loc.names[outliers.summary$bayescan]))
metada_outliers <- cbind(neutral = obj.glx.neutral.1@loc.names,
    outlier_pcadapt = c(out.all.pcadapt, rep("", length(obj.glx.neutral.1@loc.names)-length(out.all.pcadapt))),
      outlier_outflanck = c(out.all.outflanck, rep("", length(obj.glx.neutral.1@loc.names)-length(out.all.outflanck))),
      outlier_bayescan = c(out.all.bayescan, rep("", length(obj.glx.neutral.1@loc.names)-length(out.all.bayescan))))

write.table(metada_outliers,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t",
            file = "metadata_outliers_all.txt")

tuamotus <- c("AHE", "ANA", "ARA", "KAT", "KAU", "MAN", "MOT", "RAR", "TAH", "TAK", "TEA", "TKP")
marquesas <- c("NHV-S", "NHV-P", "UAH")
gambier <- c("GMBW", "MOR", "MARS", "PEI", "MAT", "RIK")
society <- c("SCI", "MOP")
australe <- c("RVV")
indo <- c("INDO")

archipelago <- as.character(pop(obj.glx.neutral.1))
for (i in tuamotus) {
  archipelago[archipelago == i] <- "Tuamotus"}

for (i in gambier) {
  archipelago[archipelago == i] <- "Gambier"}

for (i in society) {
  archipelago[archipelago == i] <- "Society"}

for (i in australe) {
  archipelago[archipelago == i] <- "Austral"}

for (i in marquesas) {
  archipelago[archipelago == i] <- "Marquesas"
}

for (i in indo) {
  archipelago[archipelago == i] <- "Indonesia"
}


#### ACP ####

obj.gi.1n <- gl2gi(obj.glx.neutral.1)
obj.gi.1o <- gl2gi(obj.glx.outliers.1)


#neutre
Xn.1<- scaleGen(obj.gi.1n, NA.method ="mean")
pcan.1 <- dudi.pca(Xn.1, scannf = FALSE, nf = 2)
new<-cbind(pcan.1$li,archipelago)
str(new)
names(new)<-c("PC1","PC2","pop")

centroids <- aggregate(cbind(PC1, PC2)~pop,new,mean)
test <- merge(new, centroids, by = "pop")
colnames(test) <- c("pop", "PC1", "PC2", "PC1_ctr", "PC2_ctr")

coly <- funky(length(unique(archipelago)))
a <- ggplot(data=test,aes(x=PC1,y=PC2,col=pop)) + 
  geom_point(alpha=0.5,size=2) + 
  scale_color_manual(values = c("#64ADF5","#00BA38", "#F8766D", "#e0ac5e", "#9a72aB", "#c27e42"))+
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") +
  stat_ellipse() +
  geom_point(data=centroids,size=3, alpha = 1) +
  geom_segment(data = test, aes(x = PC1, y = PC2, xend = PC1_ctr, yend = PC2_ctr), alpha = 0.3)+
  geom_label_repel(data=centroids, label = centroids$pop, size = 5, )+
  labs(x="PC1 (1.72%)", y="PC2 (0.80%)")
plot(a)

pcan.1$eig[1]/ sum(pcan.1$eig)
pcan.1$eig[2]/ sum(pcan.1$eig)


#outliers
Xo.1<- scaleGen(obj.gi.1o, NA.method ="mean")
pcao.1 <- dudi.pca(Xo.1, scannf = FALSE, nf = 2)
new<-cbind(pcao.1$li,archipelago)
str(new)
names(new)<-c("PC1","PC2","pop")

centroids <- aggregate(cbind(PC1, PC2)~pop,new,mean)
test <- merge(new, centroids, by = "pop")
colnames(test) <- c("pop", "PC1", "PC2", "PC1_ctr", "PC2_ctr")

#coly <- funky(length(unique(archipelago)))
b <- ggplot(data=test,aes(x=PC1,y=PC2,col=pop)) + 
  geom_point(alpha=0.5,size=2) + 
  scale_color_manual(values = c("#64ADF5","#00BA38", "#F8766D", "#e0ac5e", "#9a72aB", "#c27e42"))+
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") +
  stat_ellipse() +
  geom_point(data=centroids,size=3, alpha = 1) +
  geom_segment(data = test, aes(x = PC1, y = PC2, xend = PC1_ctr, yend = PC2_ctr), alpha = 0.3)+
  geom_label_repel(data=centroids, label = centroids$pop, size = 5, )+
  labs(x="PC1 (9.96%)", y="PC2 (2.68%)")
plot(b)

#### stat de base & find clusters ####

#neutre
grp.1n <- find.clusters(obj.gi.1n, stat = "AIC", n.pca = 600, max.n.clust = 25)
table(grp.1n$grp, archipelago)

grp.1n.bic <- find.clusters(obj.gi.1n, stat = "BIC", n.pca = 600, max.n.clust = 25)
table(grp.1n.bic$grp, pop(obj.gi.1n))

grp.1n.wss <- find.clusters(obj.gi.1n, stat = "WSS", n.pca = 600, max.n.clust = 25)
table(grp.1n.wss$grp, pop(obj.gi.1o))

FST.1 <- gl.fst.pop(obj.glx.neutral.1, verbose = 0, nboots = 999)
het.1n <- gl.report.heterozygosity(obj.glx.neutral.1, plot.out = FALSE)

upperTriangle(fstn.1$Fsts) <- lowerTriangle(fstn.1$Fsts, byrow = TRUE)
upperTriangle(fstn.1$Pvalues) <- lowerTriangle(fstn.1$Pvalues, byrow = TRUE)


as.data.frame(fstn.1$Fsts) -> df
df$Colonne_ref <- colnames(df)
df$Colonne_ref <- factor(df$Colonne_ref, levels = ordre_ligne, ordered = TRUE)
df_organise <- df %>% arrange(Colonne_ref)
fst_pannel_neutre <- df_organise[, ordre_colonne]



obj.glx.neutral.1.arch <- obj.glx.neutral.1
pop(obj.glx.neutral.1.arch) <- archipelago
FST.1arch <- gl.fst.pop(obj.glx.neutral.1.arch, verbose = 0)
het.1narch <- gl.report.heterozygosity(obj.glx.neutral.1.arch, plot.out = FALSE)

#outliers
grp.1o <- find.clusters(obj.gi.1o, stat = "AIC", n.pca = 600, max.n.clust = 25)
table(grp.1o$grp, archipelago)

grp.1o.bic <- find.clusters(obj.gi.1o, stat = "BIC", n.pca = 600, max.n.clust = 25)
table(grp.1o.bic$grp, pop(obj.gi.1o))

grp.1o.wss <- find.clusters(obj.gi.1o, stat = "WSS", n.pca = 600, max.n.clust = 25)
table(grp.1o.wss$grp, pop(obj.gi.1o))

fsto.1 <- gl.fst.pop(obj.glx.outliers.1, verbose = 0, nboots = 999)
het.1o <- gl.report.heterozygosity(obj.glx.outliers.1, plot.out = FALSE)
cbind(het.1o$pop, het.1o$Ho, het.1o$HoSD, het.1o$He, het.1o$HeSD, het.1o$FIS)

upperTriangle(fsto.1$Fsts) <- lowerTriangle(fsto.1$Fsts, byrow = TRUE)
upperTriangle(fsto.1$Pvalues) <- lowerTriangle(fsto.1$Pvalues, byrow = TRUE)

save(fsto.1, file = "/fst_outlier.Rdata")

as.data.frame(fsto.1$Fsts) -> df
df$Colonne_ref <- colnames(df)
df$Colonne_ref <- factor(df$Colonne_ref, levels = ordre_ligne, ordered = TRUE)
df_organise <- df %>% arrange(Colonne_ref)
fst_pannel_outlier <- df_organise[, ordre_colonne]

obj.glx.outliers.1.arch <- obj.glx.outliers.1
pop(obj.glx.outliers.1.arch) <- archipelago
fsto.1arch <- gl.fst.pop(obj.glx.outliers.1.arch, verbose = 0, nboots = 999)
het.1oarch <- gl.report.heterozygosity(obj.glx.outliers.1.arch, plot.out = FALSE)
cbind(het.1oarch$pop, het.1oarch$Ho, het.1oarch$HoSD, het.1oarch$He, het.1oarch$HeSD, het.1oarch$FIS)

save(fsto.1, file = "/fst_outlier.Rdata")
save(fstn.1, file = "/fst_neutral.Rdata")
save(fsto.1arch, file = "/fst_outlier_arch.Rdata")
save(fstn.1arch, file = "/fst_neutral.arch.Rdata")

#### focus on french polynesia ####
obj.glx.pol <- gl.drop.pop(obj.glx, c(indo, marquesas))
obj.gi.pol <- gl2gi(obj.glx.pol)

#Create correct format for pcadapt, outflanks and bayescan:
tidy.pol<-genind2df(obj.gi.pol)
dim(tidy.pol)

#Extract position and chromosomes
geno.pol <- tidy.pol
position.pol<- obj.glx.pol@loc.names # Positions in bp

#Get the populations
population <- pop(obj.gi.pol)

#Create a matrix to fit all the info for PCADAPT
ts.pol.pcadapt <- matrix(NA, nrow = nrow(geno.pol), ncol = ncol(geno.pol))

ts.pol.pcadapt[geno.pol == "AA"] <- 0
ts.pol.pcadapt[geno.pol == c("AC", "CA")] <- 1
ts.pol.pcadapt[geno.pol == "CC"] <- 2
ts.pol.pcadapt<-ts.pol.pcadapt[,2:ncol(ts.pol.pcadapt)]

#load the data
input.pcadapt.pol <- read.pcadapt(t(ts.pol.pcadapt), type = "pcadapt")

save(input.pcadapt.pol, file = "/input.pcadapt.pol.RData")

#perform the analyses
x <- pcadapt(input = input.pcadapt.pol, K = 10)
plot(x, option = "screeplot") # 2 is enough

x.pol <- pcadapt(input = input.pcadapt.pol, K = 2)
colypol<- c("gold","gold",
         "light blue","light blue","light blue","light blue","light blue", "light blue", rep("green", 12), "orange")
obj.gi.pol@pop<-factor(obj.gi.pol@pop, levels = c("SCI", "MOP",
                                          "MARS", "GMBW","PEI","MAT","RIK","MOR",
                                          "RAR","ANA","TAH","TAK","AHE","MAN","TKP","KAT","ARA","TEA","KAU","MOT",
                                          "RVV"))

plot(x.pol, option = "scores", pop = obj.gi.pol@pop,col=colypol)
plot(x.pol , option = "manhattan")

qval <- qvalue(x.pol$pvalues)$qvalues
alpha <- 0.05
outliers.pcadapt.pol <- which(qval < alpha)
length(outliers.pcadapt.pol) # 0.5 105,  0.1 54, 0.001 22

plot(x.pol, option = "qqplot")
par(mfrow = c(1, 2))
for (i in 1:2)
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
par(mfrow = c(1, 1))

#Write the outliers:
PC.adapt.pol<-as.character(outliers.pcadapt.pol)
write.table(PC.adapt.pol,"/outliers.pcadapt.pol.txt",quote=F,sep=",")

##now look at outflanks
ts.pol.outflanks<-ts.pol.pcadapt
ts.pol.outflanks[is.na(ts.pol.outflanks)] <- 9
colnames(ts.pol.outflanks)<-seq(1:ncol(ts.pol.outflanks))
pop<-obj.gi.pol@pop

my_fst <- MakeDiploidFSTMat(ts.pol.outflanks, locusNames = colnames(ts.pol.outflanks), popNames = pop)

hist(my_fst$FST,breaks=50)

OF.pol <- OutFLANK(my_fst,LeftTrimFraction=0.05,RightTrimFraction=0.05,NumberOfSamples=25,qthreshold=0.05)


OutFLANKResultsPlotter(OF.pol,withOutliers=T,
                       NoCorr=T,binwidth=0.005,
                       Zoom=F,RightZoomFraction=0.05,titletext=NULL)

P1.pol <- pOutlierFinderChiSqNoCorr(my_fst,Fstbar=OF$FSTNoCorrbar,
                                dfInferred=OF$dfInferred,qthreshold=0.05)
outliers.outflanks.pol <- P1.pol$OutlierFlag==TRUE #which of the SNPs are outliers?
table(outliers.outflanks.pol) #2 outliers


plot(P1.pol$LocusName,P1.pol$FST,xlab="Position",ylab="FST",col=rgb(0,0,0,alpha=0.1))
points(P1.pol$LocusName[outliers.outflanks.pol],P1.pol$FST[outliers.outflanks.pol],col="magenta")

OutFlanks.pol<-P1.pol$LocusName[outliers.outflanks.pol][1:2] ##attention 162 a adapté
write.table(OutFlanks,"outliers.outflanks.pol.txt",quote=F,sep=",")


####stats on polynesia centered dataset ####
out.bayescan.pol <- read.table("/Bayescan_pol_selected_locus_id.txt", sep = "\t", header = TRUE)
out.pcadapt.pol <- read.table("/article_1/outliers.pcadapt.pol.txt", header = TRUE, sep = ",")
out.outflank.pol <- read.table("/outliers.outflanks.pol.txt", header = TRUE, sep = ",")

obj.glx.pol <- gl.drop.pop(obj.glx, c(indo, marquesas))

atleast1.pol <- unique(c(out.outflank.pol$x, out.pcadapt.pol$x, out.bayescan.pol$index))

obj.glx.neutral.1.pol <- gl.drop.loc(obj.glx.pol, loc.list = na.omit(obj.glx.pol@loc.names[as.integer(atleast1.pol)]))
obj.glx.outliers.1.pol <- gl.keep.loc(obj.glx.pol, loc.list = na.omit(obj.glx.pol@loc.names[as.integer(atleast1.pol)]))

pol.pcadapt <- na.omit(obj.glx.pol@loc.names[as.integer(out.pcadapt.pol$x)])
pol.outflank <- na.omit(obj.glx.pol@loc.names[as.integer(out.outflank.pol$x)])
pol.bayescan <- na.omit(obj.glx.pol@loc.names[as.integer(out.bayescan.pol$index)])

metadata_outliers_pol <- cbind(neutral = obj.glx.neutral.1.pol@loc.names,
                               outliers_pcadapt = c(pol.pcadapt, rep("", length(obj.glx.neutral.1.pol@loc.names)-length(pol.pcadapt))),
                               outlier_outflank = c(pol.outflank, rep("", length(obj.glx.neutral.1.pol@loc.names)-length(pol.outflank))),
                               outlier_bayescan = c(pol.bayescan, rep("", length(obj.glx.neutral.1.pol@loc.names)-length(pol.bayescan)))
                               )

write.table(metada_outliers,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t",
            file = "/metadata_outliers_poll.txt")



#### ACP ####

obj.gi.1n.pol <- gl2gi(obj.glx.neutral.1.pol)
obj.gi.1o.pol <- gl2gi(obj.glx.outliers.1.pol)


tuamotus <- c("AHE", "ANA", "ARA", "KAT", "KAU", "MAN", "MOT", "RAR", "TAH", "TAK", "TEA", "TKP")
gambier <- c("GMBW", "MOR", "MARS", "PEI", "MAT", "RIK")
society <- c("SCI", "MOP")
australe <- c("RVV")

archipelago.pol <- as.character(pop(obj.gi.1n.pol))
for (i in tuamotus) {
  archipelago.pol[archipelago.pol == i] <- "Tuamotus"}

for (i in gambier) {
  archipelago.pol[archipelago.pol == i] <- "Gambier"}

for (i in society) {
  archipelago.pol[archipelago.pol == i] <- "Society"}

for (i in australe) {
  archipelago.pol[archipelago.pol == i] <- "Austral"}

#stats
#neutre
Xn.1.pol<- scaleGen(obj.gi.1n.pol, NA.method ="mean")
pcan.1.pol <- dudi.pca(Xn.1.pol, scannf = FALSE, nf = 2)
newn<-cbind(pcan.1.pol$li,archipelago.pol)
str(newn)
names(newn)<-c("PC1","PC2","pop")

centroidsn <- aggregate(cbind(PC1, PC2)~pop,newn,mean)
testn <- merge(newn, centroidsn, by = "pop")
colnames(testn) <- c("pop", "PC1", "PC2", "PC1_ctr", "PC2_ctr")

#coly <- funky(length(unique(archipelago.pol)))
c <- ggplot(data=testn,aes(x=PC1,y=PC2,col=pop)) + 
  geom_point(alpha=0.5,size=2) + 
  scale_color_manual(values = c("#64ADF5","#00BA38", "#9a72aB", "#c27e42"))+
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") +
  stat_ellipse() +
  geom_point(data=centroidsn,size=3, alpha = 1) +
  geom_segment(data = testn, aes(x = PC1, y = PC2, xend = PC1_ctr, yend = PC2_ctr), alpha = 0.3)+
  geom_label_repel(data=centroidsn, label = centroidsn$pop, size = 5, )+
  labs(x="PC1 (0.28%)", y="PC2 (0.25%)")
plot(c)

pcan.1.pol$eig[1]/ sum(pcan.1.pol$eig)
pcan.1.pol$eig[2]/ sum(pcan.1.pol$eig)

#outliers
Xo.1.pol<- scaleGen(obj.gi.1o.pol, NA.method ="mean")
pcao.1.pol <- dudi.pca(Xo.1.pol, scannf = FALSE, nf = 2)
newo<-cbind(pcao.1.pol$li,archipelago.pol)
str(newo)
names(newo)<-c("PC1","PC2","pop")

centroidso <- aggregate(cbind(PC1, PC2)~pop,newo,mean)
testo <- merge(newo, centroidso, by = "pop")
colnames(testo) <- c("pop", "PC1", "PC2", "PC1_ctr", "PC2_ctr")

#coly <- funky(length(unique(archipelago.pol)))
d <- ggplot(data=testo,aes(x=PC1,y=PC2,col=pop)) + 
  geom_point(alpha=0.5,size=2) + 
  scale_color_manual(values = c("#64ADF5","#00BA38", "#9a72aB", "#c27e42"))+
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") +
  stat_ellipse() +
  geom_point(data=centroidso,size=3, alpha = 1) +
  geom_segment(data = testo, aes(x = PC1, y = PC2, xend = PC1_ctr, yend = PC2_ctr), alpha = 0.3)+
  geom_label_repel(data=centroidso, label = centroidso$pop, size = 5, )+
  labs(x="PC1 (3.63%)", y="PC2 (2.53%)")
plot(d)

pcao.1.pol$eig[1]/ sum(pcao.1.pol$eig)
pcao.1.pol$eig[2]/ sum(pcao.1.pol$eig)

ggarrange(a,b,c,d, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), font.label = list(size = 16))
#### stat de base et find clusters ####


#neutre


grp.1n.pol <- find.clusters(obj.gi.1n.pol, stat = "AIC", n.pca = 400, max.n.clust = 21)
table(grp.1n.pol$grp, pop(obj.gi.1n.pol))

grp.1n.pol.bic <- find.clusters(obj.gi.1n.pol, stat = "BIC", n.pca = 400, max.n.clust = 21)
table(grp.1n.pol.bic$grp, pop(obj.gi.1n.pol))

grp.1n.pol.wss <- find.clusters(obj.gi.1n.pol, stat = "WSS", n.pca = 400, max.n.clust = 21)
# avec criterion diffNgroup il choisi n=3
table(grp.1n.pol.wss$grp, pop(obj.gi.1n.pol))

fst.1.pol <- gl.fst.pop(obj.glx.neutral.1.pol, verbose = 0, nboots = 999)
upperTriangle(fst.1.pol$Pvalues) <- lowerTriangle(fst.1.pol$Pvalues, byrow = TRUE)
het.1n.pol <- gl.report.heterozygosity(obj.glx.neutral.1.pol, plot.out = FALSE)

obj.glx.neutral.1.pol.arch <- obj.glx.neutral.1.pol
pop(obj.glx.neutral.1.pol.arch) <- archipelago.pol
FST.1arch.pol <- gl.fst.pop(obj.glx.neutral.1.pol.arch, verbose = 0)
het.1narch.pol <- gl.report.heterozygosity(obj.glx.neutral.1.pol.arch, plot.out = FALSE)

#outliers
grp.1o.pol <- find.clusters(obj.gi.1o.pol, stat = "AIC", n.pca = 400, max.n.clust = 21)
table(grp.1o.pol$grp, pop(obj.gi.1o.pol))

grp.1o.pol.bic <- find.clusters(obj.gi.1o.pol, stat = "BIC", n.pca = 400, max.n.clust = 21)
table(grp.1o.pol.bic$grp, pop(obj.gi.1o.pol))

grp.1o.pol.wss <- find.clusters(obj.gi.1o.pol, stat = "WSS", n.pca = 400, max.n.clust = 21)
table(grp.1o.pol.wss$grp, pop(obj.gi.1o.pol))

fsto.1.pol <- gl.fst.pop(obj.glx.outliers.1.pol, verbose = 0, nboots = 999)
het.1o.pol <- gl.report.heterozygosity(obj.glx.outliers.1.pol, plot.out = FALSE)

obj.glx.outliers.1.pol.arch <- obj.glx.outliers.1.pol
pop(obj.glx.outliers.1.pol.arch) <- archipelago.pol
fsto.1arch.pol <- gl.fst.pop(obj.glx.outliers.1.pol.arch, verbose = 0)
het.1oarch.pol <- gl.report.heterozygosity(obj.glx.outliers.1.pol.arch, plot.out = FALSE)

#### preparation amova ####
#all
obj.glx.neutral.amova <- obj.glx.neutral.1.arch
obj.glx.outlier.amova <- obj.glx.outliers.1.arch

cluster <- pop(obj.glx.neutral.amova)
levels(cluster) <- c(levels(cluster), "Polynesia")
cluster[cluster == "Austral"] <- "Polynesia"
cluster[cluster == "Tuamotus"] <- "Polynesia"
cluster[cluster == "Gambier"] <- "Polynesia"
cluster[cluster == "Society"] <- "Polynesia"

pop(obj.glx.neutral.amova) <- cluster
pop(obj.glx.outlier.amova) <- cluster

save(obj.glx.neutral.amova, file = "~/Desktop/obj.glx.neutral.amova.RData")
save(obj.glx.outlier.amova, file = "~/Desktop/obj.glx.outlier.amova.RData")

strata.amova <- as.data.frame(cbind(archipelago_island=paste(cluster, pop(obj.glx.neutral.1),sep="_"),
                                 archipelago = as.character(cluster),
                                 island = as.character(pop(obj.glx.neutral.1))))

save(strata.amova, file = "~/Desktop/strata.amova.RData")

obj.genclone.neutral <- as.genclone(gl2gi(obj.glx.neutral.amova))
strata(obj.genclone.neutral) <- strata.amova
obj.genclone.neutral

obj.genclone.outlier <- as.genclone(gl2gi(obj.glx.outlier.amova))
strata(obj.genclone.outlier) <- strata.amova
obj.genclone.outlier


save(obj.genclone.neutral, file = "~/Desktop/obj.genclone.neutral.RData")
save(obj.genclone.outlier, file = "~/Desktop/obj.genclone.outlier.RData")

# polynesia centered
obj.glx.neutral.pol.amova <- gl.drop.pop(obj.glx.neutral.1.pol.arch, "Austral")
obj.glx.outlier.pol.amova <- gl.drop.pop(obj.glx.outliers.1.pol.arch, "Austral")
islands.pol.wo.austral <- gl.drop.pop(obj.glx.neutral.1.pol, "RVV")

save(obj.glx.neutral.pol.amova, file = "~/Desktop/obj.glx.neutral.pol.amova.RData")
save(obj.glx.outlier.pol.amova, file = "~/Desktop/obj.glx.outlier.pol.amova.RData")

strata.pol.amova <- as.data.frame(cbind(archipelago_island=paste(pop(obj.glx.neutral.pol.amova), pop(islands.pol.wo.austral),sep="_"),
                                    archipelago = as.character(pop(obj.glx.neutral.pol.amova)),
                                    island = as.character(pop(islands.pol.wo.austral))))

save(strata.pol.amova, file = "~/Desktop/strata.pol.amova.RData")

obj.genclone.neutral.pol <- as.genclone(gl2gi(obj.glx.neutral.pol.amova))
strata(obj.genclone.neutral.pol) <- strata.pol.amova
obj.genclone.neutral.pol

obj.genclone.outlier.pol <- as.genclone(gl2gi(obj.glx.outlier.pol.amova))
strata(obj.genclone.outlier.pol) <- strata.pol.amova
obj.genclone.outlier.pol

save(obj.genclone.neutral.pol, file = "~/Desktop/obj.genclone.neutral.pol.RData")
save(obj.genclone.outlier.pol, file = "~/Desktop/obj.genclone.outlier.pol.RData")










#### centrality plotting ####
centrality <- read.table("/centrality_data.csv", 
                         sep = ";", dec = ".", header = TRUE)
centrality$Archipelago <- as.factor(centrality$Archipelago)
centrality$sampled[centrality$sampled == ""] <- "No"
centrality$sampled <- as.factor(centrality$sampled)

str(centrality)

library(sf)
library(ggplot2)
library(ggrepel)

chemin_shapefile <- "~/Downloads/ne_10m_admin_0_scale_rank_minor_islands/ne_10m_admin_0_scale_rank_minor_islands.shp"

centrality$betw_norm <-   (centrality$betweenness - min(centrality$betweenness, na.rm = TRUE)) / (max(centrality$betweenness, na.rm = TRUE) - min(centrality$betweenness, na.rm = TRUE))
centrality$harm_norm <-   (centrality$harmonic - min(centrality$harmonic, na.rm = TRUE)) / (max(centrality$harmonic, na.rm = TRUE) - min(centrality$harmonic, na.rm = TRUE))
centrality$eig_norm <- (centrality$eigenvector - min(centrality$eigenvector, na.rm = TRUE)) / (max(centrality$eigenvector, na.rm = TRUE) - min(centrality$eigenvector, na.rm = TRUE))


# Créez un objet sf directement à partir du DataFrame
sf_echantillons <- st_as_sf(
  centrality,
  coords = c("Longitude", "Latitude"),  # Spécifiez les colonnes pour les coordonnées
  crs = 4326  # Définissez le système de coordonnées (WGS84 ici)
)

# Chargez le shapefile
world <- st_read(chemin_shapefile)

# Trouvez les limites de la zone où se trouvent vos points
bbox <- st_bbox(sf_echantillons)


ggplot() +
  geom_sf(data = world, fill = "white", color = "black") +
  geom_sf(data = sf_echantillons, aes(color = betweenness, size = betweenness)) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +
  labs(x = "Longitude", y = "Latitude") +  # Ajout des noms personnalisés pour les axes
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),  # Ajustez la taille du texte de l'axe x
        axis.title.y = element_text(size = 16),
        legend.position = "right")+
  scale_color_continuous(type = "viridis", name = "Betweeness
Centrality")+
  scale_size_continuous(guide = "none")

top_iles <- sf_echantillons %>%
  arrange(desc(betw_norm)) %>%
  slice_head(n = 3) %>%
  mutate(nudge_x = c(-5, 2, 5),  # Décalages spécifiques en x
         nudge_y = c(-4, -7, 5))  # Décalages spécifiques en y

top_ilesb <- top_iles %>%
  mutate(
    x = st_coordinates(.)[,1],  # Extraction des coordonnées x
    y = st_coordinates(.)[,2],
    nudge_x = c(-5, 2, 5),  # Décalages spécifiques en x
    nudge_y = c(-4, -7, 5)) 

betweeness_plot <- ggplot() +
  geom_sf(data = world, fill = "white", color = "black") +
  geom_sf(data = na.omit(sf_echantillons), aes(color = betw_norm), size = 3) +
  geom_sf(data = subset(sf_echantillons, sampled == "Yes"), 
          color = "white", size = 1.2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +
  labs(x = "Longitude", y = "Latitude") +  # Ajout des noms personnalisés pour les axes
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),  # Ajustez la taille du texte de l'axe x
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        legend.title = element_text(size = 14),
        title = element_text(size = 16))+
  ggtitle("B")+
  scale_color_continuous(type = "viridis", name = "Betweeness
Centrality")+
 scale_size_continuous(breaks = c(1.00,0.75,0.50,0.25,0.00))+
  geom_segment(
    data = top_ilesb,
    aes(x = x, y = y, xend = x + nudge_x, yend = y + nudge_y),
    color = "black",
    linewidth = 0.5)+
  geom_label(
    data = top_iles,
    aes(label = Island, geometry = geometry, x = ..x.. + top_iles$nudge_x, y = ..y.. + top_iles$nudge_y),
    stat = "sf_coordinates",# Décalage des étiquettes
    size = 5, fill = "white", color = "black", label.size = 0.5)
  
plot(betweeness_plot)  

  
harmonic_plot <- ggplot() +
  geom_sf(data = world, fill = "white", color = "black") +
  geom_sf(data = na.omit(sf_echantillons), aes(color = harm_norm), size = 3) +
  geom_sf(data = subset(sf_echantillons, sampled == "Yes"), 
          color = "white", size = 1.2) +
  geom_label_repel(data = na.omit(sf_echantillons), 
                 #  aes(geometry = geometry, label = ifelse(harm_norm > 0.95, Island, NA)), 
                  aes(geometry = geometry, label = ifelse(Island == "Katiu" |
                                                            Island == "Kauehi" |
                                                            Island == "Mangareva" |
                                                            harm_norm > 0.95 &
                                                            Island != "RavahereMarokau", Island, NA )),
                   stat = "sf_coordinates",  # Permet d'afficher les étiquettes correctement sur une carte
                   box.padding = 5, 
                   point.padding = 0.5,
                   segment.color = "black",
                   segment.size = 0.5,
                   size = 5, 
                   min.segment.length = 0,
                   max.overlaps = 10)+
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +
  labs(x = "Longitude", y = "Latitude") +  # Ajout des noms personnalisés pour les axes
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),  # Ajustez la taille du texte de l'axe x
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        legend.title = element_text(size = 14),
        title = element_text(size = 16))+
  scale_color_continuous(type = "viridis", name = "Harmonic
Centrality")+
  scale_size_continuous(guide = "none")+
ggtitle("A")
plot(harmonic_plot)

eigen_plot <- ggplot() +
  geom_sf(data = world, fill = "white", color = "black") +
  geom_sf(data = na.omit(sf_echantillons), aes(color = eig_norm), size =3) +
  geom_sf(data = subset(sf_echantillons, sampled == "Yes"), 
          color = "white", size = 1.2) +
  geom_label_repel(data = na.omit(sf_echantillons), 
                   #  aes(geometry = geometry, label = ifelse(harm_norm > 0.95, Island, NA)), 
                   aes(geometry = geometry, label = ifelse(Island == "Takume" |
                                                             Island == "Tahanea" |
                                                             Island == "Hikueru" |
                                                             Island == "Manihi" |
                                                             Island == "Hereheretue" , Island, NA )),
                   stat = "sf_coordinates",  # Permet d'afficher les étiquettes correctement sur une carte
                   box.padding = 5, 
                   point.padding = 0.5,
                   segment.color = "black",
                   segment.size = 0.5,
                   size = 5, 
                   min.segment.length = 0,
                   max.overlaps = 10)+
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +
  labs(x = "Longitude", y = "Latitude") +  # Ajout des noms personnalisés pour les axes
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),  # Ajustez la taille du texte de l'axe x
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        legend.title = element_text(size = 14),
        title = element_text(size = 16))+
  scale_color_continuous(type = "viridis", name = "Eigenvector
Centrality")+
  scale_size_continuous(guide = "none")+
  ggtitle("C")
plot(eigen_plot)

centrality$Archipelago <- as.character(centrality$Archipelago)

centrality$Archipelago[centrality$Archipelago == "society"] <- "Society"
centrality$Archipelago[centrality$Archipelago == "austral"] <- "Austral"
centrality$Archipelago[centrality$Archipelago == "marquesas"] <- "Marquesas"
centrality$Archipelago[centrality$Archipelago == "tuamotuWest"] <- "Tuamotu"
centrality$Archipelago[centrality$Archipelago == "tuamotuEast_gambier"] <- "Gambier"

combine_plot <- ggplot(data = na.omit(centrality))+
  geom_point(aes(x= harm_norm, y = betw_norm, color = Archipelago, size = eig_norm)) +
  geom_point(data = na.omit(centrality[centrality$sampled == "Yes",]), 
             aes(x= harm_norm, y = betw_norm), 
          color = "white", size = 1) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),  # Ajustez la taille du texte de l'axe x
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        title = element_text(size = 16),
        legend.text = element_text(size = 14))+
  ggtitle("D")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_size_continuous(name = "Eigenvector
Centrality", breaks = c(1.00,0.75,0.50,0.25,0.00))+
  xlab("Harmonic Centrality")+
  ylab("Betweeness Centrality")+
  scale_color_manual(values = c("#00BA38", "#9a72aB", "#c27e42"))#+
#  geom_label_repel(aes(x = harm_norm, y = betw_norm, label = ifelse(betw_norm > 0.75|
 #                                                                     Island == "Mangareva", Island, "")), 
  #                 box.padding = 2,  # Espace entre le point et l'étiquette
  #                 point.padding = 0,  # Espace entre le texte et la bulle
  #                 segment.color = "black",  # Couleur du trait reliant le point
  #                 size = 5,                # Taille du texte
  #                 min.segment.length = 0,
  #                 max.overlaps = 10)
plot(combine_plot)

ggarrange(harmonic_plot, betweeness_plot, eigen_plot, combine_plot, nrow = 2, ncol = 2)

#### IBD ####
library(vcfR)
library(dartR)
library(geosphere)
library(vegan)
library(tidyverse)

#matrix data generation
read.table("/gps_puce.txt", sep = "\t", header = T) -> data_gps
row.names(data_gps) <- data_gps$X
data_gps <- data_gps[,2:3]

tmp_total <- rep(0,25)
for (i in 1:25) {
  distGeo(data_gps, data_gps[i,]) -> tmp
  tmp_total <- cbind(tmp_total, tmp)
}
geo_df <- as.data.frame(tmp_total[1:25, 2:26])
rownames(geo_df) <- colnames(geo_df) <- rownames(data_gps)
geo_df
geo_df_order <- geo_df[rownames(IBD_mat_n), colnames(IBD_mat_n)]
geo_df_order


geo_df_order_pol_marq <- geo_df_order[rownames(geo_df_order) != "INDO",
                                      colnames(geo_df_order)!= "INDO"]
identical(rownames(IBD_mat_n_pol_marq), rownames(geo_df_order_pol_marq))
identical(colnames(IBD_mat_n_pol_marq), colnames(geo_df_order_pol_marq))


geo_df_order_pol <- geo_df_order_pol_marq[rownames(geo_df_order_pol_marq) != "UAH",
                                          colnames(geo_df_order_pol_marq) != "UAH"]
geo_df_order_pol <- geo_df_order_pol[rownames(geo_df_order_pol) != "NHV-P",
                                     colnames(geo_df_order_pol) != "NHV-P"]
geo_df_order_pol <- geo_df_order_pol[rownames(geo_df_order_pol) != "NHV-S",
                                     colnames(geo_df_order_pol) != "NHV-S"]

identical(rownames(IBD_mat_n_pol), rownames(geo_df_order_pol))
identical(colnames(IBD_mat_n_pol), colnames(geo_df_order_pol))


fst_mat_n <- fstn.1[["Fsts"]]
IBD_mat_n <- fst_mat_n/(1-fst_mat_n)
IBD_mat_n_pol_marq <- IBD_mat_n[rownames(IBD_mat_n) != "INDO", colnames(IBD_mat_n) != "INDO"]

fst_mat_n_pol <- fst.1.pol[["Fsts"]]
IBD_mat_n_pol <- fst_mat_n_pol/(1-fst_mat_n_pol)

IBD_mat_n_tuam <- IBD_mat_n_pol[rownames(IBD_mat_n_pol) != "RIK", 
                            colnames(IBD_mat_n_pol) != "RIK"]
IBD_mat_n_tuam <- IBD_mat_n_tuam[rownames(IBD_mat_n_tuam) != "PEI", 
                                colnames(IBD_mat_n_tuam) != "PEI"]
IBD_mat_n_tuam <- IBD_mat_n_tuam[rownames(IBD_mat_n_tuam) != "MAT", 
                                 colnames(IBD_mat_n_tuam) != "MAT"]
IBD_mat_n_tuam <- IBD_mat_n_tuam[rownames(IBD_mat_n_tuam) != "GMBW", 
                                 colnames(IBD_mat_n_tuam) != "GMBW"]
IBD_mat_n_tuam <- IBD_mat_n_tuam[rownames(IBD_mat_n_tuam) != "MOR", 
                                 colnames(IBD_mat_n_tuam) != "MOR"]
IBD_mat_n_tuam <- IBD_mat_n_tuam[rownames(IBD_mat_n_tuam) != "MARS", 
                                 colnames(IBD_mat_n_tuam) != "MARS"]
IBD_mat_n_tuam <- IBD_mat_n_tuam[rownames(IBD_mat_n_tuam) != "RVV", 
                                 colnames(IBD_mat_n_tuam) != "RVV"]
IBD_mat_n_tuam <- IBD_mat_n_tuam[rownames(IBD_mat_n_tuam) != "MOP", 
                                 colnames(IBD_mat_n_tuam) != "MOP"]
IBD_mat_n_tuam <- IBD_mat_n_tuam[rownames(IBD_mat_n_tuam) != "SCI", 
                                 colnames(IBD_mat_n_tuam) != "SCI"]



geo_df_order_TUAM <- geo_df[rownames(IBD_mat_n_tuam), colnames(IBD_mat_n_tuam)]

IBD_mat_n_gamb <- IBD_mat_n_pol[(gambier),(gambier)]
geo_df_order_GAMB <- geo_df[rownames(IBD_mat_n_gamb), colnames(IBD_mat_n_gamb)]


#statistic all, 0.480 p= 0.0033
mantel(IBD_mat_n, geo_df_order, permutations = 9999)

#polynesia 0.345 p = 0.0018
mantel(IBD_mat_n_pol_marq, geo_df_order_pol_marq, permutations = 9999)

#polynesian cluster 0.506 p< 0.001
mantel(IBD_mat_n_pol, geo_df_order_pol, permutations = 9999)

#tuamotu 0.114 p = 0.2708
mantel(IBD_mat_n_tuam, geo_df_order_TUAM, permutations = 9999)

#gambier 0.8045 p=0.0653
mantel(IBD_mat_n_gamb, geo_df_order_GAMB, permutations = 9999)


#### admixture ####

populations_order <- "INDO,NHV-S,NHV-P,UAH,SCI,MOP,AHE,MAN,TKP,ARA,KAU,RAR,ANA,KAT,TAH,MOT,TEA,TAK,RVV,MARS,MOR,RIK,GMBW,PEI,MAT"
labels <- as.data.frame(cbind(ind= indNames(obj.glx), pop = as.character(pop(obj.glx))))

labels$n<-factor(labels$pop,levels=unlist(strsplit(populations_order,",")))
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))
rep<-as.vector(table(labels$n))
rep_pop <-0
for (i in 1:length(levels(as.factor(labels$pop)))) {
  rep_pop <- c(rep_pop,
               rep((as.vector(strsplit(populations_order,",")))[[1]][i],
                   rep[i]))
}

rep_pop_K1 <- rep(rep_pop[-1],1)
rep_pop_K2 <- rep(rep_pop[-1],2)
rep_pop_K3 <- rep(rep_pop[-1],3)
rep_pop_K4 <- rep(rep_pop[-1],4)
rep_pop_K5 <- rep(rep_pop[-1],5)

# read in the different admixture output files
directory_prefix <- "~/all/neutral_dataset_all"
minK=1
maxK=5
tbl<-lapply(minK:maxK, function(x) read.table(paste0(directory_prefix,".",x,".Q")))

#K1
K1_df <- as.data.frame(cbind(tbl[[1]][order(labels$n),],
                             Categorie=1:dim(tbl[[1]])[1]))
K1.df.long <- tidyr:: gather(K1_df, "Cluster", "Proportion", -Categorie)
K1.df.long.pop <- as.data.frame(cbind(K1.df.long,
                                      Population = factor(rep_pop_K1,
                                                          levels = (as.vector(strsplit(populations_order,",")))[[1]])))


K1 <- ggplot(K1.df.long.pop, aes(x=Categorie, y=Proportion, fill=Cluster))  + ylab("Assignment") +xlab("Population")
K1 <- K1 + geom_bar(stat='identity',width=1)
K1 <- K1 + scale_fill_manual(values = c("red","green","blue"),name="cluster") 
K1 <- K1 + facet_grid(~K1.df.long.pop$Population, scales = "free")
K1 <- K1 + theme_gray(base_size = 12) + guides(fill=guide_legend(ncol=1))
K1 <- K1 + theme(axis.text.x = element_blank(),
                 axis.ticks = element_blank(),
                 panel.spacing.x=unit(0.1, "lines"),
                 strip.text.x = element_text(angle = 90))

#K2
K2_df <- as.data.frame(cbind(tbl[[2]][order(labels$n),],
                             Categorie=1:dim(tbl[[2]])[1]))
K2.df.long <- tidyr:: gather(K2_df, "Cluster", "Proportion", -Categorie)
K2.df.long.pop <- as.data.frame(cbind(K2.df.long,
                                      Population = factor(rep_pop_K2,
                                                          levels = (as.vector(strsplit(populations_order,",")))[[1]])))


K2 <- ggplot(K2.df.long.pop, aes(x=Categorie, y=Proportion, fill=Cluster))  + ylab("Assignment") +xlab("Population")
K2 <- K2 + geom_bar(stat='identity',width=1)
K2 <- K2 + scale_fill_manual(values = c("#F8766D","#E0AC5E"),name="cluster") 
K2 <- K2 + facet_grid(~K2.df.long.pop$Population, scales = "free")
K2 <- K2 + theme_gray(base_size = 12) + guides(fill=guide_legend(ncol=1))
K2 <- K2 + theme(axis.text.x = element_blank(),
                 axis.ticks = element_blank(),
                 panel.spacing.x=unit(0.1, "lines"),
                 strip.text.x = element_text(angle = 90))

#K3
K3_df <- as.data.frame(cbind(tbl[[3]][order(labels$n),],
                             Categorie=1:dim(tbl[[3]])[1]))
K3.df.long <- tidyr:: gather(K3_df, "Cluster", "Proportion", -Categorie)
K3.df.long.pop <- as.data.frame(cbind(K3.df.long,
                                      Population = factor(rep_pop_K3,
                                                          levels = (as.vector(strsplit(populations_order,",")))[[1]])))


K3 <- ggplot(K3.df.long.pop, aes(x=Categorie, y=Proportion, fill=Cluster))  + ylab("Assignment") +xlab("Population")
K3 <- K3 + geom_bar(stat='identity',width=1)
K3 <- K3 + scale_fill_manual(values = c("#E0AC5E","#00BA38","#F8766D"),name="cluster") 
K3 <- K3 + facet_grid(~K3.df.long.pop$Population, scales = "free")
K3 <- K3 + theme_gray(base_size = 12) + guides(fill=guide_legend(ncol=1))
K3 <- K3 + theme(axis.text.x = element_blank(),
                 axis.ticks = element_blank(),
                 panel.spacing.x=unit(0.1, "lines"),
                 strip.text.x = element_text(angle = 90))

#K4
K4_df <- as.data.frame(cbind(tbl[[4]][order(labels$n),],
                             Categorie=1:dim(tbl[[4]])[1]))
K4.df.long <- tidyr:: gather(K4_df, "Cluster", "Proportion", -Categorie)
K4.df.long.pop <- as.data.frame(cbind(K4.df.long,
                                      Population = factor(rep_pop_K4,
                                                          levels = (as.vector(strsplit(populations_order,",")))[[1]])))


K4 <- ggplot(K4.df.long.pop, aes(x=Categorie, y=Proportion, fill=Cluster))  + ylab("Assignment") +xlab("Population")
K4 <- K4 + geom_bar(stat='identity',width=1)
K4 <- K4 + scale_fill_manual(values = c("#F8766D","#E0AC5E","#00BA38", "#C27E42"),name="cluster") 
K4 <- K4 + facet_grid(~K4.df.long.pop$Population, scales = "free")
K4 <- K4 + theme_gray(base_size = 12) + guides(fill=guide_legend(ncol=1))
K4 <- K4 + theme(axis.text.x = element_blank(),
                 axis.ticks = element_blank(),
                 panel.spacing.x=unit(0.1, "lines"),
                 strip.text.x = element_text(angle = 90))

#K5
K5_df <- as.data.frame(cbind(tbl[[5]][order(labels$n),],
                             Categorie=1:dim(tbl[[5]])[1]))
K5.df.long <- tidyr:: gather(K5_df, "Cluster", "Proportion", -Categorie)
K5.df.long.pop <- as.data.frame(cbind(K5.df.long,
                                      Population = factor(rep_pop_K5,
                                                          levels = (as.vector(strsplit(populations_order,",")))[[1]])))


K5 <- ggplot(K5.df.long.pop, aes(x=Categorie, y=Proportion, fill=Cluster))  + ylab("Assignment") +xlab("Population")
K5 <- K5 + geom_bar(stat='identity',width=1)
K5 <- K5 + scale_fill_manual(values = c("#00BA38","#C27E42","#9A72AB", "#E0AC5E", "#F8766D"),name="cluster") 
K5 <- K5 + facet_grid(~K5.df.long.pop$Population, scales = "free")
K5 <- K5 + theme_gray(base_size = 12) + guides(fill=guide_legend(ncol=1))
K5 <- K5 + theme(axis.text.x = element_blank(),
                 axis.ticks = element_blank(),
                 panel.spacing.x=unit(0.1, "lines"),
                 strip.text.x = element_text(angle = 90))


admixture_1 <- ggarrange(K2,K3, ncol = 1, nrow = 2, labels = c("K=2",
                                                               "K=3"))
admixture_2 <- ggarrange(K4,K5, ncol = 1, nrow = 2, labels = c("K=4",
                                                               "K=5"))

admixture_3 <- ggarrange(K2, K3, K4,K5, ncol = 1, nrow = 4, labels = c("K=2","K=3", "K=4","K=5"))
admixture_1
admixture_2
admixture_3
