vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))


### Setting input paths for RC, mappability, GC-content and target ###
ProgramFolder <- split.vars[1]
PathInVec <- split.vars[2]
ExpName <- split.vars[3]
TargetName <- split.vars[4]
Assembly <- split.vars[5]

TargetFolder<-file.path(ProgramFolder,"data","targets",Assembly,TargetName)

### Setting paths for RC, mappability, GC-content e target ###
PathGC <- file.path(TargetFolder,"GCC")
PathMAP <- file.path(TargetFolder,"MAP")
GCFiles <- list.files(PathGC)
MAPFiles <- list.files(PathMAP)


### Loading target chromosomes ###
TargetChrom <- file.path(TargetFolder,paste(TargetName,"_chromosome.txt",sep=""))
CHR<-readLines(con = TargetChrom, n = -1L, ok = TRUE, warn = TRUE,encoding = "unknown")
unique.chrom<-strsplit(CHR," ")[[1]]



Path2ExomeRC <- file.path(ProgramFolder,"lib","R","LibraryGenomeRC.R")
source(Path2ExomeRC)


### Loading target, mappability and GC content Files ###
TargetOut <- loadTarget(TargetFolder,unique.chrom,TargetName)
MyTarget <- TargetOut[[1]]
GCContentT <- TargetOut[[2]]*100
MAPT <- TargetOut[[3]]*100


start <- as.integer(MyTarget[,2])
end <- as.integer(MyTarget[,3])
chrom <- as.character(MyTarget[,1])
GeneName <- as.character(MyTarget[,4])
Position <- (start+end)/2
L<-end-start


GCContentL <- GCContentT


### Normalization for GC Content and mappability ###
library(Hmisc)
PathRC <- file.path(PathInVec,"RC")
RCFiles <- list.files(PathRC)
RCTL <- loadRC(PathRC,RCFiles,unique.chrom)


PathNorm <- file.path(PathInVec,"RCNorm")
PathImages <- file.path(PathInVec,"Images")
FileOut <- file.path(PathNorm,paste(ExpName,".NRC.RData",sep=""))



### Mappability normalization ###
step <- 5
RCMAPNormList <- CorrectMAP(RCTL,MAPT,step)
RCMAPNorm <- RCMAPNormList$RCNorm


### GC-content Normalization ###
step <- 5
RCGCNormList <- CorrectGC(RCMAPNorm,GCContentL,step)
RCGCNorm <- RCGCNormList$RCNorm


### Images generation ###

### Quantiles for mappability ###
step <- 5
QuantileMapList <- QuantileMAP(RCTL,MAPT,step)
MedianVecMAP <- QuantileMapList$Median
QuantileVecMAP <- QuantileMapList$Quantile

MAPBias <- file.path(PathImages,"MAPBias.pdf")
pdf(MAPBias,height=10,width=15)
yup <- median(RCTL,na.rm=T)*3
stepseq <- seq(0,100,by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "Mappability"
labely <- "Read Count"
errbar(xstep,MedianVecMAP,QuantileVecMAP[,2],QuantileVecMAP[,1],axes=F,frame.plot=TRUE,xlab=labelx,ylab=labely,add=F,ylim=c(0,yup),lty=3,cex.axis=1.2,cex.lab=1.3)
title("Mappability Bias",cex.main=2)
text((length(stepseq)-3),2000,"a",cex=3)
ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"))
axis(2)
dev.off()

#### Quantiles for GC-content ###
step <- 5
QuantileGCList <- QuantileGC(RCTL,GCContentL,step)
MedianVecGC <- QuantileGCList$Median
QuantileVecGC <- QuantileGCList$Quantile

GCBias <- file.path(PathImages,"GCBias.pdf")
pdf(GCBias,height=10,width=15)
yup <- median(RCTL,na.rm=T)*2
#par(mfrow=c(3,5),oma=c(2,0,2,0))
stepseq <- seq(0,100,by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "GC-Content"
labely <- "Read Count"
errbar(xstep,MedianVecGC,QuantileVecGC[,2],QuantileVecGC[,1],axes=F,frame.plot=TRUE,xlab=labelx,ylab=labely,add=F,ylim=c(0,yup),lty=3,cex.axis=1.2,cex.lab=1.3)
#title("Mappability Normalization",cex.main=2.7,outer=T)
title("GC-Content Bias",cex.main=2)
text((length(stepseq)-3),2000,"a",cex=3)
ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"))
axis(2)
dev.off()

### Normalized data ###


### Quantile of normalized mapping ###
step <- 5
QuantileMapListN <- QuantileMAP(RCMAPNorm,MAPT,step)
MedianVecMAPN <- QuantileMapListN$Median
QuantileVecMAPN <- QuantileMapListN$Quantile
MAPBiasNorm <- file.path(PathImages,"MAPBiasNorm.pdf")
pdf(MAPBiasNorm,height=10,width=15)
yup <- median(RCMAPNorm,na.rm=T)*2
#par(mfrow=c(3,5),oma=c(2,0,2,0))
stepseq <- seq(0,100,by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "Mappability"
labely <- "Read Count"
errbar(xstep,MedianVecMAPN,QuantileVecMAPN[,2],QuantileVecMAPN[,1],axes=F,frame.plot = TRUE,xlab=labelx,ylab=labely,add=F,ylim=c(0,yup),lty=3,cex.axis=1.2,cex.lab=1.3)
#title("Mappability Normalization",cex.main=2.7,outer=T)
title("Normalized Read Count Mappability Bias",cex.main=2)
text((length(stepseq)-3),2000,"a",cex=3)
ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"))
axis(2)
dev.off()


#### Quantile of normalized GC-content ###
step <- 5
QuantileGCListN <- QuantileGC(RCGCNorm,GCContentL,step)
MedianVecGCN <- QuantileGCListN$Median
QuantileVecGCN <- QuantileGCListN$Quantile

GCBiasNorm <- file.path(PathImages,"GCBiasNorm.pdf")
pdf(GCBiasNorm,height=10,width=15)

yup <- median(RCGCNorm,na.rm=T)*2
#par(mfrow=c(3,5),oma=c(2,0,2,0))
stepseq <- seq(0,100,by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "GC-Content"
labely <- "Read Count"
errbar(xstep,MedianVecGCN,QuantileVecGCN[,2],QuantileVecGCN[,1],axes=F,frame.plot = TRUE,xlab=labelx,ylab=labely,add=F,ylim=c(0,yup),lty=3,cex.axis=1.2,cex.lab=1.3)
#title("Mappability Normalization",cex.main=2.7,outer=T)
title("Normalizard Read Count GC-Content Bias",cex.main=2)
text((length(stepseq)-3),2000,"a",cex=3)
ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"))
axis(2)
dev.off()

##### Saving normalized RC file ###
RCGCNorm[which(RCGCNorm==0)] <- min(RCGCNorm[which(RCGCNorm!=0)])
RCNorm <- RCGCNorm
MatrixNorm <- cbind(chrom,Position,start,end,GeneName,RCNorm)
save(MatrixNorm,file=FileOut)





