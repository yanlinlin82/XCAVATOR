vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))
ProgramFolder <- split.vars[1]
target.name <- split.vars[2]
assembly <- split.vars[3]
step <- as.numeric(split.vars[4])
FastaIn <- split.vars[5]

options("scipen"=20)

if (assembly=="hg19"){
  CoordIn <- paste(ProgramFolder,"data/support",assembly,"ChromosomeCoordinate_HG19.txt",sep="/")
  FileGap <- paste(ProgramFolder,"data/support",assembly,"GapHg19.UCSC.txt",sep="/") 
}

if (assembly=="hg38"){
  CoordIn <- paste(ProgramFolder,"data/support",assembly,"ChromosomeCoordinate_HG38.txt",sep="/")
  FileGap <- paste(ProgramFolder,"data/support",assembly,"GapHg38.UCSC.txt",sep="/") 
}

FastaInParseString<-paste("head -1 ",FastaIn," | awk '{print $1}'| tr -d '>'",sep="")
ChrLab<-system(FastaInParseString,intern =TRUE)

CoordTable<-read.table(CoordIn,sep="\t",quote="\"",fill=T,header=F)
ChrCoord<-as.character(CoordTable[,1])
StartCoord<-as.numeric(CoordTable[,2])
EndCoord<-as.numeric(CoordTable[,3])



if (nchar(ChrLab)<3)
{
  ChrVec<-c(1:22,"X","Y")
  ChrCoord<-c(1:22,"X","Y")
}
if (nchar(ChrLab)>3)
{
  ChrVec<-paste("chr",c(1:22,"X","Y"),sep="")
  ChrCoord<-paste("chr",c(1:22,"X","Y"),sep="")
}



#### Creating Bed File ####
GapCoordTable<-read.table(FileGap,sep="\t",quote="\"",fill=T,header=F)

ChrTable<-as.character(GapCoordTable[,2])
StartTable<-as.numeric(GapCoordTable[,3])
EndTable<-as.numeric(GapCoordTable[,4])

ChrVecSupp<-paste("chr",c(1:22,"X","Y"),sep="")


TotalMatrixCoord<-c()
for (i in 1:length(ChrVecSupp))
{
  indC<-which(ChrTable==ChrVecSupp[i])
  StartTableC<-StartTable[indC]
  EndTableC<-EndTable[indC]
  indS<-sort(StartTableC,index.return=T)$ix
  StartTableCS<-StartTableC[indS]
  EndTableCS<-EndTableC[indS]
  
  CoordVecC<-cbind(EndTableCS[-length(EndTableCS)]+1,StartTableCS[-1])
  WindowCoord<-c()
  for (k in 1:nrow(CoordVecC))
  {
    if ((CoordVecC[k,2]-CoordVecC[k,1])>step)
    {
      CoordNew<-seq(CoordVecC[k,1],CoordVecC[k,2],by=step)
      if (length(CoordNew)>1)
      {
        WindowCoord<-rbind(WindowCoord,cbind(CoordNew[-length(CoordNew)],CoordNew[-1]-1))
        
      }
    }
  }
  TotalMatrixCoord<-rbind(TotalMatrixCoord,cbind(rep(ChrVec[i],nrow(WindowCoord)),WindowCoord))
}


#new format
a<-paste("a",seq(1,nrow(TotalMatrixCoord)),sep="")
FinalTargetF<-cbind(TotalMatrixCoord,a)








MyTarget <- FinalTargetF
MyChr <- unique(FinalTargetF[,1])


dir.create(file.path(ProgramFolder,"data","targets",assembly))
dir.create(file.path(ProgramFolder,"data","targets",assembly,target.name))
dir.create(file.path(ProgramFolder,"data","targets",assembly,target.name,"GCC"))
dir.create(file.path(ProgramFolder,"data","targets",assembly,target.name,"MAP"))
dir.create(file.path(ProgramFolder,"data","targets",assembly,target.name,"FRB"))
write.table(data.frame(rbind(MyChr)),file.path(ProgramFolder,"data","targets",assembly,target.name,paste(target.name,"_chromosome.txt",sep="")),col.names=F,row.names=F,quote=F)
save(MyTarget,file=file.path(ProgramFolder,"data","targets",assembly,target.name,paste(target.name,".RData",sep="")))
write.table(MyTarget,file.path(ProgramFolder,"data","targets",assembly,target.name,paste("Filtered.txt",sep="")),col.names=F,row.names=F,sep="\t",quote=F)






