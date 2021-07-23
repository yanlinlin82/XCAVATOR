###################### Starting Parameter ####################
ParamEstSeq <- function(DataMatrix,omega)
{
	T=ncol(DataMatrix)
	NExp<-nrow(DataMatrix)
	
	sigmax<-c()
	mi<-c()
	
	for (i in 1:NExp)
	{
	  mi[i]<-0
	  sigmax[i]<-mad(DataMatrix[i,])^2
	}
	
	smu<-sqrt(omega*sigmax)
	sepsilon<-sqrt((1-omega)*sigmax)
	Results<-list()
	Results$mi<-mi
	Results$smu<-smu
	Results$sepsilon<-sepsilon
	
	Results
}


###################### MUK Estimation ####################
MukEst <- function(DataMatrix,mw)
{
	NExp<-dim(DataMatrix)[1]
	if (NExp==1)
	{
		muk<-rbind(seq(-1,1,by=0.1))
	}
	if (NExp>1)
	{
		DataMatrix[which(DataMatrix>1)]<-1
		DataMatrix[which(DataMatrix< -1)]<- -1
		
		DataMatrixA<-c()
		for (i in 1:NExp)
		{
			DataMatrixA<-rbind(DataMatrixA,SMA(DataMatrix[i,], n=mw))
		}
		
		DataMatrixA<-DataMatrixA[,2:length(DataMatrixA[1,])]
		
		
		binsize=0.2
		binVec<-c(seq(-1,-0.2,by=binsize),0,seq(0.2,1,by=binsize))
		binmean<-c(seq(-0.9,-0.3,by=binsize),0,0,seq(0.3,0.9,by=binsize))
		
		
		DataQuant<-DataMatrixA
		
		for (i in 1:(length(binVec)-1))
		{
			DataQuant[which(DataMatrixA>binVec[i] & DataMatrixA<=binVec[i+1])]<-binmean[i]
		}
		
		muk<-unique(DataQuant,  MARGIN = 2)
		muk<-muk[,-1]
	}
	muk
}

###################### Segmentation Core ####################
JointSeg <- function(DataMatrix,eta,omega,muk,mi,smu,sepsilon)
{
	
	T=ncol(DataMatrix)
	K0<-ncol(muk)
	etav<-log(rep(1,K0)*(1/K0))
	NExp<-nrow(DataMatrix)
	
	
	
####  Transition and Emission Calculation ####
	P<-matrix(data=0,nrow=K0,ncol=K0)
	G<-matrix(data=0,nrow=K0,ncol=K0)
	emission<-matrix(data=0,nrow=K0,ncol=T)
	out<-.Fortran("transemis",as.vector(muk),as.vector(mi),as.double(eta),as.matrix(DataMatrix),as.integer(K0),as.integer(NExp),as.vector(smu),as.vector(sepsilon),as.integer(T),as.matrix(G),as.matrix(P),as.matrix(emission))
	
	P<-out[[11]]
	emission<-out[[12]]
	
	
	
####### Viterbi algorithm #########
	psi<-matrix(data=0,nrow=K0,ncol=T)
	path<-c(as.integer(rep(0,T)))
	out2<-.Fortran("bioviterbi",as.vector(etav),as.matrix(P),as.matrix(emission),as.integer(T),as.integer(K0),as.vector(path),as.matrix(psi))
	s<-out2[[6]]
	
	
	
	sortResult <- SortState(s)
	TotalPredBrek<-sortResult[[3]]
	TotalPredBrek
}


################# State Sorting #########
SortState <- function(s)
{
	l<-1
	seg<-c()
	brek<-c()
	t<-1
	for (k in 1:(length(s)-1))
	{
		if (s[k]!=s[k+1])
		{
			brek[t]<-k
			t<-t+1
			if (length(which(seg==s[k]))==0)
			{
				seg[l]<-s[k]
				l<-l+1
			}
		}
	}
	brek<-c(0,brek,length(s))
	if (length(which(seg==s[length(s)]))==0)
	{
		seg<-c(seg,s[length(s)])
	}
	
	s0<-c()
	for (k in 1:length(seg))
	{
		s0[which(s==seg[k])]<-k
	}
	
	SortResult<-list()
	SortResult[[1]]<-s0
	SortResult[[2]]<-seg
	SortResult[[3]]<-brek
	SortResult
	
}


###### Filtering small shifts #####
FilterSeg <- function(TotalPredBreak,FW)
{
  controllength<-diff(TotalPredBreak)
  
  indF<-which(controllength<=FW)
  if (length(indF)!=0)
  {
    if (indF[1]==1)
    {
      indF[1]<-2
      indF<-unique(indF)
      TotalPredBreak1 <- TotalPredBreak[-(indF)]
    }
    if (indF[1]!=1)
    {
      TotalPredBreak1 <- TotalPredBreak[-(indF)]
    }
  }
if (length(indF)==0)
{
  TotalPredBreak1<-TotalPredBreak
}
TotalPredBreak1
}



################# Function for calculating segment mean after segmentation #########
SegResults <- function(DataSeq,TotalPredBreak)
{
	
	TotalPred<-c()
	NExp<-nrow(DataSeq)
	for (j in 1:NExp)
	{
		s<-rep(0,ncol(DataSeq))
		for (i in 1:(length(TotalPredBreak)-1))
		{
			s[(TotalPredBreak[i]+1):TotalPredBreak[i+1]]<-median(DataSeq[j,(TotalPredBreak[i]+1):TotalPredBreak[i+1]])
		}
		TotalPred<-rbind(TotalPred,s)
	}
	
	Result<-TotalPred
	Result
}


################# Function for calculating statistics on log2ratio data #########
DataStatistics <- function(LogDataNorm,DataFolder,ExpLabelOut)
{
  minD<-floor(min(LogDataNorm))
  maxD<-ceiling(max(LogDataNorm))
  breaks<-seq(minD,maxD,by=0.05)
  Miohist <- hist(LogDataNorm, breaks=breaks,plot=FALSE)
  MioQuantile<-quantile(LogDataNorm,probs = c(0.01,0.1,0.25,0.5,0.75,0.9,0.99))
  Counts<-Miohist$counts
  colors = c("black","red","orange","yellow","green","blue","purple")
  DataNoise<-LogDataNorm[which(LogDataNorm>MioQuantile[1] & LogDataNorm<MioQuantile[7])]
  SNR3<-0.58/sd(DataNoise)
  SNR1<-1/sd(DataNoise)
  SNRVec<-c(SNR1,SNR3)
  names(SNRVec)<-c("SNR 1-copy","SNR 3-copy")
  ColorData<-rep(colors[1],length(breaks))
  ColorData[which(breaks<MioQuantile[1])]<-colors[2]
  ColorData[which(breaks>MioQuantile[1] & breaks<MioQuantile[2])]<-colors[3]
  ColorData[which(breaks>MioQuantile[2] & breaks<MioQuantile[3])]<-colors[4]
  
  ColorData[which(breaks>MioQuantile[5] & breaks<MioQuantile[6])]<-colors[5]
  ColorData[which(breaks>MioQuantile[6] & breaks<MioQuantile[7])]<-colors[6]
  ColorData[which(breaks>MioQuantile[7])]<-colors[7]
  xlimL<-tail(which(breaks< -2.5),1)
  xlimU<-head(which(breaks> 2.5),1)
  
  
  
  
  FilePercentileOut<-file.path(DataFolder,"Results",ExpLabelOut,paste("Percentile_",ExpLabelOut,".txt",sep=""))
  FileNoiseOut<-file.path(DataFolder,"Results",ExpLabelOut,paste("Noise_",ExpLabelOut,".txt",sep=""))
  FileDistributionOut<-file.path(DataFolder,"Results",ExpLabelOut,paste("PercentileDistribution_",ExpLabelOut,".pdf",sep=""))
  
  pdf(FileDistributionOut)
  LegendLabel<-c("50th Percentile","1th Percentile","10th Percentile","25th Percentile","75th Percentile","90th Percentile","99th Percentile")
  barplot(Counts[c(xlimL:xlimU)],col=ColorData[c(xlimL:xlimU)],names.arg=breaks[c(xlimL:xlimU)],xlab="log2Ratio")
  legend(2.3,0.9*max(Counts),legend=LegendLabel,fill=colors)
  dev.off()
  
  write.table(rbind(MioQuantile),FilePercentileOut,col.names=T,row.names=F,sep="\t",quote=F)
  write.table(rbind(SNRVec),FileNoiseOut,col.names=T,row.names=F,sep="\t",quote=F)
  
}