vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
Sample_Name <- split.vars [1]
Output_Folder <- split.vars[2]
Program_Folder <- split.vars[3]
chr <- split.vars[4]
target.name <- split.vars[5]
FileBAMIn <- split.vars[6]

DepthString<-"/home/Workspace/Tools/samtools/bin/samtools bedcov "

FileBEDIn<-file.path(Output_Folder,".tmp",".temp.bed")

SamtoolsDepthString<-paste(DepthString,FileBEDIn," ",FileBAMIn,sep="")

TableDepthOut<-system(SamtoolsDepthString,intern = TRUE)

TableDepthOutSplit<-unlist(strsplit(TableDepthOut,"\t"))
RC<-as.numeric(TableDepthOutSplit[seq(5,length(TableDepthOutSplit),by=5)])


save(RC,file = file.path(Output_Folder,"RC",paste(Sample_Name,".RC.",chr,".RData",sep="")))
