###########################################
#  - TFMir  Project                       #
#  - prepare DBs of TF-miRNA interactions #
#  - 2014-10-1                            #
#  - Copyright: Mohamed Hamed             #
###########################################

#R programming environments:
#- R studio IDE
#- R version 2.12.0 (2010-10-15)
#-Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
#-locale:   [1] C/en_US.UTF-8/C/C/C/C


## ===============================
## read the transmir database file 
## ===============================
## 1- convert the transmir xls file  to tab delimited file cause it causes problems with wild characters
transmir <- read.delim(config$transmir.file)
transmir["organism2"]=""
transmir$organism2=apply(as.data.frame(transmir$organism),1, function(x) trim(as.character(x)))
transmir=transmir[transmir$organism2 %in% c("human","Human"),]
transmir["mirna"]=""
transmir$mirna=apply(as.data.frame(transmir$mir),1, function(x)  paste("hsa-",x,sep="")  )
transmir["evidence"]="Experimental"
transmir["source"]="Transmir"
transmir=transmir[,c("gene","mirna","evidence","source")]


## ======================================
## read the PMID20584335 interaction file 
## ======================================
pmid <- read.delim(config$PMID20584335.file)
pmid=pmid[pmid$category=="TF-miRNA",]
pmid=pmid[! duplicated(pmid),]
pmid=pmid[,c(2:5)]
names(pmid)=c("gene","mirna","evidence","source")

## ===============================
## read the Chipbase database file 
## ===============================

chipbase.df <- read.delim(config$chipbase.file)
chipbase=data.frame("mirna"=character(),"gene"=character()) 
for(i in 1: dim(chipbase.df)[1])
{
  s=as.character(chipbase.df[i,]$miRNAs)
  strs=strsplit(s,",")
  df=data.frame("mirna"= strs[[1]][1:length(strs[[1]])])
  df["gene"]=as.character(chipbase.df[i,]$tfName)
  chipbase=rbind(chipbase,df)
}
chipbase=chipbase[,c(2,1)]
chipbase["evidence"]="Predicted"
chipbase["source"]="Chipbase"
names(chipbase)
chipbase=chipbase[! duplicated(chipbase),]




tf.mirna.db=rbind(transmir,pmid,chipbase)
d=tf.mirna.db[,c(1,2)]
tf.mirna.db=tf.mirna.db[! duplicated(d),]
tf.mirna.db["category"]="tf-mirna"
tf.mirna.db=tf.mirna.db[,c(5,1,2,3,4)]














# transmir["tumor_gene"]=""
# transmir[which(transmir$tumor!="Y"),]$tumor_gene  ="UNKNOWN"
# transmir[is.na(transmir$tumor),]$tumor_gene  ="UNKNOWN"
# transmir[which(transmir$tumor=="Y"),]$tumor_gene  ="TRUE"
# transmir["tumor_mirna"]=""
# transmir[which(transmir$tumor_mir!="Y"),]$tumor_mirna  ="UNKNOWN"
# transmir[is.na(transmir$tumor_mir),]$tumor_mirna  ="UNKNOWN"
# transmir[which(transmir$tumor_mir=="Y"),]$tumor_mirna  ="TRUE"
# transmir=transmir[,c(1,2,13,12,14,6,7,8,9,11)]
# names(transmir)