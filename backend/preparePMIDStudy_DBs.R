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


## ===================================
## read the PMID TF-mirna interactions
## ===================================

pmid <- read.delim("databases/tf-mirna/PMID20584335.txt")
pmid$miRNA=paste("hsa-",tolower(pmid$miRNA),sep="")
pmid["gene"]=as.vector(unlist(mget(as.character(pmid$EntrezID), envir=org.Hs.egSYMBOL, ifnotfound=NA)))
pmid=pmid[,c(4,2)]
names(pmid)=c("FROM","TO")
pmid=pmid[! duplicated(pmid),]
pmid.tf.mirna=pmid
pmid.tf.mirna["category"]="TF-miRNA"
pmid.tf.mirna=pmid.tf.mirna[,c(3,1,2)]

## ===================================
## read the PMID mirna-gene interactions
## ===================================

pmid <- read.delim("databases/mirna-gene/PMID20584335.txt")
pmid$mirna=paste("hsa-",tolower(pmid$mirna),sep="")
names(pmid)=c("FROM","TO")
pmid=pmid[! duplicated(pmid),]
pmid.mirna.gene=pmid
pmid.mirna.gene["category"]="miRNA-gene"
pmid.mirna.gene=pmid.mirna.gene[,c(3,1,2)]



## ===================================
## read the PMID TF-gene interactions
## ===================================

pmid <- read.delim("databases/tf-gene/PMID20584335.txt")
names(pmid)=c("FROMTF","TOGENE")
pmid=pmid[! duplicated(pmid),]
pmid["FROM"]=as.vector(unlist(mget(as.character(pmid$FROMTF), envir=org.Hs.egSYMBOL, ifnotfound=NA)))
pmid["TO"]=as.vector(unlist(mget(as.character(pmid$TOGENE), envir=org.Hs.egSYMBOL, ifnotfound=NA)))
pmid.tf.gene=pmid
pmid.tf.gene["category"]="TF-gene"
pmid.tf.gene=pmid.tf.gene[,c("category","FROM","TO")]


pmid=rbind(pmid.tf.mirna,pmid.tf.gene,pmid.mirna.gene)
pmid["evidence"]="Experimental"
pmid["source"]="PMID20584335"


write.table(pmid,"databases/PMID20584335.txt",quote=F,col.names=T,row.names=F,sep="\t")
