###########################################
#  - TFMir  Project                       #
#  - prepare DBs of miRNA-genes interacti #
#  - 2014-10-1                            #
#  - Copyright: Mohamed Hamed             #
###########################################

#R programming environments:
#- R studio IDE
#- R version 2.12.0 (2010-10-15)
#-Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
#-locale:   [1] C/en_US.UTF-8/C/C/C/C



## ===============================
## read the Experimental databases  
## ===============================
mirtarbase = read.xls (config$miRNA.Targets.DBs.file, sheet = 4, header = TRUE)
mirtarbase=mirtarbase[mirtarbase$Species..miRNA.=="Homo sapiens",]
mirtarbase=mirtarbase[- agrep("weak",mirtarbase$Support.Type),]
mirtarbase=mirtarbase[- agrep("weaK",mirtarbase$Support.Type),]
mirtarbase=mirtarbase[,c("miRNA","Target.Gene")]
names(mirtarbase)=c("mirna","gene")
mirtarbase$mirna=tolower(mirtarbase$mirna)
mirtarbase=mirtarbase[!duplicated(mirtarbase),]
mirtarbase["evidence"]="Experimental"
mirtarbase["source"]="mirtarbase"



tarbase = read.xls (config$miRNA.Targets.DBs.file, sheet = 5, header = TRUE)
tarbase=tarbase[tarbase$Organism=="Human",]
tarbase=tarbase[tarbase$Support_Type != FALSE,]
tarbase=tarbase[,c("miRNA","Gene")]
names(tarbase)=c("mirna","gene")
tarbase$mirna=paste("hsa-",tolower(tarbase$mirna),sep="")
tarbase=tarbase[!duplicated(tarbase),]
tarbase.df=data.frame("mirna"=character(),"gene"=character()) 
for(i in 1: dim(tarbase)[1])
{
  s=as.character(tarbase[i,]$gene)
  strs=strsplit(s," / ")
  df=data.frame(gene= strs[[1]][1:length(strs[[1]])])
  df["mirna"]=tarbase[i,]$mirna
  tarbase.df=rbind(tarbase.df,df)
}
tarbase.df=tarbase.df[,c(2,1)]
tarbase=tarbase.df
tarbase["evidence"]="Experimental"
tarbase["source"]="tarbase"




mirrecord = read.xls (config$miRNA.Targets.DBs.file, sheet = 7, header = TRUE)
mirrecord=mirrecord[mirrecord$miRNA_species=="Homo sapiens",]
mirrecord=mirrecord[mirrecord$Target.gene_species_scientific=="Homo sapiens",]
mirrecord=mirrecord[,c("miRNA_mature_ID","Target.gene_name")]
names(mirrecord)=c("mirna","gene")
mirrecord$mirna=tolower(mirrecord$mirna)
mirrecord=mirrecord[!duplicated(mirrecord),]
mirrecord["evidence"]="Experimental"
mirrecord["source"]="mirRecord"



mirtarbase.adstudy = read.xls (config$miRNA.Targets.DBs.file, sheet = 2, header = TRUE)
mirtarbase.adstudy=mirtarbase.adstudy[,c(1,2)]
names(mirtarbase.adstudy)=c("mirna","gene")
mirtarbase.adstudy$mirna=tolower(mirtarbase.adstudy$mirna)
mirtarbase.adstudy["evidence"]="Experimental"
mirtarbase.adstudy["source"]="mirtarbase.adstudy"



tarbase.adstudy = read.xls (config$miRNA.Targets.DBs.file, sheet = 3, header = TRUE)
tarbase.adstudy=tarbase.adstudy[,c(1,2)]
names(tarbase.adstudy)=c("mirna","gene")
tarbase.adstudy$mirna=tolower(tarbase.adstudy$mirna)
tarbase.adstudy["evidence"]="Experimental"
tarbase.adstudy["source"]="tarbase.adstudy"


mirrecords.adstudy = read.xls (config$miRNA.Targets.DBs.file, sheet = 1, header = TRUE)
mirrecords.adstudy=mirrecords.adstudy[,c(1,2)]
names(mirrecords.adstudy)=c("mirna","gene")
mirrecords.adstudy$mirna=tolower(mirrecords.adstudy$mirna)
mirrecords.adstudy["evidence"]="Experimental"
mirrecords.adstudy["source"]="mirrecords.adstudy"




## ======================================
## read the PMID20584335 interaction file 
## ======================================
pmid <- read.delim(config$PMID20584335.file)
pmid=pmid[pmid$category=="miRNA-gene",]
pmid=pmid[! duplicated(pmid),]
pmid=pmid[,c(2:5)]
names(pmid)=c("mirna","gene","evidence","source")


## ===========================================
## read the predicted interactions by Starbase
## ===========================================

starbase <- read.delim(config$miRNA.Starbase.file)
starbase=starbase[,c(1,2)]
names(starbase)=c("mirna","gene")
starbase$mirna=tolower(starbase$mirna)
starbase["evidence"]="Predicted"
starbase["source"]="Starbase"



mirna.gene.db=rbind(mirtarbase.adstudy,tarbase.adstudy,mirrecords.adstudy,mirtarbase,tarbase,mirrecord,pmid,starbase)
d=mirna.gene.db[,c(1,2)]
mirna.gene.db=mirna.gene.db[! duplicated(d),]
mirna.gene.db["category"]="mirna-gene"
mirna.gene.db=mirna.gene.db[,c(5,1,2,3,4)]



