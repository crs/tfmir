###########################################
#  - TFMir  Project                       #
#  - prepare DBs of TF-genes interactions #
#  - 2014-10-1                            #
#  - Copyright: Mohamed Hamed             #
###########################################

#R programming environments:
#- R studio IDE
#- R version 2.12.0 (2010-10-15)
#-Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
#-locale:   [1] C/en_US.UTF-8/C/C/C/C




## ==================================
## read the Transfac interaction file 
## ==================================
transfac = read.xls (config$transfac.file, sheet = 1, header = TRUE)
transfac=transfac[,c(1,4)]
transfac["evidence"]="Experimental"
transfac["source"]="Transfac"
transfac=transfac[!duplicated(transfac),]
names(transfac)=c("tf","gene","evidence","source")


## ==================================
## read the ITFP interaction file 
## ==================================
files= list.files(path=config$ITFP.folder,pattern = "txt")
itfb=data.frame("tf"=character(),"gene"=character()) 
for (i in 1:length(files))
{
  path=paste(config$ITFP.folder ,"/", files[i],sep="")
  tfrows=(read.delim(path, header=FALSE))[,c(1,2)]
  names(tfrows)=c("tf","gene")
  itfb=rbind(itfb,tfrows)
}
itfb["evidence"]="Predicted"
itfb["source"]="ITFB"


## =================================================================================================
## read the Oreganno interaction file  # we consider only the interactions with experiemental suport
## =================================================================================================
oreganno <- read.delim(config$oreganno.file)
oreganno=oreganno[oreganno$Species=="Homo sapiens",]
oreganno=oreganno[oreganno$Evidence.Subtypes!="UNKNOWN",]  
oreganno=oreganno[oreganno$Evidence.Subtypes!="N/A",] 
oreganno=oreganno[oreganno$Evidence.Subtypes!="Literature derived",] 
oreganno=oreganno[oreganno$TF.name!="UNKNOWN",]
oreganno=oreganno[oreganno$TF.name!="N/A",]
oreganno=oreganno[-agrep("UNKNOWN",oreganno$TF.name),]
oreganno["wordsinTFname"]=apply(as.data.frame(oreganno[,"TF.name"]),1,function(x)  length(unlist(strsplit(as.character(x)," "))))
oreganno=oreganno[oreganno$wordsinTFname<2,]
oreganno=oreganno[oreganno$Gene.name!="UNKNOWN",]
oreganno=oreganno[,c("TF.name","Gene.name")]
names(oreganno)=c("tf","gene")
oreganno["evidence"]="Experimental"
oreganno["source"]="Oreganno"



## ======================================
## read the PMID20584335 interaction file 
## ======================================
pmid <- read.delim(config$PMID20584335.file)
pmid=pmid[pmid$category=="TF-gene",]
pmid=pmid[! duplicated(pmid),]
pmid=pmid[,c(2:5)]
names(pmid)=c("tf","gene","evidence","source")



## ================================
## read the MSIGDB interaction file  (predicted)
## ================================
con <- file(config$msigdb.file, "r", blocking = FALSE, )
res=readLines(con) # empty
close(con)
msigdb=data.frame("gene"=character(),"tf"=character())
for(i in 1: length(res))
{
  s=res[i]
  strs=strsplit(s,"\t")
  tfbs=strs[[1]][1]
  df=data.frame(gene= strs[[1]][3:length(strs[[1]])])
  df["tf"]= strsplit(tfbs,"_") [[1]] [2]
  msigdb=rbind(msigdb,df)
}
msigdb=msigdb[,c(2,1)]
msigdb=msigdb[msigdb$tf !="UNKNOWN",]
msigdb["tf2"]= apply(as.data.frame(msigdb$tf),1, function(x)  substr(as.character(x), 3, nchar(as.character(x))))
msigdb=msigdb[,c(3,2)]
names(msigdb)=c("tf","gene")
msigdb["evidence"]="Predicted"
msigdb["source"]="MSigDB"




## ================================
## read the TRED interaction file  (predicted)
## ================================

con <- file(config$tred.file, "r", blocking = FALSE, )
res=readLines(con) # empty
close(con)
tred=data.frame("gene"=character(),"tf"=character())
for(i in 1: length(res))
{
  s=res[i]
  strs=strsplit(s,"\t")
  tfbs=strs[[1]][1]
  df=data.frame(gene= strs[[1]][3:length(strs[[1]])])
  df["tf"]=as.character(tfbs)
  tred=rbind(tred,df)
}
tred=tred[,c(2,1)]
tred=tred[- agrep( "http://" , tred$gene),]
names(tred)=c("tf","gene")
tred["evidence"]="Predicted"
tred["source"]="TRED"




# extra.dbs=rbind(tred,msigdb)
# extra.dbs["category"]="tf-gene"
# extra.dbs=extra.dbs[,c(5,1,2,3,4)]
# names(extra.dbs)=c("category","regulator","target","evidence","source")
# save(extra.dbs, file="temp.RDATA")






#tf.gene.db=rbind(transfac,itfb,oreganno,pmid,tred,msigdb)
tf.gene.db=rbind(transfac,oreganno,pmid,tred,msigdb)

d=tf.gene.db[,c(1,2)]
tf.gene.db=tf.gene.db[! duplicated(d),]
tf.gene.db["category"]="tf-gene"
tf.gene.db=tf.gene.db[,c(5,1,2,3,4)]

#save(tf.gene.db,file="temp.RDATA")

