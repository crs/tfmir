###########################################
#  - TFMir  Project                       #
#  - prepare DBs of miRNA-miRNA interact #
#  - 2014-10-1                            #
#  - Copyright: Mohamed Hamed             #
###########################################

#R programming environments:
#- R studio IDE
#- R version 2.12.0 (2010-10-15)
#-Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
#-locale:   [1] C/en_US.UTF-8/C/C/C/C



## ===============================
## read the PMMR database file 
## ===============================

pmmr <- read.delim(config$pmmr.file)
mirna.mirna.db=pmmr[pmmr$Score < 0.2,c(1,2)]    ## get the top 2% of the 
names(mirna.mirna.db)=c("from.mirna","to.mirna")
mirna.mirna.db["evidence"]="Predicted"
mirna.mirna.db["source"]="Pmmr"

mirna.mirna.db$from.mirna=trim(mirna.mirna.db$from.mirna)
mirna.mirna.db$to.mirna=trim(mirna.mirna.db$to.mirna)

mirna.mirna.db=mirna.mirna.db[! duplicated(mirna.mirna.db), ]
mirna.mirna.db["category"]="mirna-mirna"
mirna.mirna.db=mirna.mirna.db[,c(5,1,2,3,4)]

# 
# # 
# names(tf.gene.db)=c("category","regulator","target","evidence","source")
# names(tf.mirna.db)=c("category","regulator","target","evidence","source")
# names(mirna.mirna.db)=c("category","regulator","target","evidence","source")
# names(mirna.gene.db)=c("category","regulator","target","evidence","source")
# dbs=rbind(tf.gene.db,tf.mirna.db,mirna.mirna.db,mirna.gene.db)
# 
# save(dbs,tf.mirna.db,mirna.gene.db,tf.gene.db,mirna.mirna.db,file="databases.RDATA")


##### correct some issues in the dbs  that was saveed  in databases.RDATA file
# 1- correct the ] and [ in the mirna-genes interactions
dbs.part1= dbs[- grep("\\[",dbs$regulator), ]
x=dbs[grep("\\[",dbs$regulator), ]
x$regulator=gsub("\\[","",x$regulator)
x$regulator=gsub("\\]","",x$regulator)
x$regulator=gsub("has","hsa",x$regulator)
x[- grep("hsa-",x$regulator),]$regulator=paste("hsa-",x[- grep("hsa-",x$regulator),]$regulator,sep="")
dbs=rbind(dbs.part1,x)
# 2- correct the empty spaces before and after
dbs$regulator=trim(dbs$regulator)
dbs$target=trim(dbs$target)
# 3- remove duplicates
dbs=dbs[! duplicated(dbs),]
d=dbs[,c(2,3,4)]
dbs=dbs[! duplicated(d),]
#4-remove MSIGDB and keep only tred...... if MSIgDB there it returns many links 
dbs=dbs[! (dbs$category=="tf-gene" & dbs$source=="MSigDB"),]
#save(dbs,tf.mirna.db,mirna.gene.db,tf.gene.db,mirna.mirna.db,file="databases2.RDATA")

