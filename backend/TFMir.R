###########################################
#  - TFMir  Project                       #
#  - TFMir  main function                 #
#  - Main function to be called           #
#  - 2014-10-1                            #
#  - Copyright: Mohamed Hamed             #
###########################################

#R programming environments:
#- R studio IDE
#- R version 2.12.0 (2010-10-15)
#-Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
#-locale:   [1] C/en_US.UTF-8/C/C/C/C


## =============================================================================
## Initialize calling required packages and install them if they are not and log
## =============================================================================

source("loadpackages.R")


## ==================================================
## read the configuration file and get the parameters 
## ==================================================

source("readconfig.R")
config

## ===================================================
## Load the statistics ,color venn, and grapth scripts
## ===================================================
source("statistics.R")
source("ColorVenn.R")
source("graph.R")

## =============================================================================================================================
## Intialize and Read all files in config file (tranmir db , mirna disease association file and mirna function association file)
## =============================================================================================================================

source("initialize.R")
#getMIRnaCategory("hsa-mir-212","function")
#getMIRnaCategory("hsa-mir-212","disease")

## ===============================
## call the main function of TFMir
## ===============================

readInput=function(path)
{
  molecule.df = read.delim(path, header=FALSE)
  molecule.df=molecule.df[!duplicated (molecule.df),]
  molecule.input=unique(as.character(unlist(molecule.df[1]))) 
  return(molecule.input)
}


TFMir =function(tf.path, mirna.path,pval.cutoff=0.05,evidence,disease="",output.path)
{
  
#     tf.path="tf.sample2.txt"
#     mirna.path="mirna.sample.txt"
#     mirna.path=""
#     pval.cutoff=0.05
#     disease="Melanoma"
# #     disease="Neoplasm"
# #     disease="Alzheimer Disease"
# #     disease=""
#     output.folder="output"
#     evidence="Experimental"
#     evidence="Predicted"
#     evidence="both"
# 
#   
  ## ========================
  ## log the input parameters 
  ## ========================  
  writeToLog(tf.path)
  writeToLog(mirna.path)
  writeToLog(pval.cutoff)
  writeToLog(disease)
  writeToLog(evidence)
  writeToLog(output.path)
  writeToLog("=====================================")
  
  if(tolower(evidence)=="both") { evidence=c("Experimental","Predicted") }

  ## ==================================================================
  ## read the input files and test automatically which scenario will be
  ## ==================================================================
  tfs.df=data.frame("gene"=character(),"gene.reg"=numeric())
  if(! is.na(tf.path) & tf.path !="")
  {
    tfs.df = read.delim(tf.path, header=FALSE)    
  }else
  {
    mirnas.input=unique(tolower(readInput(mirna.path)))
    tf.pval=as.double(config$pval.cutoffgene.targets_regulators.fromMiRNA.inputlist)
    tfs.regulatorsofMiRNA= getTFsRegulatorsofMiRNAs(mirnas.input,tf.pval,evidence)
    tfs.targetsofMiRNA= getTFsTargetsofMiRNAs(mirnas.input,tf.pval,evidence)
    tfs.list=unique(c(tfs.targetsofMiRNA,tfs.regulatorsofMiRNA))
    if(length(tfs.list) > 0){
      tfs.df=data.frame("gene"=tfs.list,"gene.reg"=0) }
  }
  names(tfs.df)=c("gene","gene.reg")
  tfs.df=tfs.df[!duplicated (tfs.df),]
  tfs.input=toupper(unique(as.character(unlist(tfs.df$gene))))
  printGeneEntrezIDsMap(tfs.input,output.path)

  
  
  
  
  mirnas.df=data.frame("mirna"=character(),"mirna.reg"=numeric())
  if(! is.na(mirna.path) & mirna.path !="")
  {
    mirnas.df = read.delim(mirna.path, header=FALSE)    
  }else
  {
    
    mirna.pval=as.double(config$pval.cutoffmirna.regulators_targets.fromTFS.inputlist)
    tfstargets.mirna= getMiRNAsTargetsofTFs(tfs.input,mirna.pval,evidence)
    tfsregulators.mirna= getMiRNAsRegulatorsofTFs(tfs.input,mirna.pval,evidence)
    mirnas.list=unique(c(tfstargets.mirna,tfsregulators.mirna))
    if(length(mirnas.list) > 0){
    mirnas.df=data.frame("mirna"=mirnas.list,"mirna.reg"=0) }
  }
  names(mirnas.df)=c("mirna","mirna.reg")
  mirnas.df$mirna=tolower(mirnas.df$mirna)
  mirnas.df=mirnas.df[!duplicated (mirnas.df),]
  mirnas.input=unique(as.character(unlist(mirnas.df$mirna)))  
  
  
  
  ## ==================================
  ## get the four kinds of interactions 
  ## ==================================
  tf.mirna.res=getInteractions(category="tf-mirna",reg.input=tfs.input,target.input=mirnas.input,disease=disease,evidence=evidence,output.path=output.path,pval.cutoff=pval.cutoff)    
  tf.gene.res=getInteractions(category="tf-gene",reg.input=tfs.input,target.input=tfs.input,disease=disease,evidence=evidence,output.path=output.path,pval.cutoff=pval.cutoff)  
  #mirna.gene.res=tf.mirna.res[tf.mirna.res$category=="mirna-mirna",]  # return initial empty structure
  mirna.gene.res=getInteractions(category="mirna-gene",reg.input=mirnas.input,target.input=tfs.input,disease=disease,evidence=evidence,output.path=output.path,pval.cutoff=pval.cutoff)
  mirna.mirna.res=tf.mirna.res[tf.mirna.res$category=="mirna-mirna",]  # return initial empty structure
  if("Predicted" %in% evidence)   ## cause mirna-mirna interacctions are only predictions
    mirna.mirna.res=getInteractions(category="mirna-mirna",reg.input=mirnas.input,target.input=mirnas.input,disease=disease,evidence=evidence,output.path=output.path,pval.cutoff=pval.cutoff)
  
  ## ======================================================================================================================
  ## Combine these interactions and get those related to disease only (disease speccific network) (if disease is specified)
  ## ======================================================================================================================
  input=list( tf.genes=names(tf.gene.res), mirna.genes=names(mirna.gene.res),tf.mirna=names(tf.mirna.res),mirna.mirna=names(mirna.mirna.res))
  columns=Reduce(intersect,input)
  all.res=rbind(tf.mirna.res[,columns],mirna.mirna.res[,columns],mirna.gene.res[,columns],tf.gene.res[,columns])
  
  names(mirnas.df)=c("node","regulation")
  names(tfs.df)=c("node","regulation")
  nodes.input=rbind(mirnas.df,tfs.df)
  names(nodes.input)=c("target","target.reg")
  all.res=merge(all.res,nodes.input,by="target")
  names(nodes.input)=c("regulator","regulator.reg")
  all.res=merge(all.res,nodes.input,by="regulator")
  
  all.res.disease=all.res[all.res$is_regulator_in_disease==TRUE | all.res$is_target_in_disease==TRUE,]
  
  if(dim(all.res)[1] > 0)
    exportNetworkProperties (all.res,file.path(output.path,"all"), disease,pval.cutoff)
  
  if(dim(all.res.disease)[1] > 0)
    exportNetworkProperties (all.res.disease,file.path(output.path,"disease"),disease,pval.cutoff)
  
  write("finished", file=file.path(output.path,"finished.txt"),append=F,sep="\n")
  
}


ExportMotifs =function(net.path,output.path,evidence)
{
  ######### take care christian : for test only #####
  #net.path="output/disease/res.txt"
  #output.path="output/disease"
  ##################################################
  
  net=read.delim(net.path, header=TRUE) 
  if(tolower(evidence)=="both") { evidence=c("Experimental","Predicted") }
  
  #### extract all putative Tf-mirna paris who share target genes
  tfmir.pairs=getPutativeTFmiRPairs(net)
  #### extract all significant TF - miRNA pairs
  #tfmir.pairs=getSignificantTFmiRpairs(tfmir.pairs,evidence)
  
  if(dim(tfmir.pairs)[1] >0) 
  {
    #### relax and message the TF mir pairs who have more than one target
    tfmir.pairs=relaxMultipleTargetsForTFmiRPairs(tfmir.pairs)
    #### get motif type 1 : composite-FFL
    motifs.composite= getMotifs.composite(tfmir.pairs,net)
    #### get motif type 2 : TF-FFL
    motifs.TF.FFL= getMotifs.TF.FFL(tfmir.pairs,net)
    #### get motif type 3 : miRNA-FFL
    motifs.miRNA.FFL= getMotifs.miRNA.FFL(tfmir.pairs,net)
    #### get motif type 4 : Coregulation-FFL
    motifs.coregulation= getMotifs.coregulation(tfmir.pairs,net)
    motifs=rbind(motifs.composite,motifs.TF.FFL,motifs.miRNA.FFL,motifs.coregulation)
    if( dim(motifs)[1] > 0 )
    {
      motifs.ids=paste("motif",seq(1:dim(motifs)[1]), sep="" )
      motifs=cbind(motifs.ids,motifs)
      write.table(motifs,file=file.path(output.path,"motifs.txt"),quote=F,row.names=F,col.names=T,sep="\t")
    }
  }
  write("finished motifs", file=file.path(output.path,"finishedmotifs.txt"),append=F,sep="\n")  
}


PlotFunctionalSimilarity=function(genes,output.path)
{
###### @christian : these commented lines for testing only. u can try them urself
#   genes="ESR1, TP53, GRIN2D, AGER, AKT1, TERT, NCOA2, BBC3"
#   genes="CREB1, LTC4S, TLR9, IL5RA, MCAM, RPL10, RPS3A, ME2, CXCR4, SLC6A4, ERF, ID1, FLII, TGFB1, FLI1, UBE2I, PPRC1, CDC37, LRRFIP1, TGIF1, JAG1, TP53BP2, MSH6, MSH2"
#   genes="CREB1, RPL10, CXCR4, ID1, TGFB1, UBE2I, LRRFIP1, TGIF1, JAG1, TP53BP2, MSH6, MSH2"
#   genes="ESR1, TP53, GRIN2D, AGER, AKT1, TERT, NCOA2, BBC3"
#   genes="SPI1, BACH1, GNA13, SACM1L, FLI1, RAB23, POLE4, MSH2, SERTAD2, SKI, PHC2, ATP6V1C1, MSH6, DHX40, DPP7, RCN2, CHAF1A, PKN2, MECP2, ARL5B, MYO1E, B2M, TYROBP, FLII, MSR1, P2RY10, WAS"
#   genes="SPI1, BACH1, GNA13, SACM1L, FLI1, RAB23, POLE4, MSH2, SERTAD2, SKI, PHC2, ATP6V1C1, MSH6, DHX40, DPP7, RCN2, CHAF1A, PKN2, MECP2, ARL5B, MYO1E, B2M, TYROBP, FLII, MSR1, P2RY10, WAS"
#   output.path="output/disease/funsim.png"
#   ############################################
  
  print(output.path)

  genes=as.vector(unlist(strsplit(genes,",")))
  dput(genes)
  genes.entrez=unique(as.vector(unlist(mget(as.character(genes), envir=org.Hs.egALIAS2EG, ifnotfound=NA))))
  gosem=mgeneSim(genes.entrez,organism="human",measure="Wang")#,ont="BP"
  gosem=gosem[upper.tri(gosem)]
  
  all.entrez.genes <- mappedkeys(org.Hs.egACCNUM)
  #pvals.ks=c()
  #pvals.t=c()
  #pvals.wc=c()
  gosem.random.vector=c()
  for(i in 1: as.integer(config$NO_OF_random_permutations_for_functional_similarity))
  {
    genes.random=sample(all.entrez.genes,length(genes.entrez),replace = FALSE)
    gosem.random=mgeneSim(genes.random,organism="human",measure="Wang")#,ont="BP"
    gosem.random=gosem.random[upper.tri(gosem.random)]
    gosem.random.vector=c(gosem.random.vector,gosem.random)
    #     if(length(gosem.random)>1)
    #     {
    #       pvals.ks=c(pvals.ks,ks.test(gosem,gosem.random,alternative="l")$p.value)
    #       pvals.wc=c(pvals.wc,wilcox.test(gosem,gosem.random,alternative="g")$p.value)
    #       pvals.t=c(pvals.t,t.test(gosem,gosem.random,alternative="g")$p.value)      
    #     }
  }  
  #   pval.t.final= (length(pvals.t[pvals.t > 0.05]) /  length(pvals.t))
  #   pval.ks.final= (length(pvals.ks[pvals.ks > 0.05]) /  length(pvals.ks))
  #   pval.wc.final= (length(pvals.wc[pvals.wc > 0.05]) /  length(pvals.wc))
  #   pval=min(median(pvals.t),median(pvals.wc),median(pvals.ks))
  gosem.random.forplot=sample(gosem.random.vector,length(gosem))
  pval=ks.test(gosem,gosem.random.forplot,alternative="l")$p.value
  
  
  CairoPNG(bg="transparent",output.path,width=as.integer(config$funsimilarity.diagram.width),height=as.integer(config$funsimilarity.diagram.height))
  plot(ecdf(gosem),col="red", xlim=range(c(gosem, gosem.random.forplot)) , main="",xlab="Pair-wise similarity score", ylab="Cumulative distribution")
  #lines(ecdf(gosem),col="red",type="l")
  grid()
  lines(ecdf(gosem.random.forplot))
  #text(0.9,0.05, col="blue", paste("P-value < ",round(pval,4) ,sep=""),cex=1, adj = c(0.5, 0.5))
  text(0.9,0.05, col="blue", paste("P-value < ",format(pval, scientific = TRUE,digits=2) ,sep=""),cex=0.8, adj = c(0.5, 0.5))
  #mtext(paste("P-value < ",round(pval,3) ,sep=""), adj = 1,col="blue")
  legend(bty="n","topleft",c("Motif genes CDF","Random genes CDF") ,pch=c(19,19), col=c("red","black") ,cex=1)
  dev.off()
  
}


TFMir_old =function(tf.path, mirna.path,pval.cutoff=0.05,evidence,disease="",output.folder)
{
  
  #   tf.path="tf.sample2.txt"
  #   mirna.path="mirna.sample.txt"
  #   pval.cutoff=0.05
  #   disease="Melanoma"
  #   disease="Neoplasm"
  #   disease="Alzheimer Disease"
  #   disease=""
  #   output.folder="user8"
  #   evidence="Experimental"
  #   evidence="Predicted"
  #   evidence="both"
  #   
  
  ## ========================
  ## log the input parameters 
  ## ========================  
  writeToLog(tf.path)
  writeToLog(mirna.path)
  writeToLog(pval.cutoff)
  writeToLog(disease)
  writeToLog(evidence)
  writeToLog(output.folder)
  writeToLog("=====================================")
  
  ## ===========================
  ## Create the output directory 
  ## ===========================
  #   output.path= file.path(config$output.directory,output.folder)
  #   if  (! file.exists(output.path)){
  #     dir.create(file.path(output.path))
  #   }
  
  output.path= output.folder
  ## ================================================
  ## read the input files and intersect with transmir 
  ## ================================================
  tfs.df = read.delim(tf.path, header=FALSE)
  names(tfs.df)=c("gene","gene.reg")
  tfs.df=tfs.df[!duplicated (tfs.df),]
  tfs.input=toupper(unique(as.character(unlist(tfs.df$gene))))
  
  
  mirnas.df = read.delim(mirna.path, header=FALSE)
  names(mirnas.df)=c("mirna","mirna.reg")
  mirnas.df$mirna=tolower(mirnas.df$mirna)
  mirnas.df=mirnas.df[!duplicated (mirnas.df),]
  mirnas.input=unique(as.character(unlist(mirnas.df$mirna)))
  if(tolower(evidence)=="both") { evidence=c("Experimental","Predicted") }
  
  ## ==================================
  ## get the four kinds of interactions 
  ## ==================================
  tf.mirna.res=getInteractions(category="tf-mirna",reg.input=tfs.input,target.input=mirnas.input,disease=disease,evidence=evidence,output.path=output.path,pval.cutoff=pval.cutoff)  
  mirna.gene.res=getInteractions(category="mirna-gene",reg.input=mirnas.input,target.input=tfs.input,disease=disease,evidence=evidence,output.path=output.path,pval.cutoff=pval.cutoff)
  tf.gene.res=getInteractions(category="tf-gene",reg.input=tfs.input,target.input=tfs.input,disease=disease,evidence=evidence,output.path=output.path,pval.cutoff=pval.cutoff)  
  mirna.mirna.res=tf.mirna.res[tf.mirna.res$category=="mirna-mirna",]  # return initial empty structure
  if("Predicted" %in% evidence)   ## cause mirna-mirna interacctions are only predictions
    mirna.mirna.res=getInteractions(category="mirna-mirna",reg.input=mirnas.input,target.input=mirnas.input,disease=disease,evidence=evidence,output.path=output.path,pval.cutoff=pval.cutoff)
  
  ## ======================================================================================================================
  ## Combine these interactions and get those related to disease only (disease speccific network) (if disease is specified)
  ## ======================================================================================================================
  input=list( tf.genes=names(tf.gene.res), mirna.genes=names(mirna.gene.res),tf.mirna=names(tf.mirna.res),mirna.mirna=names(mirna.mirna.res))
  columns=Reduce(intersect,input)
  all.res=rbind(tf.mirna.res[,columns],mirna.mirna.res[,columns],mirna.gene.res[,columns],tf.gene.res[,columns])
  
  names(mirnas.df)=c("node","regulation")
  names(tfs.df)=c("node","regulation")
  nodes.input=rbind(mirnas.df,tfs.df)
  names(nodes.input)=c("target","target.reg")
  all.res=merge(all.res,nodes.input,by="target")
  names(nodes.input)=c("regulator","regulator.reg")
  all.res=merge(all.res,nodes.input,by="regulator")
  
  all.res.disease=all.res[all.res$is_regulator_in_disease==TRUE | all.res$is_target_in_disease==TRUE,]
  
  if(dim(all.res)[1] > 0)
    exportNetworkProperties (all.res,file.path(output.path,"all"), disease,pval.cutoff)
  
  if(dim(all.res.disease)[1] > 0)
    exportNetworkProperties (all.res.disease,file.path(output.path,"disease"),disease,pval.cutoff)
  
  write("finished", file=file.path(output.path,"finished.txt"),append=F,sep="\n")
  
}




# Let Graph G(V,E) be a connected graph, n=|V|, adj is the adjacency matrix of the graph G, and adj(i,i)=0, X is a binary array of size n, such that X(i)=1 if node I was marked as a key node, and 0 otherwise.
# Objective Function:  Min Σi  X(i)         Subjected to: ∀i  Σi adj(i,j) . X(j)>= 1
# 
# Let \ Graph \ G(V,E) \ be \ a \ connected\  graph, n=|V|,\\\  adj\  is \ the \ adjacency\  matrix \ of\  G, and \ adj(i,i)=0,\\  X is \ a\  binary\  array\  ofsize\  n, such\  that \ X(i)=1\  if \\ node \ i \ was \ marked \ as \ a \ key\  node, and \ 0 \ otherwise. \\ Then,\   The \ 
# Objective \ Function:  Min\sum _{i=1}^{ n} {X(i)}        \ \ \ \ \ \ \\       Subjected \ to:\forall i \sum_{i}^{n} adj(i,j)\ast X(j) >=1 


# # \textrm{Let Graph } G(V,E) \textrm{ be a connected graph, } n=|V|,  \\ \textrm{adj is the adjacency matrix of } G, \textrm{and adj}(i,i)=0,\\  X \textrm{ is a binary array of size } n, \textrm{such that }X(i)=1 \textrm{ if node}\\ i \textrm{ was marked as a key node, and 0  otherwise. Then, the}\\ 
# \textrm{Objective Function:} \\
# \min\sum _{i=1}^{ n} {X(i)}, \textrm{subjected to: }\forall i \sum_{i}^{n} \textrm{adj}(i,j)\cdot X(j) >=1 
# \\
# CR = \frac{N_d}{N_t}
# P = 1 - \sum^x_{i=0} \frac{\binom{k}{i}\binom{M-k}{N-i}}{\binom{M}{N}}\\
# Z = \frac{N_o - N_m}{\sigma}