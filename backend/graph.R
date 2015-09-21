###########################################
#  - TFMir  Project                       #
#  - TFMir  Graph code statistics         #
#  - Main function to be called           #
#  - 2014-10-1                            #
#  - Copyright: Mohamed Hamed             #
###########################################

#R programming environments:
#- R studio IDE
#- R version 2.12.0 (2010-10-15)
#-Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
#-locale:   [1] C/en_US.UTF-8/C/C/C/C


## ====================================
## analyze the disease specific network
## ====================================
getEdgeIndexforDisease=function(net)
{
  return( (dim(net[net$is_regulator_in_disease==TRUE | net$is_target_in_disease==TRUE,])[1]) / (dim(net)[1] ) * 100)
}

getNodeIndexforDisease=function(net)
{
  nodes.all.net=unique(c(as.character(net$regulator) , as.character(net$target) ))
  dis.net=net[net$is_regulator_in_disease==TRUE | net$is_target_in_disease==TRUE,]
  nodes.dis.net= unique(c(as.character(dis.net$regulator) , as.character(dis.net$target) ))
  index= (length(nodes.dis.net) / length(nodes.all.net)) *100
  return(index)
}


plotDegreeDistribution_old=function (net,disease.net,output.path,disease)
{
  #### the whole network
  nodes.df=as.data.frame( unique(c(as.character(unlist(net$regulator)), as.character(unlist(net$target)))))  
  g=graph.data.frame(net,directed=T,vertices=nodes.df)
  degree.dist.net= degree.distribution(g, cumulative = T)
  
  #### the disease network
  nodes.df=as.data.frame( unique(c(as.character(unlist(disease.net$regulator)), as.character(unlist(disease.net$target)))))  
  g=graph.data.frame(disease.net,directed=T,vertices=nodes.df)
  degree.dist.disease.net= degree.distribution(g, cumulative = T)
  
  ## ==========================================
  ## Degree distribution of the both network
  ## ==========================================
  CairoPNG(output.path,width=as.integer(config$degreediagram.width),height=as.integer(config$degreediagram.hight))
  plot(degree.dist.net,log="xy", col=2, pch=20,type="p",main="",xlab="Degree", ylab="Frequency")
  if(! is.na(disease) & disease !=""){
    points(degree.dist.disease.net,col=3,pch=18)
    legend("bottomleft", c("Whole network links ","Disease specific network"),cex=0.75 ,col=2:3, pch=c(20,18), ncol=1, yjust=0, lty=0)
  }
  dev.off()
}

plotDegreeDistribution=function (net,output.path)
{
  nodes.df=as.data.frame( unique(c(as.character(unlist(net$regulator)), as.character(unlist(net$target)))))  
  g=graph.data.frame(net,directed=T,vertices=nodes.df)
  degree.dist.net= degree.distribution(g, cumulative = T)
  ## ==========================================
  ## Degree distribution of the network
  ## ==========================================
  CairoPNG(output.path,width=as.integer(config$degreediagram.width),height=as.integer(config$degreediagram.hight))
  plot(degree.dist.net,log="xy", col=2, pch=20,type="p",main="",xlab="Degree", ylab="Frequency")
  dev.off()
}

getnetworkProperties=function (df)
{
  
  nodes.df=as.data.frame( unique(c(as.character(unlist(df$regulator)), as.character(unlist(df$target)))))  
  g=graph.data.frame(df,directed=T,vertices=nodes.df)
  ## ====================================
  ## Network measures and some properites
  ## ====================================
  graph.nodes.no=length(V(g))
  graph.edges.no=length(E(g))
  graph.density=round(graph.density(g, loops=T),digits=3)
  graph.diameter=diameter(g)
  graph.avg.path.length=  round( average.path.length(g) , digits=3)
  graph.transitivity= round (transitivity(g) , digits=3)
  
  #plot(g,vertex.label = V(g)$name,layout=layout.fruchterman.reingold(g),edge.arrow.size =0.4)
  #get the hot spots from the top 5 % of the highest degree nodes 
  index=rev(order(degree(g)))
  degree.hotspots=V(g)$name[index] [1:ceiling( length(V(g)) * as.integer(config$hotspot.percentage) /100) ]
  
  index=rev(order(closeness(g)))
  closeness.hotspots=V(g)$name[index] [1:ceiling( length(V(g)) * as.integer(config$hotspot.percentage) /100) ]
  
  
  index=rev(order(betweenness(g)))
  betweenness.hotspots=V(g)$name[index] [1:ceiling( length(V(g)) * as.integer(config$hotspot.percentage) /100) ]
  
  eigenvector.hotspots=c("")
  if(length(E(g)) > 1 ){
  index=rev(order(evcent(g)$vector))
  eigenvector.hotspots=V(g)$name[index] [1:ceiling( length(V(g)) * as.integer(config$hotspot.percentage) /100) ]
  }
    
  common.hotspots=unique(Reduce(intersect,  list(v1 = degree.hotspots, 
                                v2 = closeness.hotspots, 
                                v3 = betweenness.hotspots,
                                v4 = eigenvector.hotspots)))
  union.hotspots=unique(c(degree.hotspots,closeness.hotspots,betweenness.hotspots,eigenvector.hotspots))
  
  return(c( graph.nodes.no,graph.edges.no,graph.density,graph.diameter,graph.avg.path.length,graph.transitivity,
            toString(degree.hotspots),toString(closeness.hotspots),toString(betweenness.hotspots),
            toString(eigenvector.hotspots),toString(common.hotspots),toString(union.hotspots)))
}




exportNetworkProperties=function (net,net.output.path,disease,pval.cutoff)
{
  
  if  (! file.exists(net.output.path)){
    dir.create(file.path(net.output.path))
  }
  write.table(net,file=file.path(net.output.path,"res.txt"), quote=F,row.names=F,col.names=T,sep="\t")
  
  
  ## =======================================================================
  ## DO functional, miRNA enrichment, and statistics for the final node list
  ## =======================================================================
  allnodes.disease=c()
  if(! is.na(disease) & disease !="")
    {
      allnodes.disease=c( as.character(getGenesforDisease(disease)),as.character(getmiRNAforDisease(disease)))  
    }
  
  nodes.net= unique(c(as.character(net$regulator) , as.character(net$target) ))
  pval.hypergeom.disease.node=overlapSignificance_By_HyperGEOM(total = as.integer(config$Total.No.of.miRNA.in.human) + as.integer(config$Total.No.of.genes.in.human),
                                                               numgA=length(nodes.net),numgB=length(allnodes.disease),
                                                               overlap=length(intersect(nodes.net,allnodes.disease)))
  
  nodes.net.mirna=  nodes.net[ grep("hsa-",nodes.net)]
  nodes.net.gene=  setdiff(nodes.net,nodes.net.mirna)
  nodes.net.gene.etrezID=toString(as.vector(unlist(mget(as.character(nodes.net.gene), envir=org.Hs.egALIAS2EG, ifnotfound=NA))))
  david.BP=paste("http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=",nodes.net.gene.etrezID,",&tool=chartReport&annot=GOTERM_BP_ALL",sep="")
  david.KEGG=paste("http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=",nodes.net.gene.etrezID,",&tool=chartReport&annot=KEGG_PATHWAY",sep="")
  david.OMIM=paste("http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=",nodes.net.gene.etrezID,",&tool=chartReport&annot=OMIM_DISEASE",sep="")    
  david.functional.clust=paste("http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=",nodes.net.gene.etrezID,",&tool=term2term&annot=GOTERM_BP_ALL",sep="")    
  
  if(length(nodes.net.mirna) > 0)
  {
    mirna.ora.function=DO_ORA_FOR_MIRNA(nodes.net.mirna,category="function",pval.cutoff=pval.cutoff)
    mirna.ora.disease=DO_ORA_FOR_MIRNA(nodes.net.mirna,category="disease",pval.cutoff = pval.cutoff)
    mirna.ora=rbind(mirna.ora.function, mirna.ora.disease)
    write.table(mirna.ora,file=file.path(net.output.path,"mirna.ora.txt"), quote=F,row.names=F,col.names=T,sep="\t")
  }
  
  if(length(nodes.net.gene) > 0)
  {
    summary.file=file.path (net.output.path,"genes.ora.txt")
    write(nodes.net.gene.etrezID,      file=summary.file,append=F,sep="\n")
    
    #write(paste("david.BP.link=",david.BP,sep=""),      file=summary.file,append=F,sep="\n")
    #write(paste("david.KEGG.link=",david.KEGG,sep=""),      file=summary.file,append=T,sep="\n")
    #write(paste("david.OMIM.link=",david.OMIM,sep=""),      file=summary.file,append=T,sep="\n")
    #write(paste("david.functional.clust.link=",david.functional.clust,sep=""),      file=summary.file,append=T,sep="\n")
  }
  
  
  ## ==============================================================
  ## Perform network analysis and statistics on the disease network 
  ## ==============================================================
  disease.net.coverage.rate.edge=round(getEdgeIndexforDisease(net),digits=2)
  disease.net.coverage.rate.node=round(getNodeIndexforDisease(net),digits=2)  
  degree.image.path=file.path(net.output.path,"degree.png")
  plotDegreeDistribution(net,degree.image.path)
  network.properties=getnetworkProperties(net)
  
  ## =========================================
  ## output the summary of the disease network
  ## =========================================
  
  
  summary.file=file.path (net.output.path,"summary.txt")
  
  write(paste("pval.hypergeom.disease.node=",pval.hypergeom.disease.node,sep=""),      file=summary.file,append=F,sep="\n")
  write(paste("disease.net.coverage.rate.edge=",disease.net.coverage.rate.edge,sep=""),      file=summary.file,append=T,sep="\n")
  write(paste("disease.net.coverage.rate.node=",disease.net.coverage.rate.node,sep=""),  file=summary.file,append=T,sep="\n")
  write(paste("degree.image.path=",degree.image.path,sep=""),      file=summary.file,append=T,sep="\n")
  
  write(paste("graph.nodes.no=",network.properties[1],sep=""),      file=summary.file,append=T,sep="\n")
  write(paste("graph.edges.no=",network.properties[2],sep=""),      file=summary.file,append=T,sep="\n")
  write(paste("graph.density=",network.properties[3],sep=""),      file=summary.file,append=T,sep="\n")
  write(paste("graph.diameter=",network.properties[4],sep=""),      file=summary.file,append=T,sep="\n")
  write(paste("graph.avg.path.length=",network.properties[5],sep=""),      file=summary.file,append=T,sep="\n")
  write(paste("graph.transitivity=",network.properties[6],sep=""),      file=summary.file,append=T,sep="\n")
  #http://horicky.blogspot.de/2012/04/basic-graph-analytics-using-igraph.html
  write(paste("degree.hotspots=",network.properties[7],sep=""),      file=summary.file,append=T,sep="\n")
  write(paste("closeness.hotspots=",network.properties[8],sep=""),      file=summary.file,append=T,sep="\n")
  write(paste("betweenness.hotspots=",network.properties[9],sep=""),      file=summary.file,append=T,sep="\n")
  write(paste("eigenvector.hotspots=",network.properties[10],sep=""),      file=summary.file,append=T,sep="\n")
  write(paste("common.hotspots=",network.properties[11],sep=""),      file=summary.file,append=T,sep="\n")
  write(paste("union.hotspots=",network.properties[12],sep=""),      file=summary.file,append=T,sep="\n")
  
}




## =======================
## Motif detection methods 
## =======================

getPutativeTFmiRPairs=function (net)
{
  regulators=unique(as.character(unlist(net$regulator)))
  mirnas=regulators[ grep("hsa-",regulators)]
  tfs=setdiff(regulators, mirnas)
  
  tfmir.pairs=data.frame(tf=character(), mirna=character(),targets=character(),pvals.tfmirpair=numeric())
  
  for(i in 1 : length(tfs))
  {
    for(j in 1 :length(mirnas))
    {
      targets.mirna= unlist(net[net$regulator ==mirnas[j], ]$target)
      targets.mirna= setdiff(targets.mirna, targets.mirna [ grep("hsa-", targets.mirna)])
      
      targets.tf= unlist(net[net$regulator ==tfs[i], ]$target)
      targets.tf= setdiff(targets.tf, targets.tf [ grep("hsa-", targets.tf)])
      
      
      common.targets=intersect(targets.mirna,targets.tf )
      if( length(common.targets) > 0 )
      {
        row=data.frame (tf=tfs[i], mirna=mirnas[j], targets=toString(common.targets),pvals.tfmirpair=0)
        tfmir.pairs=rbind(tfmir.pairs,row)
      }
    }
  }
  
  return(tfmir.pairs)
}

getSignificantTFmiRpairs=function (net,evidence)
{
  
  ### total is the number of common genes between all human genes targetted by human miRNA and all hman genes targeted by TFs
  targets.bymirna=unique(as.character(unlist(dbs[dbs$category=="mirna-gene" & dbs$evidence %in% evidence, ]$target)))
  targets.tfs=unique(as.character(unlist(dbs[dbs$category=="tf-gene" & dbs$evidence %in% evidence, ]$target)))
  total=length(intersect(targets.bymirna, targets.tfs))
  
  pvals=c()
  for(i in 1 : dim(net)[1])
  {
    tf=as.character(net$tf[i])
    targets.tf= unique(as.character(unlist(dbs[dbs$category=="tf-gene" & dbs$evidence %in% evidence & dbs$regulator == tf, ]$target)))  
    mirna=as.character(net$mirna[i])
    targets.mirna= unique(as.character(unlist(dbs[dbs$category=="mirna-gene" & dbs$evidence %in% evidence & dbs$regulator == mirna, ]$target)))  
    overlap=length(intersect(tolower(targets.mirna),tolower(targets.tf)))
    numgA=length(targets.tf)
    numgB=length(targets.mirna)
    pvals=c(pvals,phyper(overlap-1, numgA, total - numgA, numgB,lower.tail=FALSE) )       
  }
  net["pvals.tfmirpair"]=pvals
#   pvals.BH=mt.rawp2adjp(pvals,proc="BH")
#   net=net[pvals.BH$index,]
#   pvals.BH=(pvals.BH$adjp)[,2]
#   net["pvals.BH"]=pvals.BH
  net=net[net$pvals < as.double(config$pval.cutoff.significant.TF.MIR.PAIRS),]
  return(net)
}


relaxMultipleTargetsForTFmiRPairs=function(tfmir.pairs)
{
    tfmir.pairs.new=tfmir.pairs[tfmir.pairs$tf=="SISI IBn KELAB",]  
    for(i in 1: dim(tfmir.pairs)[1])
    {
      s=as.character(tfmir.pairs[i,]$targets)
      strs=strsplit(s,", ")
      df=data.frame(targets= strs[[1]][1:length(strs[[1]])])
      df["mirna"]=tfmir.pairs[i,]$mirna
      df["tf"]=tfmir.pairs[i,]$tf
      df["pvals.tfmirpair"]=tfmir.pairs[i,]$pvals.tfmirpair
      df=df[,c(3,2,1,4)]
      tfmir.pairs.new=rbind(tfmir.pairs.new,df)
    }
  tfmir.pairs.new
}


getMotifs.composite=function (tfmir.pairs,net)
{ 
  motifs.empty=data.frame("motif.type"=character(), "tf"=character(),"mirna"=character(),"targets"=character(), "pvals.tfmirpair"=numeric(),
             "pval.random"=numeric(),"zscore"=numeric(),"ENTREZ.IDs.coreg.targets"=character(),
             "ENTREZ.IDs.motiftype.subnetwork"=character())
  index=c()
  for (i  in 1 : dim(tfmir.pairs)[1])
  {
    tf   = as.character( tfmir.pairs$tf[i])
    mirna= as.character( tfmir.pairs$mirna[i] )
    ## check for composite motifs condition
    if( dim(net[net$regulator == tf & net$target == mirna, ])[1] > 0   & dim(net[net$regulator == mirna & net$target == tf, ])[1] > 0 )
    {index=c(index,i)}
  }
  motifs=tfmir.pairs[index,]
  motifs
  
  if(dim(motifs)[1] > 0)
  {
  ######## Randomize the network and get significance of this type
  nodes.df=as.data.frame( unique(c(as.character(unlist(net$regulator)), as.character(unlist(net$target)))))  
  g=graph.data.frame(net[,c("regulator","target")],directed=T,vertices=nodes.df)
  tmpres=c()
  for( i in 1 : as.integer(config$NO_OF_Times_FOR_Network_Randomizations))
  {
    #g.random=rewire(g, mode = "simple", niter = 2*length(E(g)))
    g.random=rewire(g, keeping_degseq(loops=TRUE, niter = 2*length(get.edgelist(g))))
    net.random=as.data.frame(get.edgelist(g.random,names=T))
    names(net.random)=c("regulator","target")
    tfmir.pairs.random=getPutativeTFmiRPairs(net.random)    
    #tfmir.pairs.random=getSignificantTFmiRpairs(tfmir.pairs.random,evidence)
    tfmir.pairs.random["pvals.tfmirpair"]=1
    tfmir.pairs.random=relaxMultipleTargetsForTFmiRPairs(tfmir.pairs.random)
    index=c()
    for ( j  in 1 : dim(tfmir.pairs.random)[1])
    {
      tf   = as.character( tfmir.pairs.random$tf[j])
      mirna= as.character( tfmir.pairs.random$mirna[j] )
      ## check for composite motifs condition
      if( dim(net.random[net.random$regulator == tf & net.random$target == mirna, ])[1] > 0   & dim(net.random[net.random$regulator == mirna & net.random$target == tf, ])[1] > 0 )
      {index=c(index,j)}
    }
    motifs.random=tfmir.pairs.random[index,]
    tmpres= c(tmpres,dim(motifs.random)[1])    
  }
  pval<-mean(tmpres > dim(motifs)[1])
  zscore=   ((dim(motifs)[1]) - (mean(tmpres))) / sd(tmpres)
  motifs["pval.random"]=pval
  motifs["zscore"]=zscore
  ##################################################################
  ###### get the entrez ids of co-regulated taget genes of each motif
  motifs["Symbol.coreg.targets"]=""
  motifs["ENTREZ.IDs.coreg.targets"]=""
  for( i in 1 : dim(motifs)[1])
  {
    tfmirpair=c(as.character(motifs[i,]$tf), as.character(motifs[i,]$mirna))
    targets=as.character( unlist(net[net$regulator %in% tfmirpair,]$target))
    targets.mir=targets[grep("hsa-",targets)]
    targets.gene=unique(c(as.character(motifs[i,]$tf), setdiff(targets,targets.mir)))
    motifs[i, ]$Symbol.coreg.targets=toString(targets.gene)
    motifs[i, ]$ENTREZ.IDs.coreg.targets=toString(as.vector(unlist(mget(as.character(targets.gene), envir=org.Hs.egALIAS2EG, ifnotfound=NA))))
  }
  
  ###### get the entrez ids of co-targetted genes by the same Tf and miRNA
  motifs["Symbol.cotargeted.genes"]=""
  motifs["ENTREZ.IDs.cotargeted.genes"]=""
  for( i in 1 : dim(motifs)[1])
  {
    targets.mirna= unlist(net[net$regulator == as.character(motifs[i,]$mirna) , ]$target)
    targets.mirna= setdiff(targets.mirna, targets.mirna [ grep("hsa-", targets.mirna)])  
    targets.tf= unlist(net[net$regulator == as.character(motifs[i,]$tf), ]$target)
    targets.tf= setdiff(targets.tf, targets.tf [ grep("hsa-", targets.tf)])
    common.targets=unique(intersect(targets.mirna,targets.tf))
    motifs[i, ]$Symbol.cotargeted.genes=toString(common.targets)
    motifs[i, ]$ENTREZ.IDs.cotargeted.genes=toString(as.vector(unlist(mget(as.character(common.targets), envir=org.Hs.egALIAS2EG, ifnotfound=NA))))
  }
  #   ###### get the entrez ids of subnetworks of each motif type
  #   tfs= as.character(unlist(motifs$tf))
  #   genes=as.character(unlist(motifs$targets))
  #   tfs.genes=unique(c(tfs,genes))
  #   tfs.genes.entrez=as.vector(unlist(mget(as.character(tfs.genes), envir=org.Hs.egALIAS2EG, ifnotfound=NA)))
  #   motifs["ENTREZ.IDs.motiftype.subnetwork"]=toString(tfs.genes.entrez)
  ######  add the type of the motif
  motif.type=replicate(dim(motifs)[1],"COMPOSITE-FFL")
  motifs=cbind(motif.type,motifs)
  return (motifs)
  }
  
  return (motifs.empty)
}

getMotifs.TF.FFL=function (tfmir.pairs,net)
{ 
  motifs.empty=data.frame("motif.type"=character(), "tf"=character(),"mirna"=character(),"targets"=character(), "pvals.tfmirpair"=numeric(),
                          "pval.random"=numeric(),"zscore"=numeric(),"ENTREZ.IDs.coreg.targets"=character(),
                          "ENTREZ.IDs.motiftype.subnetwork"=character())
  index=c()
  for (i  in 1 : dim(tfmir.pairs)[1])
  {
    tf   = as.character( tfmir.pairs$tf[i])
    mirna= as.character( tfmir.pairs$mirna[i] )
    ## check for composite motifs condition
    if( dim(net[net$regulator == tf & net$target == mirna, ])[1] > 0   & dim(net[net$regulator == mirna & net$target == tf, ])[1] == 0 )
    {index=c(index,i)}
  }
  motifs=tfmir.pairs[index,]
  motifs
  
  if(dim(motifs)[1] > 0)
  {
    ######## Randomize the network and get significance of this type
    nodes.df=as.data.frame( unique(c(as.character(unlist(net$regulator)), as.character(unlist(net$target)))))  
    g=graph.data.frame(net[,c("regulator","target")],directed=T,vertices=nodes.df)
    tmpres=c()
    for( i in 1 : as.integer(config$NO_OF_Times_FOR_Network_Randomizations))
    {
      g.random=rewire(g, keeping_degseq(loops=TRUE, niter = 2*length(get.edgelist(g))))
      net.random=as.data.frame(get.edgelist(g.random,names=T))
      names(net.random)=c("regulator","target")
      tfmir.pairs.random=getPutativeTFmiRPairs(net.random)    
      #tfmir.pairs.random=getSignificantTFmiRpairs(tfmir.pairs.random,evidence)
      tfmir.pairs.random["pvals.tfmirpair"]=1
      tfmir.pairs.random=relaxMultipleTargetsForTFmiRPairs(tfmir.pairs.random)
      index=c()
      for ( j  in 1 : dim(tfmir.pairs.random)[1])
      {
        tf   = as.character( tfmir.pairs.random$tf[j])
        mirna= as.character( tfmir.pairs.random$mirna[j] )
        ## check for composite motifs condition
        if( dim(net.random[net.random$regulator == tf & net.random$target == mirna, ])[1] > 0   & dim(net.random[net.random$regulator == mirna & net.random$target == tf, ])[1] == 0 )
        {index=c(index,j)}
      }
      motifs.random=tfmir.pairs.random[index,]
      tmpres= c(tmpres,dim(motifs.random)[1])    
    }
    pval<-mean(tmpres > dim(motifs)[1])
    zscore=   ((dim(motifs)[1]) - (mean(tmpres))) / sd(tmpres)
    motifs["pval.random"]=pval
    motifs["zscore"]=zscore
    ##################################################################
    ###### get the entrez ids of co-regulated taget genes of each motif
    motifs["Symbol.coreg.targets"]=""
    motifs["ENTREZ.IDs.coreg.targets"]=""
    for( i in 1 : dim(motifs)[1])
    {
      tfmirpair=c(as.character(motifs[i,]$tf), as.character(motifs[i,]$mirna))
      targets=as.character( unlist(net[net$regulator %in% tfmirpair,]$target))
      targets.mir=targets[grep("hsa-",targets)]
      targets.gene=unique(c(as.character(motifs[i,]$tf), setdiff(targets,targets.mir)))
      motifs[i, ]$Symbol.coreg.targets=toString(targets.gene)
      motifs[i, ]$ENTREZ.IDs.coreg.targets=toString(as.vector(unlist(mget(as.character(targets.gene), envir=org.Hs.egALIAS2EG, ifnotfound=NA))))
    }
    
    ###### get the entrez ids of co-targetted genes by the same Tf and miRNA
    motifs["Symbol.cotargeted.genes"]=""
    motifs["ENTREZ.IDs.cotargeted.genes"]=""
    for( i in 1 : dim(motifs)[1])
    {
      targets.mirna= unlist(net[net$regulator == as.character(motifs[i,]$mirna) , ]$target)
      targets.mirna= setdiff(targets.mirna, targets.mirna [ grep("hsa-", targets.mirna)])  
      targets.tf= unlist(net[net$regulator == as.character(motifs[i,]$tf), ]$target)
      targets.tf= setdiff(targets.tf, targets.tf [ grep("hsa-", targets.tf)])
      common.targets=unique(intersect(targets.mirna,targets.tf))
      motifs[i, ]$Symbol.cotargeted.genes=toString(common.targets)
      motifs[i, ]$ENTREZ.IDs.cotargeted.genes=toString(as.vector(unlist(mget(as.character(common.targets), envir=org.Hs.egALIAS2EG, ifnotfound=NA))))
    }
    #   ###### get the entrez ids of subnetworks of each motif type
    #   tfs= as.character(unlist(motifs$tf))
    #   genes=as.character(unlist(motifs$targets))
    #   tfs.genes=unique(c(tfs,genes))
    #   tfs.genes.entrez=as.vector(unlist(mget(as.character(tfs.genes), envir=org.Hs.egALIAS2EG, ifnotfound=NA)))
    #   motifs["ENTREZ.IDs.motiftype.subnetwork"]=toString(tfs.genes.entrez)
    ######  add the type of the motif
    motif.type=replicate(dim(motifs)[1],"TF-FFL")
    motifs=cbind(motif.type,motifs)
    return (motifs)
  }
  
  return (motifs.empty)
}

getMotifs.miRNA.FFL=function (tfmir.pairs,net)
{ 
  motifs.empty=data.frame("motif.type"=character(), "tf"=character(),"mirna"=character(),"targets"=character(), "pvals.tfmirpair"=numeric(),
                          "pval.random"=numeric(),"zscore"=numeric(),"ENTREZ.IDs.coreg.targets"=character(),
                          "ENTREZ.IDs.motiftype.subnetwork"=character())
  index=c()
  for (i  in 1 : dim(tfmir.pairs)[1])
  {
    tf   = as.character( tfmir.pairs$tf[i])
    mirna= as.character( tfmir.pairs$mirna[i] )
    ## check for composite motifs condition
    if( dim(net[net$regulator == tf & net$target == mirna, ])[1] == 0   & dim(net[net$regulator == mirna & net$target == tf, ])[1] > 0 )
    {index=c(index,i)}
  }
  motifs=tfmir.pairs[index,]
  motifs
  
  if(dim(motifs)[1] > 0)
  {
    ######## Randomize the network and get significance of this type
    nodes.df=as.data.frame( unique(c(as.character(unlist(net$regulator)), as.character(unlist(net$target)))))  
    g=graph.data.frame(net[,c("regulator","target")],directed=T,vertices=nodes.df)
    tmpres=c()
    for( i in 1 : as.integer(config$NO_OF_Times_FOR_Network_Randomizations))
    {
      g.random=rewire(g, keeping_degseq(loops=TRUE, niter = 2*length(get.edgelist(g))))
      net.random=as.data.frame(get.edgelist(g.random,names=T))
      names(net.random)=c("regulator","target")
      tfmir.pairs.random=getPutativeTFmiRPairs(net.random)    
      #tfmir.pairs.random=getSignificantTFmiRpairs(tfmir.pairs.random,evidence)
      tfmir.pairs.random["pvals.tfmirpair"]=1
      tfmir.pairs.random=relaxMultipleTargetsForTFmiRPairs(tfmir.pairs.random)
      index=c()
      for ( j  in 1 : dim(tfmir.pairs.random)[1])
      {
        tf   = as.character( tfmir.pairs.random$tf[j])
        mirna= as.character( tfmir.pairs.random$mirna[j] )
        ## check for composite motifs condition
        if( dim(net.random[net.random$regulator == tf & net.random$target == mirna, ])[1] == 0   & dim(net.random[net.random$regulator == mirna & net.random$target == tf, ])[1] > 0 )
        {index=c(index,j)}
      }
      motifs.random=tfmir.pairs.random[index,]
      tmpres= c(tmpres,dim(motifs.random)[1])    
    }
    pval<-mean(tmpres > dim(motifs)[1])
    zscore=   ((dim(motifs)[1]) - (mean(tmpres))) / sd(tmpres)
    motifs["pval.random"]=pval
    motifs["zscore"]=zscore
    ##################################################################
    ###### get the entrez ids of co-regulated taget genes of each motif
    motifs["Symbol.coreg.targets"]=""
    motifs["ENTREZ.IDs.coreg.targets"]=""
    for( i in 1 : dim(motifs)[1])
    {
      tfmirpair=c(as.character(motifs[i,]$tf), as.character(motifs[i,]$mirna))
      targets=as.character( unlist(net[net$regulator %in% tfmirpair,]$target))
      targets.mir=targets[grep("hsa-",targets)]
      targets.gene=unique(c(as.character(motifs[i,]$tf), setdiff(targets,targets.mir)))
      motifs[i, ]$Symbol.coreg.targets=toString(targets.gene)
      motifs[i, ]$ENTREZ.IDs.coreg.targets=toString(as.vector(unlist(mget(as.character(targets.gene), envir=org.Hs.egALIAS2EG, ifnotfound=NA))))
    }
    
    ###### get the entrez ids of co-targetted genes by the same Tf and miRNA
    motifs["Symbol.cotargeted.genes"]=""
    motifs["ENTREZ.IDs.cotargeted.genes"]=""
    for( i in 1 : dim(motifs)[1])
    {
      targets.mirna= unlist(net[net$regulator == as.character(motifs[i,]$mirna) , ]$target)
      targets.mirna= setdiff(targets.mirna, targets.mirna [ grep("hsa-", targets.mirna)])  
      targets.tf= unlist(net[net$regulator == as.character(motifs[i,]$tf), ]$target)
      targets.tf= setdiff(targets.tf, targets.tf [ grep("hsa-", targets.tf)])
      common.targets=unique(intersect(targets.mirna,targets.tf))
      motifs[i, ]$Symbol.cotargeted.genes=toString(common.targets)
      motifs[i, ]$ENTREZ.IDs.cotargeted.genes=toString(as.vector(unlist(mget(as.character(common.targets), envir=org.Hs.egALIAS2EG, ifnotfound=NA))))
    }
    #   ###### get the entrez ids of subnetworks of each motif type
    #   tfs= as.character(unlist(motifs$tf))
    #   genes=as.character(unlist(motifs$targets))
    #   tfs.genes=unique(c(tfs,genes))
    #   tfs.genes.entrez=as.vector(unlist(mget(as.character(tfs.genes), envir=org.Hs.egALIAS2EG, ifnotfound=NA)))
    #   motifs["ENTREZ.IDs.motiftype.subnetwork"]=toString(tfs.genes.entrez)
    ######  add the type of the motif
    motif.type=replicate(dim(motifs)[1],"miRNA-FFL")
    motifs=cbind(motif.type,motifs)
    return (motifs)
  }
  
  return (motifs.empty)
}

getMotifs.coregulation=function (tfmir.pairs,net)
{ 
  motifs.empty=data.frame("motif.type"=character(), "tf"=character(),"mirna"=character(),"targets"=character(), "pvals.tfmirpair"=numeric(),
                          "pval.random"=numeric(),"zscore"=numeric(),"ENTREZ.IDs.coreg.targets"=character(),
                          "ENTREZ.IDs.motiftype.subnetwork"=character())
  index=c()
  for (i  in 1 : dim(tfmir.pairs)[1])
  {
    tf   = as.character( tfmir.pairs$tf[i])
    mirna= as.character( tfmir.pairs$mirna[i] )
    ## check for composite motifs condition
    if( dim(net[net$regulator == tf & net$target == mirna, ])[1] == 0   & dim(net[net$regulator == mirna & net$target == tf, ])[1] == 0 )
    {index=c(index,i)}
  }
  motifs=tfmir.pairs[index,]
  motifs
  
  if(dim(motifs)[1] > 0)
  {
    ######## Randomize the network and get significance of this type
    nodes.df=as.data.frame( unique(c(as.character(unlist(net$regulator)), as.character(unlist(net$target)))))  
    g=graph.data.frame(net[,c("regulator","target")],directed=T,vertices=nodes.df)
    tmpres=c()
    for( i in 1 : as.integer(config$NO_OF_Times_FOR_Network_Randomizations))
    {
      g.random=rewire(g, keeping_degseq(loops=TRUE, niter = 2*length(get.edgelist(g))))
      net.random=as.data.frame(get.edgelist(g.random,names=T))
      names(net.random)=c("regulator","target")
      tfmir.pairs.random=getPutativeTFmiRPairs(net.random)    
      #tfmir.pairs.random=getSignificantTFmiRpairs(tfmir.pairs.random,evidence)
      tfmir.pairs.random["pvals.tfmirpair"]=1
      tfmir.pairs.random=relaxMultipleTargetsForTFmiRPairs(tfmir.pairs.random)
      index=c()
      for ( j  in 1 : dim(tfmir.pairs.random)[1])
      {
        tf   = as.character( tfmir.pairs.random$tf[j])
        mirna= as.character( tfmir.pairs.random$mirna[j] )
        ## check for composite motifs condition
        if( dim(net.random[net.random$regulator == tf & net.random$target == mirna, ])[1] == 0   & dim(net.random[net.random$regulator == mirna & net.random$target == tf, ])[1] == 0 )
        {index=c(index,j)}
      }
      motifs.random=tfmir.pairs.random[index,]
      tmpres= c(tmpres,dim(motifs.random)[1])    
    }
    pval<-mean(tmpres > dim(motifs)[1])
    zscore=   ((dim(motifs)[1]) - (mean(tmpres))) / sd(tmpres)
    motifs["pval.random"]=pval
    motifs["zscore"]=zscore
    ##################################################################
    ###### get the entrez ids of co-regulated taget genes of each motif
    motifs["Symbol.coreg.targets"]=""
    motifs["ENTREZ.IDs.coreg.targets"]=""
    for( i in 1 : dim(motifs)[1])
    {
      tfmirpair=c(as.character(motifs[i,]$tf), as.character(motifs[i,]$mirna))
      targets=as.character( unlist(net[net$regulator %in% tfmirpair,]$target))
      targets.mir=targets[grep("hsa-",targets)]
      targets.gene=unique(c(as.character(motifs[i,]$tf), setdiff(targets,targets.mir)))
      motifs[i, ]$Symbol.coreg.targets=toString(targets.gene)
      motifs[i, ]$ENTREZ.IDs.coreg.targets=toString(as.vector(unlist(mget(as.character(targets.gene), envir=org.Hs.egALIAS2EG, ifnotfound=NA))))
    }
    
    ###### get the entrez ids of co-targetted genes by the same Tf and miRNA
    motifs["Symbol.cotargeted.genes"]=""
    motifs["ENTREZ.IDs.cotargeted.genes"]=""
    for( i in 1 : dim(motifs)[1])
    {
      targets.mirna= unlist(net[net$regulator == as.character(motifs[i,]$mirna) , ]$target)
      targets.mirna= setdiff(targets.mirna, targets.mirna [ grep("hsa-", targets.mirna)])  
      targets.tf= unlist(net[net$regulator == as.character(motifs[i,]$tf), ]$target)
      targets.tf= setdiff(targets.tf, targets.tf [ grep("hsa-", targets.tf)])
      common.targets=unique(intersect(targets.mirna,targets.tf))
      motifs[i, ]$Symbol.cotargeted.genes=toString(common.targets)
      motifs[i, ]$ENTREZ.IDs.cotargeted.genes=toString(as.vector(unlist(mget(as.character(common.targets), envir=org.Hs.egALIAS2EG, ifnotfound=NA))))
    }
    #   ###### get the entrez ids of subnetworks of each motif type
    #   tfs= as.character(unlist(motifs$tf))
    #   genes=as.character(unlist(motifs$targets))
    #   tfs.genes=unique(c(tfs,genes))
    #   tfs.genes.entrez=as.vector(unlist(mget(as.character(tfs.genes), envir=org.Hs.egALIAS2EG, ifnotfound=NA)))
    #   motifs["ENTREZ.IDs.motiftype.subnetwork"]=toString(tfs.genes.entrez)
    ######  add the type of the motif
    motif.type=replicate(dim(motifs)[1],"Coreg-FFL")
    motifs=cbind(motif.type,motifs)
    return (motifs)
  }
  
  return (motifs.empty)
}
