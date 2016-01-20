###########################################
#  - TFMir  Project                       #
#  - TFMir  statistics functions          #
#  - 2014-10-1                            #
#  - Copyright: Mohamed Hamed             #
###########################################

#R programming environments:
#- R studio IDE
#- R version 2.12.0 (2010-10-15)
#-Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
#-locale:   [1] C/en_US.UTF-8/C/C/C/C


## ===============================
## significance_overlapbysimulation
## ===============================
countOverlap<-function(total, numgA, numgB, replace=FALSE){
  groupA<-sample.int(total, numgA, replace=replace)
  groupB<-sample.int(total, numgB, replace=replace)
  return(length(intersect(groupA, groupB)))
}

overlapSignificance_By_simulation=function(total, numgA, numgB,overlap,no_of_simulation, replace=FALSE )
{
  tmpres<-replicate(no_of_simulation, countOverlap(total, numgA, numgB))
  pval<-mean(tmpres >= overlap)
  pval
}

overlapSignificance_By_HyperGEOM=function(total, numgA, numgB,overlap, replace=FALSE )
{
#   If you want to test a number overlap, then the probability of getting that number or smaller from this model is
#   phyper(overlap, sampleb, totala - sampleb, samplec) 
#   and of getting that number or larger is
#   1 - phyper(overlap - 1, sampleb, totala - sampleb, samplec)
#   or   phyper(overlap - 1, sampleb, totala - sampleb, samplec,lower.tail=FALSE)

  pval=phyper(overlap-1, numgA, total - numgA, numgB,lower.tail=FALSE) 
  pval
}

overlapSignificance_By_Python=function(total, numgA, numgB,overlap, replace=FALSE )
{
  stat=paste( "python list_overlap_p.py ",overlap," ",total," ",numgA," ",numgB,sep="")
  pval=system(stat)
  pval
}


DO_ORA_FOR_MIRNA=function(mirnas,category,pval.cutoff)
{
  association.file=""
  if(category =="function") { association.file=mirna.function.association
  } else if(category =="disease") association.file=mirna.disease.association
  terms=unique(unlist(association.file[association.file[,1] %in% mirnas,2]))
  total= length(unique(association.file$mirna))
  terms.df=data.frame("Category"=character(),"Term"=character(),"Count"=numeric(),"Mir"=character(),"Percentage"=numeric(),"Pval"=numeric())
  
  for(i in 1: length(terms))
  {
    t=terms[i]
    mirna.t=unique(unlist(association.file[association.file[,2] %in% t,1]))
    overlap=intersect(mirna.t,mirnas)
    overlap.count=length(overlap)
    numgA=length(mirnas)
    numgB=length(mirna.t)
    pval=phyper(overlap.count-1, numgA, total - numgA, numgB,lower.tail=FALSE)
    perc= (overlap.count / numgB)* 100
    #print(c("Category"=category,"Term"=t,"Count"=overlap.count,"Mir"=toString(overlap),"Percentage"=perc,"Pval"=pval))
    if ((overlap.count > 0) & !is.na(pval)) {
      row=data.frame("Category"=category,"Term"=t,"Count"=overlap.count,"Mir"=toString(overlap),"Percentage"=perc,"Pval"=pval)
      terms.df=rbind(terms.df,row)
    }
  }

  pval.BH=mt.rawp2adjp(terms.df$Pval,proc="BH")
  terms.df=terms.df[pval.BH$index,]
  terms.df["Pval.BH"]=(pval.BH$adjp)[,2]
  
  pval.Bonf=mt.rawp2adjp(terms.df$Pval,proc="Bonferroni")
  terms.df=terms.df[pval.Bonf$index,]
  terms.df["Pval.Bonf"]=(pval.Bonf$adjp)[,2]
  
  terms.df=terms.df[terms.df$Pval <= pval.cutoff,]
}




getMiRNAsRegulatorsofTFs=function( tfs, pval.cutoff,evidence)
{
  ###  mirna-gene interactions
  db=dbs[dbs$category=="mirna-gene" & dbs$evidence %in% evidence, ]
  tfs.regulators= unique(as.character(unlist(db[db$target %in%  tfs, ]$regulator)))
  #total= length(unique(as.character(unlist(db$target))))
  total= as.integer(config$Total.No.of.genes.in.human)
  pvals=c()
  for(i in 1: length(tfs.regulators))
  {
    mirna=tfs.regulators[i]
    mirna.targets=unique(as.character(unlist(db[db$regulator %in% mirna,]$target)))
    overlap=length(intersect(tolower(mirna.targets),tolower(tfs)))
    numgA=length(tfs)
    numgB=length(mirna.targets)
    #total=as.integer(config$Total.No.of.genes.in.human)
    pvals=c(pvals,phyper(overlap-1, numgA, total - numgA, numgB,lower.tail=FALSE) )    
  }
  pvals.BH=mt.rawp2adjp(pvals,proc="BH")
  tfs.regulators=tfs.regulators[pvals.BH$index]
  pvals.BH=(pvals.BH$adjp)[,2]
  return(tfs.regulators[which(pvals.BH < pval.cutoff)])
  
}

getMiRNAsTargetsofTFs=function( tfs, pval.cutoff,evidence)
{
  ###  TF-miRNA interactions
  db=dbs[dbs$category=="tf-mirna" & dbs$evidence %in% evidence, ]
  tfs.targets= unique(as.character(unlist(db[db$regulator %in%  tfs, ]$target)))
  #total= length(unique(as.character(unlist(db$regulator))))
  total= as.integer(config$Total.No.of.genes.in.human)
  
  pvals=c()
  for(i in 1: length(tfs.targets))
  {
    mirna=tfs.targets[i]
    mirna.regulators=unique(as.character(unlist(db[db$target %in% mirna,]$regulator)))
    overlap=length(intersect(tolower(mirna.regulators),tolower(tfs)))
    numgA=length(tfs)
    numgB=length(mirna.regulators)
    #total=as.integer(config$Total.No.of.genes.in.human)
    pvals=c(pvals,phyper(overlap-1, numgA, total - numgA, numgB,lower.tail=FALSE) )    
  }
  pvals.BH=mt.rawp2adjp(pvals,proc="BH")
  tfs.targets=tfs.targets[pvals.BH$index]
  pvals.BH=(pvals.BH$adjp)[,2]
  return(tfs.targets[which(pvals.BH<pval.cutoff)])  
}




##### for the third scenarion : when user inputs a list of miRNAs only
### then u get the list of Tfs(regulators) and genes(targets) which are 
## statistically eniched in the input miRNA list


getTFsRegulatorsofMiRNAs=function( mirnas, pval.cutoff,evidence)
{
  ###  mirna-gene interactions
  db=dbs[dbs$category=="tf-mirna" & dbs$evidence %in% evidence, ]
  tfs.regulators= unique(as.character(unlist(db[db$target %in%  mirnas, ]$regulator)))
  #total= length(unique(as.character(unlist(db$target))))
  total= as.integer(config$Total.No.of.miRNA.in.human)
  pvals=c()
  for(i in 1: length(tfs.regulators))
  {
    tf=tfs.regulators[i]
    mirna.targets=unique(as.character(unlist(db[db$regulator %in% tf,]$target)))
    overlap=length(intersect(tolower(mirna.targets),tolower(mirnas)))
    numgA=length(mirnas)
    numgB=length(mirna.targets)
    pvals=c(pvals,phyper(overlap-1, numgA, total - numgA, numgB,lower.tail=FALSE) )    
  }
  pvals.BH=mt.rawp2adjp(pvals,proc="BH")
  tfs.regulators=tfs.regulators[pvals.BH$index]
  pvals.BH=(pvals.BH$adjp)[,2]
  return(tfs.regulators[which(pvals.BH < pval.cutoff)])
  
}

getTFsTargetsofMiRNAs=function( mirnas, pval.cutoff,evidence)
{
  ###  TF-miRNA interactions
  db=dbs[dbs$category=="mirna-gene" & dbs$evidence %in% evidence, ]
  tfs.targets= unique(as.character(unlist(db[db$regulator %in%  mirnas, ]$target)))
  #total= length(unique(as.character(unlist(db$regulator))))
  total= as.integer(config$Total.No.of.miRNA.in.human)
  pvals=c()
  for(i in 1: length(tfs.targets))
  {
    tf=tfs.targets[i]
    mirna.regulators=unique(as.character(unlist(db[db$target %in% tf,]$regulator)))
    overlap=length(intersect(tolower(mirna.regulators),tolower(mirnas)))
    numgA=length(mirnas)
    numgB=length(mirna.regulators)
    pvals=c(pvals,phyper(overlap-1, numgA, total - numgA, numgB,lower.tail=FALSE) )    
  }
  pvals.BH=mt.rawp2adjp(pvals,proc="BH")
  tfs.targets=tfs.targets[pvals.BH$index]
  pvals.BH=(pvals.BH$adjp)[,2]
  return(tfs.targets[which(pvals.BH<pval.cutoff)])  
}



