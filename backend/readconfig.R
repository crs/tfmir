###########################################
#  - TFMir  Project                       #
#  - TFMir  read config function          #
#  - 2014-10-1                            #
#  - Copyright: Mohamed Hamed             #
###########################################

#R programming environments:
#- R studio IDE
#- R version 2.12.0 (2010-10-15)
#-Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
#-locale:   [1] C/en_US.UTF-8/C/C/C/C


## ==================================================
## read the configuration file and get the parameters 
## ==================================================

configfile="config.txt"
config=read.table(configfile)
key=apply(config, 1, function(x) strsplit(x,"=")[[1]][1])
value=apply(config, 1, function(x) strsplit(x,"=")[[1]][2])
config=as.list(value)
names(config)=key
config["transmir.file"]



## =====================================================
## create the log r file and its connection
## =====================================================
logfilename="/tmp/rlog.txt"
logfilename="rlog.txt"
append=T
#cat(paste("=========================  RUNNING FOR THE FIRST TIME AT ", Sys.time()  ," ====================",sep=""),file=logfilename,sep="\n")

writeToLog=function(msg)
{
  msg=paste("At ",Sys.time()," : ",msg)
  #cat(msg,file=logfilename,sep="\n",append=T)
}


