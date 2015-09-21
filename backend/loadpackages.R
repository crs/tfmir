###########################################
#  - TFMir  Project                       #
#  - TFMir  init function                 #
#  - 2014-10-1                            #
#  - Copyright: Mohamed Hamed             #
###########################################

#R programming environments:
#- R studio IDE
#- R version 2.12.0 (2010-10-15)
#-Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
#-locale:   [1] C/en_US.UTF-8/C/C/C/C

`%ni%`<- Negate(`%in%`)

## ==========================================
## Check for R packages if they are installed  
## ==========================================


pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    #install.packages("x", repos = "http://cran.us.r-project.org")
    source("http://bioconductor.org/biocLite.R")
    biocLite(x)
    #library("x")
    #require(x)
    sapply(x, require, character.only = TRUE)
  }

}

if (!require('Cairo')) {
	install.packages('Cairo')
}



## =====================================================
## check for bioconductor packages if they are installed
## =====================================================

pkgTest("gdata")
pkgTest("multtest")
pkgTest("org.Hs.eg.db")
pkgTest("igraph")
pkgTest("GOSemSim")

#library(c('gdata', 'multtest','org.Hs.eg.db','igraph'))




# 
# 
# pkgs <- library()$results
# pkgs <- pkgs[,1]
# if (length(which(pkgs == 'ape'))<1) {
#   stop('This package depedends on the package \'ape\',please install \'lars\' first, using command \'install.packages(\'ape\')\'.')
# } else if (length(which(pkgs == 'KEGG.db'))<1) {
#   stop('This package depedends on the package \'KEGG.db\',please install \'lars\' first, using command \'install.packages(\'KEGG.db\')\'.')
# } else if (length(which(pkgs == 'HOBBA'))<1) {
#   stop('This package depedends on the package \'HOBBA\',please install \'HOBBA\' first, using command \'install.packages(\'HOBBA\')\'.')
# }
