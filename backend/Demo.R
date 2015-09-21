
###########################################
#  - TFMir  Project                       #
#  - TFMir  Demo test for the TFMIR script#
#  - 2014-10-1                            #
#  - Copyright: Mohamed Hamed             #
###########################################

#R programming environments:
#- R studio IDE
#- R version 2.12.0 (2010-10-15)
#-Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
#-locale:   [1] C/en_US.UTF-8/C/C/C/C


tf.path="nksea.tfmir.txt"
mirna.path="mirna.sample.txt"

tf.path="genes.bc.diff.txt"
mirna.path="mirna.bc.diff.txt"

tf.path="genes.ad.diff.txt"
mirna.path="mirna.ad.diff.txt"



mirna.path=""
pval.cutoff=0.05
disease="Melanoma"
disease="Neoplasm"
disease="halabessahbelzabadymolokhia"
disease="Breast Neoplasms"
disease="alzheimer"
disease=""

output.path="output"
evidence="Experimental"
evidence="Predicted"
evidence="both"

source("TFMIR.R")
TFMir(tf.path,mirna.path,pval.cutoff,evidence,disease,output.path)

