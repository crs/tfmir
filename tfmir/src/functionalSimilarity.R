
start <- function(genes, output) {


#setwd('/Users/chris/Website/TF-site/backend');

zz <- file(file.path(output, "funcSim.Rout"), open="wt")
sink(zz, type="message")

setwd(normalizePath(file.path(getwd(),'../../../backend')));
	#print (c('parameters are ', tf, mirna))
#warning(getwd());
print(getwd());
print('Before normalizing');
tf <- normalizePath(tf)
mirna <- normalizePath(mirna)
output <- normalizePath(output)
print('Source TF Mir');
print(source('TFMir.R'));
#setwd('../tf');
print('Call Functional Similarity');
PlotFunctionalSimilarity(tf, mirna,pval,evidence,disease,output);
print('Func sim finished');
}