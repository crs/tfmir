openOpt=function(input_file, solver){
  stat=paste("python dsp.py",input_file,solver, sep = " ")
  pval=system(stat)
  result_file = paste(gsub(".txt", "", input_file), "dsp", solver, sep = "_")
  result_file = paste(result_file, "txt", sep = ".")
  result = res.table(result_file)
  return(as.character(result$V1))
}