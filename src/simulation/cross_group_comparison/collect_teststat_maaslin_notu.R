#!/usr/bin/env Rscript

rm(list = ls())
library(GetoptLong)

ntop = NULL
signal = NULL

GetoptLong(
  "notu=i", "number of OTUs",
  "WORK_DIR=s", "working directory",
  "INPUT_DIR=s",".../MaAsLin2",
  "J=i","number of subjects.",
  "nj=i","max number of samples per subject.",
  "signal=f", "strength of signal. e.g., signal=0.5 means counts +50%",
  "ntop=i", "top x OTUs to choose from when injecting signals",
  "H=s", "H0, H1_single, or H1_multi",
  "nseed=i", "number of different random seeds"
)

input_dir = paste0(INPUT_DIR, "/", H, "/J", J, "nj", nj)
result_dir = paste0(WORK_DIR, "/results/simulation/cross_group_comparison/MaAsLin2/", "otu",notu,"/")
system(paste0("mkdir -p ", result_dir))
if (H == "H0"){
  result_filename = paste0("J", J, "nj", nj, H, "_teststat.rds")
  files = paste0(input_dir, "/sim_J", J, "nj", nj, "seed", 1:nseed, H, "/all_results.tsv")
}else{
    files = paste0(input_dir, "/ntop", ntop, "signal", signal, "/sim_J",J, "nj", nj, "seed", 1:nseed, H, "/all_results.tsv")
    result_filename = paste0("J", J, "nj", nj, "signal", signal, "ntop", ntop, H, "_teststat.rds")
}
#if (!file.exists(paste0(result_dir,'/',result_filename))){
cat("save to: ", result_dir, "/", result_filename, '\n')
teststat = - sapply(files, function(x){
  if (file.exists(x)){
    res = read.csv(x,sep='\t')
    return(min(p.adjust(p = res[res$metadata == "Xtest", "pval"], method = "BH")))
  }else{
    return(NA)
  }
})
teststat = na.omit(teststat)
cat(length(teststat), '\n')
if (length(teststat)!=nseed){
  warning("less files than nseed", "\n")
}
if (length(teststat) > 0){
  saveRDS(teststat, paste0(result_dir,'/',result_filename))
}else{
  cat("no data available", "\n")
}
#}else{
#  cat(result_dir,'/',result_filename, "exists", "\n")
#}
