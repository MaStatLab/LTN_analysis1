#!/usr/bin/env Rscript

rm(list = ls())
library(GetoptLong)

ntop = NULL
signal = NULL

GetoptLong(
  "notu=i", "number of OTUs",
  "WORK_DIR=s", "working directory",
  "INPUT_DIR=s",".../LTNoutput",
  "J=i","number of subjects.",
  "nj=i","max number of samples per subject.",
  "signal=f", "strength of signal. e.g., signal=0.5 means counts +50%",
  "ntop=i", "top x OTUs to choose from when injecting signals (multiple OTUs)",
  "niter=i", "number of Gibbs iterations",
  "lambda=f", "lambda",
  "H=s", "H0, H1_single, or H1_multi",
  "nseed=i", "number of different random seeds"
)

input_dir = paste0(INPUT_DIR, "/", H, "/J", J, "nj", nj)
result_dir = paste0(WORK_DIR, "/results/simulation/cross_group_comparison/LTN/")
system(paste0("mkdir -p ", result_dir))
system(paste0("mkdir -p ", result_dir, "otu", notu))
if (H == "H0"){
  result_filename = paste0("otu",notu, "/lambda", lambda, "_niter", niter, "_J", J, "nj", nj, H, "_teststat.rds")
  files = paste0(input_dir, "/lambda", lambda, "_niter", niter, "_LTNoutput_sim_J", J, "nj", nj, "seed", 1:nseed, "H0.rds")
}else{
    files = paste0(input_dir, "/ntop", ntop, "signal", signal, "/lambda", lambda, "_niter", niter, "_LTNoutput_sim_J",J, "nj", nj, "seed", 1:nseed, H, ".rds")
    result_filename = paste0("otu",notu, "/lambda", lambda, "_niter", niter, "_J", J, "nj", nj, "signal", signal, "ntop", ntop, H, "_teststat.rds")
}
cat("save to: ", result_dir, "/", result_filename, '\n')
teststat = sapply(files, function(x){
  if (file.exists(x)){
    return(readRDS(x))
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
