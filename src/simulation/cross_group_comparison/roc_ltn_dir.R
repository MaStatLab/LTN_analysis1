#!/usr/bin/env Rscript

# library(GetoptLong)
library(ROCR)
# read arguments from command line 
# GetoptLong(
#   "WORK_DIR=s","working directory.",
#   "J=i","number of subjects.",
#   "nj=i","max number of samples per subject.",
#   "signal_single=f", "strength of signal_single signal. e.g., signal_single=2 means counts +200%",
#   "signal_multi=f", "strength of signal_multi signal. e.g., signal_multi=0.5 means counts +50%",
#   "ntop_single=i", "top x OTUs to choose from when injecting signals (single OTU)",
#   "ntop_multi=i", "top x OTUs to choose from when injecting signals (multiple OTUs)"
# )
WORK_DIR = "REPLACE_WITH_WORKDIR/"
niter = 10000
nseed = 500
J = 33
nj = 5

par(mfrow = c(1,2))
# single, lambda = 10, h0
lambda = 10
output_dir = paste0(WORK_DIR, "/cache/cross_group_comparison/sim_new/LTNoutput/pjap/")
files0 = sapply(1:nseed, function(seed)paste0("lambda", lambda, "_niter", niter, "_LTNoutput_ps_subset_J", J, "nj", nj, "seed", seed, "H0.rds"))
pjap0 = sapply(paste0(output_dir, files0), readRDS)
print(length(pjap0))
saveRDS(pjap0, "REPLACE_WITH_WORKDIR/results/simulation/cross_group_comparison/lambda10_niter10000_J33nj5ntop_multi20ntop_single20signal_multi0.5signal_single2_pjap_H0.rds")
# H1
files = sapply(1:nseed, function(seed)paste0("lambda", lambda, "_niter", niter, "_LTNoutput_sim_J", J, "nj", nj, "seed", seed, "H1_single.rds"))
pjap1 = sapply(paste0(output_dir, files), readRDS)
print(length(pjap1))
saveRDS(pjap1, "REPLACE_WITH_WORKDIR/results/simulation/cross_group_comparison/lambda10_niter10000_J33nj5ntop_multi20ntop_single20signal_multi0.5signal_single2_pjap_single.rds")
labels = c(rep(0, length(pjap0)), rep(1, length(pjap1)))
pred = c(pjap0, pjap1)
ltn_single = performance(prediction(pred, labels), "tpr", "fpr")
plot(ltn_single)

# single, lambda = 1
# H0
lambda = 1
output_dir = paste0(WORK_DIR, "/cache/cross_group_comparison/sim_new/LTNoutput/pjap/")
files0 = sapply(1:nseed, function(seed)paste0("lambda", lambda, "_niter", niter, "_LTNoutput_ps_subset_J", J, "nj", nj, "seed", seed, "H0.rds"))
pjap0 = sapply(paste0(output_dir, files0), function(x){
  if(file.exists(x)){
    return(readRDS(x))
  }else{
    return(NA)
  }
})
pjap0=na.omit(pjap0)
print(length(pjap0))
saveRDS(pjap0, "REPLACE_WITH_WORKDIR/results/simulation/cross_group_comparison/lambda1_niter10000_J33nj5ntop_multi20ntop_single20signal_multi0.5signal_single2_pjap_H0.rds")

files = sapply(1:nseed, function(seed)paste0("lambda", lambda, "_niter", niter, "_LTNoutput_sim_J", J, "nj", nj, "seed", seed, "H1_single.rds"))
pjap1 = sapply(paste0(output_dir, files), function(x){
  if(file.exists(x)){
    return(readRDS(x))
  }else{
    return(NA)
  }
})
pjap1=na.omit(pjap1)
print(length(pjap1))
saveRDS(pjap1, "REPLACE_WITH_WORKDIR/results/simulation/cross_group_comparison/lambda1_niter10000_J33nj5ntop_multi20ntop_single20signal_multi0.5signal_single2_pjap_single.rds")

labels = c(rep(0, length(pjap0)), rep(1, length(pjap1)))
pred = c(pjap0, pjap1)
ltn_single = performance(prediction(pred, labels), "tpr", "fpr")
plot(ltn_single, add = T, col = "blue")

legend("bottomright", legend = c("lambda = 10", "lambda = 1"), fill = c("black", "blue"))

# multi, lambda = 10
# H0
lambda = 10
output_dir = paste0(WORK_DIR, "/cache/cross_group_comparison/sim_new/LTNoutput/pjap/")
files0 = sapply(1:nseed, function(seed)paste0("lambda", lambda, "_niter", niter, "_LTNoutput_ps_subset_J", J, "nj", nj, "seed", seed, "H0.rds"))
pjap0 = sapply(paste0(output_dir, files0), function(x){
  if(file.exists(x)){
    return(readRDS(x))
  }else{
    return(NA)
  }
})
pjap0=na.omit(pjap0)
print(length(pjap0))

files = sapply(1:nseed, function(seed)paste0("lambda", lambda, "_niter", niter, "_LTNoutput_sim_J", J, "nj", nj, "seed", seed, "H1_multi.rds"))
pjap1 = sapply(paste0(output_dir, files), function(x){
  if(file.exists(x)){
    return(readRDS(x))
  }else{
    return(NA)
  }
})
pjap1=na.omit(pjap1)
print(length(pjap1))
labels = c(rep(0, length(pjap0)), rep(1, length(pjap1)))
pred = c(pjap0, pjap1)
ltn_single = performance(prediction(pred, labels), "tpr", "fpr")
plot(ltn_single)

# multi, lambda = 1
# h0
lambda = 1
output_dir = paste0(WORK_DIR, "/cache/cross_group_comparison/sim_new/LTNoutput/pjap/")
files0 = sapply(1:nseed, function(seed)paste0("lambda", lambda, "_niter", niter, "_LTNoutput_ps_subset_J", J, "nj", nj, "seed", seed, "H0.rds"))
pjap0 = sapply(paste0(output_dir, files0), function(x){
  if(file.exists(x)){
    return(readRDS(x))
  }else{
    return(NA)
  }
})
pjap0=na.omit(pjap0)
print(length(pjap0))
# H1
files = sapply(1:nseed, function(seed)paste0("lambda", lambda, "_niter", niter, "_LTNoutput_sim_J", J, "nj", nj, "seed", seed, "H1_multi.rds"))
pjap1 = sapply(paste0(output_dir, files), function(x){
  if (file.exists(x)){
    return(readRDS(x))
  }else{
    return(NA)
  }
})
pjap1 = na.omit(pjap1)
print(length(pjap1))
labels = c(rep(0, length(pjap0)), rep(1, length(pjap1)))
pred = c(pjap0, pjap1)
ltn_single = performance(prediction(pred, labels), "tpr", "fpr")
plot(ltn_single, add = T, col = "blue")
legend("bottomright", legend = c("lambda = 10", "lambda = 1"), fill = c("black", "blue"))


#single, top 5
lambda = 10
ntop_multi = 20
ntop_single = 5
signal_multi = 0.5
signal_single = 2
dirs = sapply(1:nseed, function(x)paste0(WORK_DIR, "/cache/cross_group_comparison/sim_new/", "J", J, "nj", nj, "ntop_multi", ntop_multi, "ntop_single", ntop_single, "seed", x, "signal_multi", signal_multi, "signal_single", signal_single, "/LTNoutput/pjap"))
sum(dir.exists(dirs))
files = sapply(1:nseed, function(seed)paste0(dirs[seed], "/lambda", lambda, "_niter", niter, "_LTNoutput_sim_J", J, "nj", nj, "seed", seed, "H1_", "single", ".rds"))
pjap_single = sapply(files, function(x){
  if(file.exists(x)){
    return(readRDS(x))
  }else{
    return(NA)
  }
})
pjap_single = na.omit(pjap_single)
print(length(pjap_single))
saveRDS(pjap_single, "REPLACE_WITH_WORKDIR/results/simulation/cross_group_comparison/lambda10_niter10000_J33nj5ntop_multi20ntop_single5signal_multi0.5signal_single2_pjap_single.rds")
labels = c(rep(0, length(pjap0)), rep(1, length(pjap_single)))
pred = c(pjap0, pjap_single)
ltn_single = performance(prediction(pred, labels), "tpr", "fpr")
plot(ltn_single)

# single, +400%
ntop_multi = 20
ntop_single = 20
signal_multi = 0.5
signal_single = 4
dirs = sapply(1:nseed, function(x)paste0(WORK_DIR, "/cache/cross_group_comparison/sim_new/", "J", J, "nj", nj, "ntop_multi", ntop_multi, "ntop_single", ntop_single, "seed", x, "signal_multi", signal_multi, "signal_single", signal_single, "/LTNoutput/pjap"))
sum(dir.exists(dirs))
files = sapply(1:nseed, function(seed)paste0(dirs[seed], "/lambda", lambda, "_niter", niter, "_LTNoutput_sim_J", J, "nj", nj, "seed", seed, "H1_", "single", ".rds"))
pjap_single = sapply(files, function(x){
  if(file.exists(x)){
    return(readRDS(x))
  }else{
    return(NA)
  }
})
pjap_single = na.omit(pjap_single)
print(length(pjap_single))
labels = c(rep(0, length(pjap0)), rep(1, length(pjap_single)))
pred = c(pjap0, pjap_single)
ltn_single = performance(prediction(pred, labels), "tpr", "fpr")
plot(ltn_single, add = T, col = 'blue')

