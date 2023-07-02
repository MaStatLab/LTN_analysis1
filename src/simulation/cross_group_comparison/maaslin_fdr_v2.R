#!/usr/bin/env Rscript

rm(list = ls())
library(GetoptLong)
library(phyloseq)

ntop = NULL
signal = NULL
fdr_target = 0.05

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
  result_filename = paste0("J", J, "nj", nj, H, "_fdr005.rds")
  files = paste0(input_dir, "/sim_J", J, "nj", nj, "seed", 1:nseed, H, "/all_results.tsv")
}else{
  files = paste0(input_dir, "/ntop", ntop, "signal", signal, "/sim_J",J, "nj", nj, "seed", 1:nseed, H, "/all_results.tsv")
  result_filename = paste0("J", J, "nj", nj, "signal", signal, "ntop", ntop, H, "_fdr005.rds")
}

fdr_i = function(i){
# load simulated dataset
if (H != 'H0'){
  H1 = readRDS(paste0(WORK_DIR, "/cache/cross_group_comparison/otu", notu, "/sim/", H, "/J", J, "nj", nj, "/ntop", ntop, "signal", signal, "/sim_J", J, "nj", nj, "seed",i, H, ".rds"))
}
H0 = readRDS(paste0(WORK_DIR, "/cache/cross_group_comparison/otu", notu, "/sim/", "H0", "/J", J, "nj", nj, "/sim_", "J", J, "nj", nj, "seed",i, "H0", ".rds"))

# obtain truth d
if (H == "H0"){
  sim = H0
}else{
  sim = H1
}
d = data.frame(d = colSums(otu_table(sim)) != colSums(otu_table(H0)))
# obtain d hat
res = read.csv(files[i],sep='\t')
res = res[res$metadata == "Xtest", c("feature", "pval")]
res$d_hat = p.adjust(p = res$pval, method = 'BH') <= fdr_target
# res$d_hat = res$qval <= fdr_target
res$d = d[gsub('X', '', res$feature),]
if (sum(res$d_hat) > 0){
  fdr = sum(res$d_hat*(1-res$d))/sum(res$d_hat)
} else{
  fdr = 0
}
return(fdr)
}

fdr_all = sapply(1:500, fdr_i)
saveRDS(fdr_all, paste0(result_dir,'/',result_filename))

