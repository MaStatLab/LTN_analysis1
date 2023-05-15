#!/usr/bin/env Rscript

library(phyloseq)
library(GetoptLong)
GetoptLong(
  "notu=i","number of OTUs.",
  "J=i","number of subjects.",
  "nj=i","max number of samples per subject.",
  "seed=i","random seed of simulation.",
  "WORK_DIR=s","working directory.",
  "scenario=s","H0, H1_single or H1_multi",
  "signal=f", "strength of signal. e.g., signal=2 means counts +200%",
  "ntop=i", "top x OTUs to choose from when injecting signals"
)

# set up cache directory
cache_dir = paste0(WORK_DIR, "/cache/cross_group_comparison/otu", notu, "/sim/")
system(paste0("mkdir -p ", cache_dir))

# output directory
output_dir = paste0(cache_dir, "/", scenario, "/J", J, "nj", nj)
if (scenario!="H0"){
  output_dir = paste0(output_dir, "/ntop", ntop, "signal", signal)
}
system(paste0("mkdir -p ", output_dir))
# output file name
output_file = paste0("sim_J", J, "nj", nj, "seed", seed, scenario, ".rds")

# check whether output file exists
if (!file.exists(paste0(output_dir, "/", output_file))){
  cat("generating data: ", paste0(output_dir, "/", output_file), "\n")
# load data
ps_otu = readRDS(paste0(WORK_DIR, "/cache", "/ps_otu_", notu,".rds"))
sampdat = data.frame(sample_data(ps_otu))
sampdat$Subject_ID = as.character(sampdat$Subject_ID)
subjects = unique(as.character(sampdat$Subject_ID))
# set random seed
set.seed(seed)
# H0 
if (scenario == "H0"){
  # draw J subjects from population
  subject_selected = sample(subjects, size = J, replace = F)
  # draw <= nj samples from each subject
  samples_selected = NULL
  for (j in 1:J){
    # obtain samples in subject j
    subject = subject_selected[j]
    sample_candidate = rownames(sampdat[sampdat$Subject_ID == subject,]) 
    # down-sample to nj samples if there are more than nj samples
    if (length(sample_candidate) > nj){
      samples = sample(sample_candidate, nj, replace = F)
    }
    # keep all samples if there are no more than nj samples
    else{
      samples = sample_candidate
    }
    samples_selected = c(samples_selected, samples)
  }
  
  # prune samples to obtain the subset of original data
  ps_subset = prune_samples(samples_selected, ps_otu)
  # obtain sample names
  samples = sample_names(ps_subset)
  # categorize half of the samples as "case"
  case_samples = sample(samples, ceiling(length(samples)/2), replace = F)
  sample_data(ps_subset)$group = as.numeric(rownames(sample_data(ps_subset)) %in% case_samples)
  # h0
  saveRDS(ps_subset, paste0(output_dir, "/", output_file))
}
if (scenario!="H0"){
  H0_file = gsub(scenario, "H0", paste0(dirname(output_dir), "/", output_file))
  ps_subset = readRDS(H0_file)
}
if (scenario=="H1_single"){
  ps_h1 = ps_subset
  OTU = otu_table(ps_subset)
  case_samples = sample_names(ps_h1)[which(sample_data(ps_h1)$group == 1)]
  # OTU candidates with high relative abundance
  relative_abundance = rowMeans(apply(OTU, 1, function(x){x/sum(x)}))
  otu_candidate = names(relative_abundance)[order(relative_abundance, decreasing = T)[1:ntop]]
  # sample 1 OTU from the top ntop
  otu_selected = sample(otu_candidate, 1)
  # increase the count of this single OTU in the case group: count +signal
  OTU[case_samples, otu_selected] = ceiling((1 + signal) * OTU[case_samples, otu_selected])
  otu_table(ps_h1) = OTU
  saveRDS(ps_h1, paste0(output_dir, "/", output_file))
}
if (scenario=="H1_multi"){
  ps_h1 = ps_subset
  OTU = otu_table(ps_subset)
  case_samples = sample_names(ps_h1)[which(sample_data(ps_h1)$group == 1)]
  # OTU candidates with high relative abundance
  relative_abundance = rowMeans(apply(OTU, 1, function(x){x/sum(x)}))
  otu_candidate = names(relative_abundance)[order(relative_abundance, decreasing = T)[1:ntop]]
  # sample 8 OTUs from the top ntop
  otu_selected = sample(otu_candidate, 8, replace = F)
  # count +signal, round to integer
  OTU[case_samples, otu_selected] = ceiling((1 + signal) * OTU[case_samples, otu_selected])
  otu_table(ps_h1) = OTU
  saveRDS(ps_h1, paste0(output_dir, "/", output_file))
}
}



