#!/usr/bin/env Rscript
rm(list = ls())

library(Maaslin2)
library(phyloseq)

# set maaslin work dir
WORK_DIR = "REPLACE_WITH_WORKDIR/"

# load config file
config = read.csv("REPLACE_WITH_WORKDIR//cache/cross_group_comparison/otu100/maaslin_confignj10_multi2.csv")


# subset files
# files = sort(grep("H", list.files(paste0(WORK_DIR, "/cache/cross_group_comparison/sim_new/"), full.names = T), val = T))

# array id equivalence
for (id in 1:nrow(config)){
  cat("maaslin id: ", id, "\n")
  # read file as in ltn
  file = as.character(config[id, "input"])
  data = readRDS(file)
  # output file name
  # dirnam = strsplit(rev(strsplit(file,'/')[[1]])[1], "\\.")[[1]][1]
  dirnam = as.character(config[id, "output"])
  # check if output file already exists 
  # if (!file.exists(paste0(WORK_DIR, "/cache/cross_group_comparison/MaAsLin2/", dirnam, "/significant_results.tsv"))){
  output_file = paste0(dirnam, "/significant_results.tsv")
  if (!file.exists(output_file)){
    cat("fit and save to: ", output_file, "\n")
    # preprocess data: explicitly set datatype as either factors or numeric
    input_data = data.frame(otu_table(data))
    sampdat = sample_data(data)
    input_metadata = data.frame(sample = rownames(sampdat))
    rownames(input_metadata) = rownames(sampdat)
    # group label: convert to factor
    input_metadata$Xtest = as.factor(sampdat$group)
    # age: standardize
    age=sampdat$Age_at_Collection
    input_metadata$age_s = (age-mean(age))/sd(age)
    # t1d and country
    input_metadata$t1d = as.factor(sampdat$Case_Control)
    input_metadata$country = as.factor(sampdat$Country)
    # random effect: Subject_ID
    input_metadata$Subject_ID = as.factor(as.character(sampdat$Subject_ID))
    # fit model
    # output_dir = paste0(WORK_DIR, "/cache/cross_group_comparison/MaAsLin2/", dirnam)
    output_dir = dirnam
    system(paste0("mkdir -p ", output_dir))
    fit1 = Maaslin2(input_data = input_data,
                    input_metadata = input_metadata, 
                    output = output_dir, 
                    fixed_effects = c("Xtest", "age_s", "t1d", "country"),
                    random_effects = "Subject_ID",
                    plot_heatmap = F,
                    plot_scatter = F)
  }else{
    cat("file already exists: ", output_file, "\n")
  }
}
