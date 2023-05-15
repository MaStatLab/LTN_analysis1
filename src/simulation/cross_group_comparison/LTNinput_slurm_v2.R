#!/usr/bin/env Rscript

library(GetoptLong)
library(phyloseq)
# input: sim_* files
# output: Xtest, Xadjust etc.
# arguments in fitting LTN: gibbs_crossgroup(N=N,p=p,g=g,r=r,YL=YL,Y=Y,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,niter=niter,adjust=adjust,reff=T,reffcov=reffcov,SEED=SEED,save_alpha_only=save_alpha_only,gprior_m=gm,pnull=pnull,a1a2='none',lambda_fixed=lambda,verbose=T)

# read arguments from command line 
GetoptLong(
  "WORK_DIR=s","working directory.",
  "config_file=s", "configuration file with input, output and resdir specification"
)
# source functions
source(paste0(WORK_DIR, '/src/LTN/utils.R'))
source(paste0(WORK_DIR, '/src/LTN/mixed_effects.R'))
id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
config = read.csv(config_file, header = T)
file = as.character(config[id, "raw"])
cat(file, "\n")
data = readRDS(file)
output_file = as.character(config[id, "input"])
if (!dir.exists(dirname(output_file))){
  system(paste0("mkdir -p ", dirname(output_file)))
}
if (!file.exists(output_file)){
  cat("output:", output_file, "\n")
  print(paste0("processing ", file, "\n"))
  # process data to the input arguments of gibbs_crossgroup()
  N = nsamples(data)
  p = ntaxa(data) - 1
  sampdat = sample_data(data)
  g = length(unique(sampdat$Subject_ID)) # g = J
  yyl = seqtab2y(otu_table(data), phy_tree(data))
  Y = yyl$Y
  YL = yyl$YL
  Xtest = matrix(sampdat$group, ncol = 1)
  # adjust for age, nationality and T1D status
  age=sampdat$Age_at_Collection
  age_s=(age-mean(age))/sd(age)
  # make model matrix for Xadjust
  covariates = data.frame(age_s = age_s, t1d = sampdat$Case_Control, country = sampdat$Country)
  Xadjust = model.matrix(~., covariates)
  rownames(Xadjust) = rownames(sampdat)
  # check for collinearity
  Xfull = cbind(Xtest, Xadjust)
  if (min(eigen(t(Xfull) %*% Xfull)$val) <= 0){
    sink(paste0(output_dir, gsub("rds", "log", output_file)))
    cat(id, "\n")
    cat("collinearity", "\n")
    sink()
  }
  # assign numeric subject labels ("grouplabel")
  grouplabel = as.numeric(as.factor(as.character(sampdat$Subject_ID))) # grouplabel is the numeric label of subject ID
  # check whether group label is 1:J
  if (!identical(sort(unique(grouplabel)), as.numeric(1:g))){
    stop("invalid group labels")
  }
  # save IMAGE as RData?
  cat("finish processing ", file, "\n")
  saveRDS(list(N = N, p = p, Y = Y, YL = YL, g = g, Xtest = Xtest, Xadjust = Xadjust, grouplabel = grouplabel), output_file)
}else{
  print(paste0("already processed ", file, "\n"))
}
