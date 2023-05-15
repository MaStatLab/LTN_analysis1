#!/usr/bin/env Rscript

library(GetoptLong)
library(phyloseq)
rm(list = ls())
# get arguments
GetoptLong(
  "WORK_DIR=s","working directory.",
  "config_file=s", "configuration file with input, output and resdir specification",
  "lambda=f","lambda.",
  "niter=i","number of MCMC iterations"
)
# load functions
source(paste0(WORK_DIR, '/src/LTN/utils.R'))
source(paste0(WORK_DIR, '/src/LTN/mixed_effects.R'))
id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
config = read.csv(config_file, header = T)
file = as.character(config[id, "input"])
data = readRDS(file)
cat("read data from: ", file, "\n")
data = readRDS(file)
N = data$N
p = data$p
Y = data$Y
YL = data$YL
g = data$g
Xtest = data$Xtest
Xadjust = data$Xadjust
grouplabel = data$grouplabel
# specify parameters
r = 0
adjust = T
reffcov = 2
SEED = 1
save_alpha_only = T
gm = 100
pnull = 0.5

cachedir = as.character(config[id, "output_dir"])
gibbsdir=paste0(cachedir,'/gibbs/')
pmapdir=paste0(cachedir,'/pmap/')
pjapdir=paste0(cachedir,'/pjap/')
timedir=paste0(cachedir,'/time/')
sdir=paste0(cachedir,'/SEED/')
try(system(paste0('mkdir -p ',cachedir)))
try(system(paste0('mkdir -p ',gibbsdir)))
try(system(paste0('mkdir -p ',pmapdir)))
try(system(paste0('mkdir -p ',pjapdir)))
try(system(paste0('mkdir -p ',timedir)))
try(system(paste0('mkdir -p ',sdir)))
output_pjap_dir = as.character(config[id, "output_pjap_dir"])
try(system(paste0('mkdir -p ',output_pjap_dir)))
filenam = paste0("lambda", lambda, "_niter", niter, "_", as.character(config[id, "output_pjap_file"]))
output_file = paste0(output_pjap_dir, "/", filenam)
cat("pjap file: ", paste0(output_file), "\n")
if ((!file.exists(output_file))){
  # fit model
  st<-system.time(t<-try(gibbs<-gibbs_crossgroup(N=N,p=p,g=g,r=r,YL=YL,Y=Y,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,niter=niter,adjust=adjust,reff=T,reffcov=reffcov,SEED=SEED,save_alpha_only=save_alpha_only,gprior_m=gm,pnull=pnull,a1a2='none',lambda_fixed=lambda,verbose=T)))
  if ("try-error" %in% class(t)) {
    # SEED=SEED+1
    warning('numerical issues')
    saveRDS(t,paste0(sdir,'/',filenam))
    # t<-try(gibbs<-gibbs_crossgroup(N=N,p=p,g=g,r=r,YL=YL,Y=Y,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,niter=niter,adjust=adjust,reff=T,reffcov=reffcov,SEED=SEED,save_alpha_only=save_alpha_only,gprior_m=gm,pnull=pnull,a1a2='none',lambda_fixed=lambda))
  }
  BETA1=lapply(gibbs$BETA1,function(x){as.vector(x)})
  BETAMAT1=do.call(rbind,BETA1)
  PMAP=matrix(apply(BETAMAT1[(niter/2+1):niter,],2,function(x){sum(x!=0)})/(niter/2),nrow = 1)
  saveRDS(PMAP,paste0(pmapdir,'/',filenam))
  all0=rowSums(BETAMAT1[(niter/2+1):niter,]!=0)
  PJAP=1-sum(all0==0)/(niter/2)
  saveRDS(PJAP,paste0(pjapdir,'/',filenam))
  saveRDS(PJAP, output_file)
  saveRDS(st,paste0(timedir,'/',filenam))
  print(st)
  saveRDS(BETAMAT1, paste0(gibbsdir,"/",filenam))
}else{
  cat("pjap file already exists: ", paste0(pjapdir, filenam), "\n")
}


