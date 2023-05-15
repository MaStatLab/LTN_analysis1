#!/usr/bin/env Rscript

library(GetoptLong)
library(phyloseq)
rm(list = ls())
# get arguments
GetoptLong(
  "WORK_DIR=s","working directory.",
  "config_file=s", "configuration file with input, output and resdir specification",
  "r=i","number of factors.",
  "niter=i","number of MCMC iterations"
)
id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
config = read.csv(config_file, header = T)
input_file = as.character(config[id, "input"])
resdir = as.character(config[id, "output"])
dirnam = paste0("r",r,"niter",niter)
# source files
workdir = paste0(WORK_DIR,'/src/simulation/cross_group_comparison/dirfactor/')
source(paste0(workdir, "/R/MCMC.R"))
source(paste0(workdir, "/R/utilities.R"))
source(paste0(workdir, '/R/sim_compare_functions.R'))

# cache directory
system(paste0('mkdir -p ',resdir,'/gibbs'))
system(paste0('mkdir -p ',resdir,'/teststat'))

# load data 
cat("input:", input_file, "\n")
data = readRDS(input_file)
# check and fit
filenam = paste0(dirnam, "_v_post.rds")
if (!file.exists(paste0(resdir,'/teststat/', filenam))){
  cat(id, "fit and save to: ", resdir,'/teststat/', filenam, "\n")
  cat("id: ", id, ", fit and save to: ", paste0(resdir,'/teststat/', filenam), "\n")
  # read data
  cnt = otu_table(data)
  data0 = t(cnt)
  n = ncol(data0)
  p = nrow(data0)
  # blocking factors (random effects)
  sampdat = sample_data(data)
  grouplabel=data.frame(id=as.character(sampdat$Subject_ID))
  sub.design = t(model.matrix( ~.-1, grouplabel ))
  n.sub = nrow(sub.design)
  # hyper: stays the same for all simulations
  # m is number of factors
  # all hyperparameters other than m=r are default
  hyper = list( nv = 30, a.er = 1, b.er = 0.3, a1 = 10, a2 = 20, m = r,
                sub.design = sub.design, alpha = 10, beta = 0, a.x.sigma=5, b.x.sigma=5 )
  sigma.value = seq(0.001,0.999,0.001)
  tmp = c(0,pbeta( sigma.value, 10/68, 1/2-10/68 ))
  sigma.prior = sapply( 1:999, function(x) tmp[x+1]-tmp[x] )
  sigma.prior[length(sigma.prior)] = sigma.prior[length(sigma.prior)] + 1-sum(sigma.prior)
  # fixed effects
  Xtest = as.numeric(sampdat$group)
  # age: standardize
  age=sampdat$Age_at_Collection
  age_s = (age-mean(age))/sd(age)
  # t1d and country
  t1d = as.factor(sampdat$Case_Control)
  country = as.factor(sampdat$Country)  
  fixed_effects = data.frame(Xtest,age_s,t1d,country)
  fixed_effects_design = model.matrix( ~., fixed_effects)
  # make Xtest the first column
  fixed_effects_design = cbind(fixed_effects_design[,-1], fixed_effects_design[,1]) 
  colnames(fixed_effects_design)[5] = "intercept"
  y.fix.inter=t(fixed_effects_design)
  # start
  start.sub = list( er = 1/rgamma( 1, shape = hyper$a.er, rate = hyper$b.er ),
                    sigma = sample( sigma.value, size = p, replace = T, prob = sigma.prior ),
                    T.aug = rgamma( n, 10, 1 ),
                    Q = matrix( 0.5, nrow = p, ncol = n ),
                    X = matrix( rnorm( hyper$m*p ), nrow = hyper$m ),
                    Y.sub = matrix( 0, nrow = hyper$m, ncol = n.sub ),
                    delta = c( rgamma( 1, shape = hyper$a1, rate = 1 ), rgamma( hyper$m-1, shape = hyper$a2, rate = 1 ) ),
                    phi = matrix( rgamma( hyper$m*n.sub, shape = 3/2, rate = 3/2 ), nrow = n.sub ),
                    y = y.fix.inter,
                    x = matrix( rnorm( p*nrow(y.fix.inter) ), nrow=nrow(y.fix.inter) ),
                    x.sigma = rgamma( nrow(y.fix.inter), 1, 1 ) )
  # make output directory
  if (!dir.exists(paste0(resdir, '/gibbs/', dirnam))){
    dir.create(paste0(resdir, '/gibbs/', dirnam), recursive = T)
  }
  nod.free.mcmc.new( data0, start.sub, hyper, sigma.value, sigma.prior,paste0(resdir, '/gibbs/', dirnam), burnin = 0.5, step = niter, thin = 10 )
  # calculate teststat
  sampdir=paste0(resdir, '/gibbs/', dirnam)
  V=sapply(list.files(sampdir),function(filename){readRDS(paste0(sampdir,"/",filename))$x[1,]})
  v_post=apply(V,1,mean)
  v_post_norm=sqrt(sum(v_post^2))
  saveRDS(v_post, paste0(resdir, '/teststat/', filenam))
  saveRDS(v_post_norm, paste0(resdir, '/teststat/', gsub("v_post", "v_post_norm", filenam)))
}else{
  cat("id: ", id, "results already exist", "\n")
}
