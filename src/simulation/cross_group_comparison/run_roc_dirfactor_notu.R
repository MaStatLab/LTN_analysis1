#!/usr/bin/env Rscript
rm(list=ls())
library(ROCR)
notu = 50
nseed = 500
roc_script = "REPLACE_WITH_WORKDIR/src/simulation/cross_group_comparison/roc_simv1.R"
result_dir = paste0("REPLACE_WITH_WORKDIR/results/simulation/cross_group_comparison/dirfactor/", "otu", notu, "/")
J_grid = c(33)
nj_grid = c(10, 20)
lambda_grid = c(10)
signal_grid = c(0.5, 0.75, 1, 2, 4)
ntop_grid = c(20)
r = 5
niter = 100000

for (J in J_grid){
  for (nj in nj_grid){
    f0 = paste0(result_dir, "r",r,"_niter",niter,"_J",J, "nj", nj, "H0_teststat.rds")
    if (file.exists(f0)){
      n0 = length(readRDS(f0))
      for (H in c('H1_single', 'H1_multi')){
      for (ntop in ntop_grid){
        for (signal in signal_grid){
          f1 = paste0(result_dir, "r",r,"_niter",niter,"_J",J, "nj", nj, "signal", signal, "ntop", ntop, H, "_teststat.rds")
          if (file.exists(f1)){
            n1 = length(readRDS(f1))
            save_to = paste0(result_dir, "r",r,"_niter",niter,"_J",J, "nj", nj, "signal", signal, "ntop", ntop, H, "_roc_n0_", n0, "_n1_", n1, ".rds")
            # if (!file.exists(save_to)){
              system(paste0(roc_script, " --f0 ", f0, " --f1 ", f1, " --save_to ", save_to))
              roc = readRDS(save_to)
            # }
          }
        }
      }
      }
    }
  }
}
