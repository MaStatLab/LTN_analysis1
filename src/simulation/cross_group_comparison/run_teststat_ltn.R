#!/usr/bin/env Rscript

rm(list=ls())
nseed = 500
niter = 10000

for (J in c(33, 10)){
  for (nj in c(5, 10, 20)){
    for (lambda in c(1, 10)){
      # H0
      system(paste0("REPLACE_WITH_WORKDIR/src/simulation/cross_group_comparison/collect_teststat_ltn.R --WORK_DIR REPLACE_WITH_WORKDIR --INPUT_DIR REPLACE_WITH_WORKDIR/cache/cross_group_comparison/otu100/LTNoutput/ --J ", J," --nj ", nj," --niter ",niter," --lambda ",lambda," --H H0 --nseed ", nseed))
      # H1 single
      for (ntop_single in c(5, 20)){
        for (signal_single in c(2, 4)){
          system(paste0("REPLACE_WITH_WORKDIR/src/simulation/cross_group_comparison/collect_teststat_ltn.R --WORK_DIR REPLACE_WITH_WORKDIR --INPUT_DIR REPLACE_WITH_WORKDIR/cache/cross_group_comparison/otu100/LTNoutput/ --J ", J," --nj ", nj," --niter ",niter," --lambda ",lambda," --H H1_single --nseed ", nseed, " --ntop_single ", ntop_single, " --signal_single ", signal_single))
        }
      }
      # H1 multi
      signal_multi = 0.5
      ntop_multi = 20
          system(paste0("REPLACE_WITH_WORKDIR/src/simulation/cross_group_comparison/collect_teststat_ltn.R --WORK_DIR REPLACE_WITH_WORKDIR --INPUT_DIR REPLACE_WITH_WORKDIR/cache/cross_group_comparison/otu100/LTNoutput/ --J ", J," --nj ", nj," --niter ",niter," --lambda ",lambda," --H H1_multi --nseed ", nseed, " --ntop_multi ", ntop_multi, " --signal_multi ", signal_multi))
    }
  }
}
