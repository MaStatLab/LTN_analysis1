#!/usr/bin/env Rscript
rm(list = ls())
nseed = 500
niter = 100000
r = 5

for (J in c(33, 10)){
  for (nj in c(5, 10, 20)){
    # H0
    try(system(paste0("REPLACE_WITH_WORKDIR/src/simulation/cross_group_comparison/collect_teststat_dirfactor.R --WORK_DIR REPLACE_WITH_WORKDIR --INPUT_DIR REPLACE_WITH_WORKDIR/cache/cross_group_comparison/otu100/dirfactor/ --J ", J," --nj ", nj," --H H0 --nseed ", nseed, " --r ", r, " --niter ", format(niter, scientific = F))))
    # H1 single
    for (ntop_single in c(5, 20)){
      for (signal_single in c(2, 4)){
        try(system(paste0("REPLACE_WITH_WORKDIR/src/simulation/cross_group_comparison/collect_teststat_dirfactor.R --WORK_DIR REPLACE_WITH_WORKDIR --INPUT_DIR REPLACE_WITH_WORKDIR/cache/cross_group_comparison/otu100/dirfactor/ --J ", J," --nj ", nj," --H H1_single --nseed ", nseed, " --ntop_single ", ntop_single, " --signal_single ", signal_single, " --r ", r, " --niter ", format(niter, scientific = F))))
      }
    }
    # H1 multi
    signal_multi = 0.5
    ntop_multi = 20
    try(system(paste0("REPLACE_WITH_WORKDIR/src/simulation/cross_group_comparison/collect_teststat_dirfactor.R --WORK_DIR REPLACE_WITH_WORKDIR --INPUT_DIR REPLACE_WITH_WORKDIR/cache/cross_group_comparison/otu100/dirfactor/ --J ", J," --nj ", nj," --H H1_multi --nseed ", nseed, " --ntop_multi ", ntop_multi, " --signal_multi ", signal_multi, " --r ", r, " --niter ", format(niter, scientific = F))))
  }
}
