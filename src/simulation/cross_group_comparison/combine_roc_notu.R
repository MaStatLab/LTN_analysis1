#!/usr/bin/env Rscript

rm(list = ls())
notu = 50
J_grid = c(33)
nj_grid = c(10, 20)
lambda_grid = c(10)
# signal_grid = c(0.5, 1, 2, 4)
ntop_grid = c(20)
models = c("LTN", "MaAsLin2", "dirfactor")
resdir = paste0("REPLACE_WITH_WORKDIR/results/simulation/cross_group_comparison/", models, "/otu", notu)
names(resdir) = models
cols = c('black', 'orange', 'blue')

# single
H = "H1_single"
signal_grid = c(0.5, 2, 4)
pdf(paste0("REPLACE_WITH_WORKDIR/results/roc_all_otu", notu, H, format(Sys.time(), "_%H:%M_%b_%e_%Y"), ".pdf"), height = 6, width = 8)
par(mfrow = c(2,3))
par(pty = 's')
for (J in J_grid){
  for (nj in nj_grid){
    for (signal in signal_grid){
      for (ntop in ntop_grid){
        # plot_title = paste0("nj=",nj,", count x",signal+1, " at 1 OTU")
        plot_title = paste0('a* = ', signal)
        str1 = paste0("J",J,"nj",nj,"signal",signal,"ntop",ntop,H,"_roc")
        plot(1, col = "white", xlim = c(0,1), ylim = c(0,1), main = plot_title, xlab = "FPR", ylab = "TPR")
        abline(a = 0, b = 1, col = "grey")
        auc = rep(NA, 3)
        if (length(grep(str1, list.files(resdir["LTN"]))) > 0){
          fs = grep(str1, list.files(resdir["LTN"], full.names = T), val = T)
          f10 = grep("lambda10\\_", fs, val = T)
          f1 = grep("lambda1\\_", fs, val = T)
          roc_auc = readRDS(f10)
          if (length(f10) > 0){
            plot(roc_auc[['roc']], add = T, col = cols[1], lty = 1)
            auc[1] = roc_auc[['auc']]
          }
          # if (length(f1) > 0){
          #   plot(readRDS(f1), add = T, col = 4)
          # }
        }
        if (length(grep(str1, list.files(resdir["MaAsLin2"]))) > 0){
          roc_auc = readRDS(grep(str1, list.files(resdir["MaAsLin2"], full.names = T), val = T))
          plot(roc_auc[['roc']], add = T, col = cols[2], lty = 2)
          auc[2] = roc_auc[['auc']]
        }
        if (length(grep(str1, list.files(resdir["dirfactor"]))) > 0){
          roc_auc = readRDS(grep(str1, list.files(resdir["dirfactor"], full.names = T), val = T))
          plot(roc_auc[['roc']], add = T, col = cols[3], lty = 3)
          auc[3] = roc_auc[['auc']]
        }
        legend("bottomright", c("LTN", "MaAsLin2", "DirFactor"), col = cols, lty = 1:3)
        if (!is.null(auc)){        
         # legend("bottomright", title = 'AUC', legend = paste(c("LTN", "MaAsLin2", "DirFactor"), round(auc, 2), sep = ":"), col = cols, lty = 1:3)
        }
      }
    }
    mtext(paste0("n* = ", nj), side = 4, line = 1, at = 0.5)
  }
}  
dev.off()



# multi
H = "H1_multi"
signal_grid = c(0.5, 0.75, 1)
pdf(paste0("REPLACE_WITH_WORKDIR/results/roc_all_otu", notu, H, format(Sys.time(), "_%H:%M_%b_%e_%Y"), ".pdf"), height = 6, width = 8)
par(mfrow = c(2,3))
par(pty = 's')
for (J in J_grid){
  for (nj in nj_grid){
    for (signal in signal_grid){
      for (ntop in ntop_grid){
        # plot_title = paste0("nj=",nj,", count x",signal+1, " at 8 OTUs")
        plot_title = paste0('a* = ', signal)
        str1 = paste0("J",J,"nj",nj,"signal",signal,"ntop",ntop,H,"_roc")
        plot(1, col = "white", xlim = c(0,1), ylim = c(0,1), main = plot_title, xlab = "FPR", ylab = "TPR")
        abline(a = 0, b = 1, col = "grey")
        auc = rep(NA, 3)
        if (length(grep(str1, list.files(resdir["LTN"]))) > 0){
          fs = grep(str1, list.files(resdir["LTN"], full.names = T), val = T)
          f10 = grep("lambda10\\_", fs, val = T)
          f1 = grep("lambda1\\_", fs, val = T)
          if (length(f10) > 0){
            roc_auc = readRDS(f10)
            plot(roc_auc[['roc']], add = T, col = cols[1], lty = 1)
            auc[1] = roc_auc[['auc']]
          }
          # if (length(f1) > 0){
          #   plot(readRDS(f1), add = T, col = 4)
          # }
        }
        if (length(grep(str1, list.files(resdir["MaAsLin2"]))) > 0){
          roc_auc = readRDS(grep(str1, list.files(resdir["MaAsLin2"], full.names = T), val = T))
          plot(roc_auc[['roc']], add = T, col = cols[2], lty = 2)
          auc[2] = roc_auc[['auc']]
        }
        if (length(grep(str1, list.files(resdir["dirfactor"]))) > 0){
          roc_auc = readRDS(grep(str1, list.files(resdir["dirfactor"], full.names = T), val = T))
          plot(roc_auc[['roc']], add = T, col = cols[3], lty = 3)
          auc[3] = roc_auc[['auc']]
        }
        legend("bottomright", c("LTN", "MaAsLin2", "DirFactor"), col = cols, lty = 1:3)
        if (!is.null(auc)){
          # legend("bottomright", title = 'AUC', legend = paste(c("LTN", "MaAsLin2", "DirFactor"), round(auc, 2), sep = ":"), col = cols, lty = 1:3)
        }
      }
    }
    mtext(paste0("n* = ", nj), side = 4, line = 1, at = 0.5)
  }
}
dev.off()
