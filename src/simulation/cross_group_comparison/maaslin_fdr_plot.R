
par(mfrow = c(1,3))
WORK_DIR = Sys.getenv('WORK_DIR')
d = paste0(WORK_DIR,"/results/simulation/cross_group_comparison/MaAsLin2/otu100/")
f = list.files(d, full.names = T)
ff = grep('fdr005_potu', f, val=T)
titles = c('K*=0', 'K*=1', 'K*=8')
for (i in 1:3){
  s = c('H0', 'single', 'multi')[i]
  files = grep(s, ff, val = T)
  fdr_vec = sapply(files, readRDS)
  hist(fdr_vec, main = titles[i], breaks = 20, freq = F, xlab = 'FDR')
  cat(titles[i], ': mean: ', mean(fdr_vec), ', standard deviation: ', sd(fdr_vec), '\n')
}


