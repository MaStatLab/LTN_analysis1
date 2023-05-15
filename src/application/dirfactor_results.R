library(glue)

# within_ind v norm -------------------------------------------------------

within_ind = 1
resdir = '/work/zw122/LTN_analysis1/cache/application/dirfactor/teststat/'
files = grep(paste0('within',within_ind,'r5iter1e\\+05_v_post_norm'), list.files(resdir, full.names = T), val = T)
v_test1 = readRDS(grep('i0',files,val=T))
v_null1 = sapply(files[-grep('i0',files,val=F)], readRDS)
mean(v_null1>=v_test1)

# across_ind v norm -------------------------------------------------------

within_ind = 0
resdir = '/work/zw122/LTN_analysis1/cache/application/dirfactor/teststat/'
files = grep(paste0('within',within_ind,'r5iter1e\\+05_v_post_norm'), list.files(resdir, full.names = T), val = T)
v_test0 = readRDS(grep('i0',files,val=T))
v_null0 = sapply(files[-grep('i0',files,val=F)], readRDS)
mean(v_null0 >= v_test0)


# plot -------------------------------------------------------
par(mfrow = c(1,2))
hist(v_null1, xlim = range(c(v_test1,v_null1)), main = glue('permutating within individuals',), xlab = 'test statistics')
abline(v = v_test1, col='red')

hist(v_null0, xlim = range(c(v_test0,v_null0)), main = glue('permutating across all samples',), xlab = 'test statistics')
abline(v = v_test0, col='red')

# within_ind v -------------------------------------------------------

within_ind = 1
resdir = '/work/zw122/LTN_analysis1/cache/application/dirfactor/teststat/'
files = grep(paste0('within',within_ind,'r5iter1e\\+05_v_post\\.rds'), list.files(resdir, full.names = T), val = T)
vec_test = readRDS(grep('i0',files,val=T))
vec_null = sapply(files[-grep('i0',files,val=F)], readRDS)
pval_otu = rowMeans(apply(vec_null,2,function(x) abs(vec_test)<=abs(x)))

taxtab = readRDS("/work/zw122/LTN_analysis1/cache/diabimmune_taxtab_otu_100.rds")
taxtab[which(pval_otu == 0),]
