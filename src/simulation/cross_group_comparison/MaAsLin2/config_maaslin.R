
rm(list=ls())
WORK_DIR = 'REPLACE_WITH_WORKDIR/'
notu = 100
batch_id = "nj10_multi2"
filedir = paste0(WORK_DIR, "/cache/cross_group_comparison/otu", notu, "/sim/") # 
# input = grep("H", list.files(paste0(filedir, "/H0"), full.names = T, recursive = T), val = T)
input = sort(c(grep("H", list.files(paste0(filedir, "/H0"), full.names = T, recursive = T), val = T),
               grep("H", list.files(paste0(filedir, "/H1_single"), full.names = T, recursive = T), val = T),
               grep("H", list.files(paste0(filedir, "/H1_multi"), full.names = T, recursive = T), val = T)
))
input = c(input[grep('multi\\/J33nj10\\/ntop20signal2', input)]
)
# input = input[-grep('old', input)]
# input = input[-grep('single', input)]
# input = input[grep('J33nj20', input)]
# # input = input[-grep('signal2', input)]
# input = input[-grep("ntop_multi20signal_multi0\\.5", input)]
# input = input[-grep('signal1', input)]
# input = input[-grep('signal2', input)]
# # input = input[-grep('signal4', input)]
# input = input[-grep('signal0\\.75', input)]
output = gsub("\\.rds", "", gsub("\\/sim\\/", "\\/MaAsLin2\\/", input))
config_file = paste0(WORK_DIR, "/cache/cross_group_comparison/otu", notu, "/maaslin_config", batch_id, ".csv")
write.csv(data.frame(input, output), file = config_file)

