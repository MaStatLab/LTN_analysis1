#!/usr/bin/env Rscript

rm(list = ls())
WORK_DIR = 'REPLACE_WITH_WORKDIR/'
notu = 100 # total number of OTUs
batch_id = "multi075" # identifier for this round of experiment
filedir = paste0(WORK_DIR, "/cache/cross_group_comparison/otu", notu, "/sim/") 

# raw = grep("H", list.files(paste0(filedir, "/H0"), full.names = T, recursive = T), val = T)
raw = sort(c(grep("H", list.files(paste0(filedir, "/H0"), full.names = T, recursive = T), val = T),
               grep("H", list.files(paste0(filedir, "/H1_single"), full.names = T, recursive = T), val = T),
               grep("H", list.files(paste0(filedir, "/H1_multi"), full.names = T, recursive = T), val = T)
))
# raw = raw[-grep('backup', raw)]
# raw = raw[-grep('old', raw)]
# raw = raw[-grep('single', raw)]
# raw = raw[grep('J33nj20', raw)]
# raw = raw[-grep('signal2', raw)]
# raw = raw[-grep("ntop_multi20signal_multi0\\.5", raw)]
# raw = raw[-grep('signal1', raw)]
# raw = raw[-grep('signal2', raw)]
# raw = raw[-grep('signal0\\.75', raw)]
raw = c(raw[grep('multi\\/J33nj10\\/ntop20signal0.75', raw)], 
        raw[grep('multi\\/J33nj20\\/ntop20signal0.75', raw)]
        )
input = gsub("\\/sim\\/", "\\/LTNinput\\/", raw)
output_dir = gsub("\\.rds", "", gsub("\\/sim\\/", "\\/LTNcache\\/", raw))
output_pjap_dir = sapply(raw, function(x){gsub("\\/sim\\/", "\\/LTNoutput\\/", dirname(x))})
output_pjap_file = sapply(raw, function(x){paste0("LTNoutput_", basename(x))})

config_file = paste0(WORK_DIR, "/cache/cross_group_comparison/otu",notu,"/ltn_config_notu",notu,"batch",batch_id,".csv")
write.csv(data.frame(raw, input, output_dir, output_pjap_dir, output_pjap_file), file = config_file)
