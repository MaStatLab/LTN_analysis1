
# setup
source("REPLACE_WITH_WORKDIR/src/LTN/utils.R")
result_dir = "REPLACE_WITH_WORKDIR/results/application/"

# group list
dietcovariates=c("BF","Solid_Food","Eggs","Fish","Soy_Prod","Rye","Barley","Buckwheat_Millet")
covariates=c("Case_Control","post_seroconversion",dietcovariates)
# PJAPs of groups
pjap_files = paste0(result_dir, "pjap/", covariates, "_lambda10_m005.rds")
pjap = sapply(pjap_files,readRDS)
names(pjap) = gsub("BF", "Breastfeeding", covariates)
cat("prior joint alternative probability: ", 1 - (1 - 0.05)^99, "\n")
print(pjap)

# PMAPs of groups
pmap_files = paste0(result_dir, "pmap/", covariates, "_lambda10_m005.rds")
pmap = sapply(pmap_files,readRDS)
colnames(pmap) = names(pjap)


