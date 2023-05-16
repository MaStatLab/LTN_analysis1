#!/usr/bin/env Rscript

library(Maaslin2)
library(phyloseq)

rm(list = ls())

WORK_DIR = Sys.getenv("WORK_DIR")
plot_pmap_temp = function(pmap,tree,main.text,alpha=NULL,label=NULL,label_nodes=NULL,tip_label=NULL, main.cex = 2.5, cols = c('white', 'red')){
  col_Pal=grDevices::colorRampPalette(cols)
  graphics::layout(t(1:2), widths=c(96,4))
  graphics::par(mar=rep(0.5, 4), oma=c(1,0.5,2,2), las=1)
  K=length(tree$tip.label)
  node_col = col_Pal(500)[as.numeric(cut(c((pmap),0,1), breaks = 500)) ]
  if (is.null(tip_label)){
    tree$tip.label=rep('',K)
  }else{
    tree$tip.label=tip_label
  }
  graphics::plot(tree,main='',show.node.label=FALSE,direction="downwards",show.tip.label = TRUE,cex.main=1,use.edge.length=F,align.tip.label=T,node.depth=2)
  graphics::mtext(main.text,side=3,cex=main.cex)
  ape::nodelabels(bg=node_col,frame='none',cex=3,pch=21,col ='black')
  if (!is.null(label)){
    ape::nodelabels(text=label,frame='none',cex=1.2)
  }
  if (!is.null(label_nodes) & is.null(label)){
    label=rep('',K-1)
    label[label_nodes]=paste0('A',1:length(label_nodes))
    ape::nodelabels(text=label,frame='none',cex=1.2)
  }
  if (is.null(alpha)){
    alpha=rep(0,K-1)
  }
  edges=data.frame(tree$edge)
  edges$leftSign=''
  edges$rightSign=''
  leftBranch=sapply((K+1):(2*K-1),function(n){min(which(edges[,1]==n))})
  rightBranch=sapply((K+1):(2*K-1),function(n){max(which(edges[,1]==n))})
  edges[leftBranch,'leftSign']=c('-','+')[as.numeric(alpha>0)+1]
  edges[rightBranch,'rightSign']=c('-','+')[as.numeric(alpha<0)+1]
  for (j in seq_along(alpha)){
    if (alpha[j]==0){
      edges[leftBranch[j],'leftSign']=''
      edges[rightBranch[j],'rightSign']=''
    }
  }
  ape::edgelabels(text=edges$leftSign,frame='none',cex=1.2,adj = 1)
  ape::edgelabels(text=edges$rightSign,frame='none',cex=1.2,adj = -0.2)
  legend_image <- grDevices::as.raster(matrix(col_Pal(500), ncol=1))
  graphics::image(z=t(1:500), col=legend_image, axes=FALSE)
  graphics::mtext('PMAP',side=3,cex=1)
  graphics::axis(side=4,cex.axis=0.8,tick=T)
  # check sign
  lr=sapply((K+1):(2*K-1),function(n){which(edges[,1]==n)})
  try(if (!all(apply(rbind(lr,alpha), 2, function(x) {
    (edges[x[1], 'leftSign'] != edges[x[2], 'rightSign'] &
     edges[x[1], 'leftSign'] != '' &
     edges[x[2], 'rightSign'] != '')|
      (edges[x[1], 'leftSign'] == edges[x[2], 'rightSign'] &
       edges[x[1], 'leftSign'] ==''&
       x[3]==0)
  })))
    stop(paste0("inconsistent left and right sign", which(!(
      apply(rbind(lr,alpha), 2, function(x) {
        (edges[x[1], 'leftSign'] != edges[x[2], 'rightSign'] &
           edges[x[1], 'leftSign'] != '' &
           edges[x[2], 'rightSign'] != '')|
          (edges[x[1], 'leftSign'] == edges[x[2], 'rightSign'] &
             edges[x[1], 'leftSign'] ==''&
             x[3]==0)
      })
    )), '\n')))
  #left and right sign being exactly opposite
  alpha_sign=c('-','+')[as.numeric(alpha>0)+1]
  alpha_sign[alpha==0]=''
  try(if (!all(sapply(1:length(alpha), function(x) {
    alpha_sign[x] == edges[which(edges[, 1] == x + K)[1], 'leftSign']
  })))
    stop(paste0('inconsistent left and node sign', which((
      sapply(1:length(alpha), function(x) {
        alpha_sign[x] != edges[which(edges[, 1] == x + K)[1], 'leftSign']
      })
    )),'\n')))
  # left sign same as alpha sign
}

ps = readRDS(paste0(WORK_DIR,"/cache/ps_otu_100.rds"))
input_data = data.frame(otu_table(ps)) # taxa are columns
covariates = c("Case_Control", "post_seroconversion", "Age_at_Collection", "Country", "Gender",
"BF", "Solid_Food", "Eggs", "Fish", "Soy_Prod", "Rye", "Barley", "Buckwheat_Millet")
# baseline = c("control",F,rep("false",length(dietcovariates)))
# covariate = "BF"
# formula = as.formula(paste0("~",covariate,"+log(Age_at_Collection)+Country+post_seroconversion+Gender+Case_Control+",paste(dietcovariates[-which(dietcovariates==covariate)],collapse = "+"),"+(1 | Subject_ID)"))
sample_data = data.frame(sample_data(ps))[,c(covariates, "Subject_ID")]
sample_data$log_age = log(sample_data$Age_at_Collection)
sample_data$Age_at_Collection = NULL

output_dir = "default_output"
fixed_effects = gsub("Age_at_Collection", "log_age", covariates)

if (!file.exists(paste0(output_dir, "/significant_results.tsv"))){
    test_fit = Maaslin2(input_data = input_data,
input_metadata = sample_data, 
output = output_dir, 
fixed_effects = fixed_effects,
random_effects = "Subject_ID",
plot_heatmap = F, plot_scatter = F) 
}


signif = read.csv(paste0(output_dir, "/significant_results.tsv"), sep = "\t")
taxtab = data.frame(tax_table(ps))
taxtab$OTU = paste0("X", rownames(taxtab))
signif_taxa = merge(x = signif, y = taxtab, by.x = "feature", by.y = "OTU")
signif_taxa = signif_taxa[order(signif_taxa$metadata),]
if (!file.exists(paste0(output_dir, "/signif_taxa.csv"))){
    write.csv(signif_taxa, paste0(output_dir, "/signif_taxa.csv"))
}
signif_otu = list()
for (x in fixed_effects){
    signif_otu[[x]] = signif_taxa[signif_taxa$metadata == x, ]
}

# visualization
library(ape)
source(paste0(WORK_DIR, "/src/LTN/utils.R"))
source("bayesian_fdr.R")
# pdf(paste0(output_dir,"/pmap_qval.pdf"), width = 14)

result_dir = paste0(WORK_DIR, "results/application/pmap/")
pdf(paste0(result_dir,"/pmap_qval_m005_new.pdf"), width = 14)
par(cex = 0.5)
rename_fixed_effects = c('Case_Control','Post_seroconversion','log_age','Country','Gender','Breastfeeding','Solid_Food','Eggs','Fish','Soy_Product','Rye','Barley','Buckwheat_Millet')
names(rename_fixed_effects) = fixed_effects
for (x in fixed_effects){
    file = paste0(result_dir, x, "_lambda10_m005.rds")
    if (file.exists(file)){
        pmap = readRDS(file)
        pjap = readRDS(paste0("results/application/pjap/", x, "_lambda10_m005.rds"))
        tree = phy_tree(ps)
        pmap_threshold = bfdr_control(pmap, 0.05)
        n_signif_nodes = sum(pmap >= pmap_threshold)
        signif_nodes = which(pmap >= pmap_threshold)
        node_label = NULL
        node_label005 = NULL
        node_label01 = NULL
        if (n_signif_nodes > 0){
            signif_nodes = order(pmap, decreasing = T)[1:n_signif_nodes]
            node_label005 = c('','*')[(1:99 %in% signif_nodes) + 1]
        }
        pmap_threshold = bfdr_control(pmap, 0.1)
        n_signif_nodes = sum(pmap >= pmap_threshold)
        signif_nodes = which(pmap >= pmap_threshold)
        if (n_signif_nodes > 0){
            signif_nodes = order(pmap, decreasing = T)[1:n_signif_nodes]
            node_label01 = c('','*')[(1:99 %in% signif_nodes) + 1]
            node_label = paste(node_label005, node_label01, sep = '')
        }
        # tip_label = sapply(paste0("X", tree$tip.label) %in% signif_otu[[x]], function(x){if(x){"o"}else{""}})
        tip_label = sapply(paste0("X", tree$tip.label), function(otu){
            if (!(otu %in% signif_otu[[x]]$feature)){
                return("")
            }
            else{
                qval = signif_otu[[x]][signif_otu[[x]]$feature==otu, "qval"]
                mark = ""
                # if (qval < 0.25){
                #     mark = paste0(mark, "*")
                # }
                if (qval < 0.1){
                    mark = paste0(mark, "*")
                }
                if (qval < 0.05){
                    mark = paste0(mark, "*")
                }
                return(mark)
            }
        })
        # main.text = paste0(x, ", PJAP:", round(pjap, 4), ", PMAP threshold: ", round(pmap_threshold, 2))
        main.text = paste0(rename_fixed_effects[x], ", PJAP=", round(pjap, 4))
        plot_pmap_temp(pmap, tree, main.text =  main.text, tip_label = tip_label, label = node_label, main.cex = 1.5, cols = c('#55b6fc', 'yellow'))
    }
    else{
        warning(paste0(x, ": file does not exist"))
    }
}
dev.off()

# seroconversion: list of taxa
fdr_level = 0.05
x = "post_seroconversion"
file = paste0(result_dir, x, "_lambda10_m005.rds")
pmap = readRDS(file)
tree = phy_tree(ps)
pmap_threshold = bfdr_control(pmap, fdr_level)
n_signif_nodes = sum(pmap >= pmap_threshold)
signif_nodes = which(pmap >= pmap_threshold)
if (n_signif_nodes > 0){
  signif_nodes = order(pmap, decreasing = T)[1:n_signif_nodes]
}
library(phytools)
otu_ltn = NULL
for (node in signif_nodes){
  otu = getDescendants(tree, node + 100)[getDescendants(tree, node + 100) <= 100]
  df = data.frame(tax_table(ps)[otu, ])
  df$node = node
  otu_ltn = rbind(otu_ltn, df)
}
otu_ltn$feature = paste0('X', rownames(otu_ltn))
otu_ltn$model = 'ltn'
otu_maaslin = signif_otu[[x]][signif_otu[[x]]$qval <= fdr_level,]
otu_maaslin$model = 'maaslin'
otu_merged = merge(otu_ltn, otu_maaslin, by = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'feature'), all.x = T, all.y = T)
saveRDS(otu_merged, paste0(result_dir, 'otu_merged.rds'))
write.csv(otu_merged, paste0(result_dir, 'otu_merged.csv'))
