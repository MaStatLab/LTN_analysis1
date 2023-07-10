#!/usr/bin/env Rscript

library(Maaslin2)
library(phyloseq)
library(ape)

rm(list = ls())

plot_pmap_temp <- function(pmap, tree, main.text,
                           alpha = NULL, label = NULL, label_nodes = NULL,
                           tip_label = NULL, main.cex = 2.5,
                           cols = c("white", "red")) {
  col_Pal <- grDevices::colorRampPalette(cols)
  graphics::layout(t(1:2), widths = c(96, 4))
  graphics::par(mar = rep(0.5, 4), oma = c(1, 0.5, 2, 2), las = 1)
  K <- length(tree$tip.label)
  node_col <- col_Pal(500)[as.numeric(cut(c((pmap), 0, 1), breaks = 500))]
  if (is.null(tip_label)) {
    tree$tip.label <- rep("", K)
  } else {
    tree$tip.label <- tip_label
  }
  graphics::plot(tree, main = "", show.node.label = FALSE, direction = "downwards", show.tip.label = TRUE, cex.main = 1, use.edge.length = F, align.tip.label = T, node.depth = 2)
  graphics::mtext(main.text, side = 3, cex = main.cex)
  ape::nodelabels(bg = node_col, frame = "none", cex = 3, pch = 21, col = "black")
  if (!is.null(label)) {
    ape::nodelabels(text = label, frame = "none", cex = 1.2)
  }
  if (!is.null(label_nodes) & is.null(label)) {
    label <- rep("", K - 1)
    label[label_nodes] <- paste0("A", 1:length(label_nodes))
    ape::nodelabels(text = label, frame = "none", cex = 1.2)
  }
  if (is.null(alpha)) {
    alpha <- rep(0, K - 1)
  }
  edges <- data.frame(tree$edge)
  edges$leftSign <- ""
  edges$rightSign <- ""
  leftBranch <- sapply((K + 1):(2 * K - 1), function(n) {
    min(which(edges[, 1] == n))
  })
  rightBranch <- sapply((K + 1):(2 * K - 1), function(n) {
    max(which(edges[, 1] == n))
  })
  edges[leftBranch, "leftSign"] <- c("-", "+")[as.numeric(alpha > 0) + 1]
  edges[rightBranch, "rightSign"] <- c("-", "+")[as.numeric(alpha < 0) + 1]
  for (j in seq_along(alpha)) {
    if (alpha[j] == 0) {
      edges[leftBranch[j], "leftSign"] <- ""
      edges[rightBranch[j], "rightSign"] <- ""
    }
  }
  ape::edgelabels(text = edges$leftSign, frame = "none", cex = 1.2, adj = 1)
  ape::edgelabels(text = edges$rightSign, frame = "none", cex = 1.2, adj = -0.2)
  legend_image <- grDevices::as.raster(matrix(col_Pal(500), ncol = 1))
  graphics::image(z = t(1:500), col = legend_image, axes = FALSE)
  graphics::mtext("PMAP", side = 3, cex = 1)
  graphics::axis(side = 4, cex.axis = 0.8, tick = T)
  # check sign
  lr <- sapply((K + 1):(2 * K - 1), function(n) {
    which(edges[, 1] == n)
  })
  try(if (!all(apply(rbind(lr, alpha), 2, function(x) {
    (edges[x[1], "leftSign"] != edges[x[2], "rightSign"] &
      edges[x[1], "leftSign"] != "" &
      edges[x[2], "rightSign"] != "") |
      (edges[x[1], "leftSign"] == edges[x[2], "rightSign"] &
        edges[x[1], "leftSign"] == "" &
        x[3] == 0)
  }))) {
    stop(paste0("inconsistent left and right sign", which(!(
      apply(rbind(lr, alpha), 2, function(x) {
        (edges[x[1], "leftSign"] != edges[x[2], "rightSign"] &
          edges[x[1], "leftSign"] != "" &
          edges[x[2], "rightSign"] != "") |
          (edges[x[1], "leftSign"] == edges[x[2], "rightSign"] &
            edges[x[1], "leftSign"] == "" &
            x[3] == 0)
      })
    )), "\n"))
  })
  # left and right sign being exactly opposite
  alpha_sign <- c("-", "+")[as.numeric(alpha > 0) + 1]
  alpha_sign[alpha == 0] <- ""
  try(if (!all(sapply(1:length(alpha), function(x) {
    alpha_sign[x] == edges[which(edges[, 1] == x + K)[1], "leftSign"]
  }))) {
    stop(paste0("inconsistent left and node sign", which((
      sapply(1:length(alpha), function(x) {
        alpha_sign[x] != edges[which(edges[, 1] == x + K)[1], "leftSign"]
      })
    )), "\n"))
  })
  # left sign same as alpha sign
}

ps <- readRDS("cache/ps_otu_100.rds")
input_data <- data.frame(otu_table(ps)) # taxa are columns
covariates <- c(
  "Case_Control", "post_seroconversion", "Age_at_Collection", "Country", "Gender",
  "BF", "Solid_Food", "Eggs", "Fish", "Soy_Prod", "Rye", "Barley", "Buckwheat_Millet"
)
# baseline = c("control",F,rep("false",length(dietcovariates)))
# covariate = "BF"
# formula = as.formula(paste0("~",covariate,"+log(Age_at_Collection)+Country+post_seroconversion+Gender+Case_Control+",paste(dietcovariates[-which(dietcovariates==covariate)],collapse = "+"),"+(1 | Subject_ID)"))
sample_data <- data.frame(sample_data(ps))[, c(covariates, "Subject_ID")]
sample_data$log_age <- log(sample_data$Age_at_Collection)
sample_data$Age_at_Collection <- NULL

output_dir <- "default_output"
fixed_effects <- gsub("Age_at_Collection", "log_age", covariates)

taxtab <- data.frame(tax_table(ps))
taxtab$OTU <- paste0("X", rownames(taxtab))
all_results <- read.csv(paste0(output_dir, "/all_results.tsv"), sep = "\t")
groups <- c(
  "Case_Control", "post_seroconversion",
  "BF", "Solid_Food", "Eggs", "Fish", "Soy_Prod",
  "Rye", "Barley", "Buckwheat_Millet"
)


# calculate BH q-value for each of the groups
signif_otu <- list()
for (i in seq_along(groups)) {
  covariate <- groups[i]
  # subset: this group only
  results_single_covariate <- all_results[all_results$metadata == covariate, ]
  # print(dim(results_single_covariate))
  # BH correction
  results_single_covariate$qval <- p.adjust(results_single_covariate$pval,
    method = "BH"
  )
  # merge results and tax table
  signif_taxa <- merge(
    x = results_single_covariate, y = taxtab,
    by.x = "feature", by.y = "OTU"
  )
  # print(sum(signif_taxa$metadata != covariate))
  # signif_taxa <- signif_taxa[order(signif_taxa$metadata), ]
  # write.csv(signif_taxa, paste0(output_dir, "/signif_taxa.csv"))
  # print(sum(signif$qval <= 0.05))
signif_otu[[covariate]] <- signif_taxa[signif_taxa$metadata == covariate, ]
}

# visualization
source("src/LTN/utils.R")
source("bayesian_fdr.R")

result_dir <- "results/application/pmap/"

# pdf(paste0("results/application/fig/pmap_qval_m005_BHgroup.pdf"), width = 14)
tree <- phy_tree(ps)
rename_fixed_effects <- c("Case_Control", "Post_seroconversion", "log_age", "Country", "Gender", "Breastfeeding", "Solid_Food", "Eggs", "Fish", "Soy_Product", "Rye", "Barley", "Buckwheat_Millet")
names(rename_fixed_effects) <- fixed_effects
for (x in groups) {
  pdf(paste0("results/application/fig/pmap_qval_m005_BHgroup",rename_fixed_effects[x],".pdf"), width = 14)
  par(cex = 0.5)
  file <- paste0(result_dir, x, "_lambda10_m005.rds")
  if (file.exists(file)) {
    pmap <- readRDS(file)
    pjap <- readRDS(paste0("results/application/pjap/", x, "_lambda10_m005.rds"))
    pmap_threshold <- bfdr_control(pmap, 0.05)
    n_signif_nodes <- sum(pmap >= pmap_threshold)
    signif_nodes <- which(pmap >= pmap_threshold)
    node_label <- NULL
    node_label005 <- NULL
    node_label01 <- NULL
    if (n_signif_nodes > 0) {
      signif_nodes <- order(pmap, decreasing = T)[1:n_signif_nodes]
      node_label005 <- c("", "*")[(1:99 %in% signif_nodes) + 1]
    }
    pmap_threshold <- bfdr_control(pmap, 0.1)
    n_signif_nodes <- sum(pmap >= pmap_threshold)
    signif_nodes <- which(pmap >= pmap_threshold)
    if (n_signif_nodes > 0) {
      signif_nodes <- order(pmap, decreasing = T)[1:n_signif_nodes]
      node_label01 <- c("", "*")[(1:99 %in% signif_nodes) + 1]
      node_label <- paste(node_label005, node_label01, sep = "")
    }
    # tip_label = sapply(paste0("X", tree$tip.label) %in% signif_otu[[x]], function(x){if(x){"o"}else{""}})
    tip_label <- sapply(paste0("X", tree$tip.label), function(otu) {
      if (!(otu %in% signif_otu[[x]]$feature)) {
        return("")
      } else {
        qval <- signif_otu[[x]][signif_otu[[x]]$feature == otu, "qval"]
        mark <- ""
        if (qval < 0.1) {
          mark <- paste0(mark, "*")
        }
        if (qval < 0.05) {
          mark <- paste0(mark, "*")
        }
        return(mark)
      }
    })
    # main.text = paste0(x, ", PJAP:", round(pjap, 4), ", PMAP threshold: ", round(pmap_threshold, 2))
    # main.text <- paste0(rename_fixed_effects[x])# , ", PJAP=", round(pjap, 4))
    main.text <- ""
    plot_pmap_temp(pmap, tree, main.text = main.text, tip_label = tip_label, label = node_label, main.cex = 1.5, cols = c("#55b6fc", "yellow"))
  } else {
    warning(paste0(x, ": file does not exist"))
  }
  dev.off()
}


# sanity check
# dev.off()
# pdf('tmp.pdf')
# for (x in groups){
# plot(rank(all_results[all_results$metadata == x, 'qval']), 
# rank(p.adjust(all_results[all_results$metadata == x, 'pval'], 
#   method = "BH"
# )), main = x)
# plot(
# rank(signif_otu[[x]][,'pval']),rank(signif_otu[[x]][,'qval'])
# )
# aaa = all_results[all_results$metadata == x, ]
# plot((signif_otu[[x]][order(signif_otu[[x]][,'feature']),'qval']), (aaa[order(aaa$feature),'qval'])) 
# plot(rank(signif_otu[[x]][order(signif_otu[[x]][,'feature']),'qval']), rank(aaa[order(aaa$feature),'qval']))
# }
# dev.off()
