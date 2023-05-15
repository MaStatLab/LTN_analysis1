#!/usr/bin/env Rscript

library(ROCR)
library(GetoptLong)

GetoptLong(
  "f0=s", "file name of test statistics under H0",
  "f1=s", "file name of test statistics under H1",
  "save_to=s", "save to this file"
)

make_roc = function(t0, t1){
  pred = c(t0, t1)
  labels = c(rep(0, length(t0)), rep(1, length(t1)))
  roc = performance(prediction(pred, labels), "tpr", "fpr")
  auc <- performance(prediction(pred, labels), measure = "auc")@y.values[[1]]
  return(list(roc = roc, auc = auc))
}

t0 = readRDS(f0)
t1 = readRDS(f1)
roc = make_roc(t0 = t0, t1 = t1)
saveRDS(roc, save_to)

