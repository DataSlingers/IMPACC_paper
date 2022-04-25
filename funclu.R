### clustering 
fit_cluster = function(mydata,K){
  d0 = dist(mydata,method = 'manhattan')
  fit <- hclust(d0,method = 'ward.D')
  cls = cutree(fit, k=K)
  return(cls)
}
get_F1=function(prediction,true){
  retrieved <- sum(prediction)
  precision <- sum(prediction & true) / retrieved
  recall <- sum(prediction & true) / sum(true)
  F1 = 2 * precision * recall / (precision + recall)
  return(F1)
}
get_F1_more=function(prediction,true){
  TP = sum(true & prediction)
  TN = sum(!true & !prediction)
  FP = sum(!true & prediction)
  FN = sum(true & !prediction)
  Recall = TP / sum(true)
  Precision = TP / (TP + FP)
  F1 = 2 * ((Precision * Recall) /
              (Precision + Recall))
  Accuracy =  sum(TP + TN) /length(true)
  FPR = FP / (FP + TN)
  FDR = FP/(FP+TP)
  return(c(Accuracy,Precision,Recall,F1,FPR,FDR))
}

  
library(matrixStats)
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- rowSds(data, center = cm)
  csd[which(csd==0)]=0.001
  (data - cm) / csd
}

