
### to calculate F1 scores
get_F1=function(prediction,true){
  retrieved <- sum(prediction)
  precision <- sum(prediction & true) / retrieved
  recall <- sum(prediction & true) / sum(true)
  F1 = 2 * precision * recall / (precision + recall)
  return(F1)
}

get_summary = function(x){
  res=list()
  for (i in 1:length(x)){
    res[[i]]=as.data.frame(rbindlist(lapply(x[[i]], function(a)  reshape2::melt(a))))
  }
  return(res)
}


#############################################
######## compare ARI and time in real data
###########################################
get_ari_time = function(res){

  aris = c(res$IMPACC$ARI_hc,res$IMPACC$ARI_spec,
           res$MPCC$ARI_hc,res$MPCC$ARI_spec,
           res$consensus$ARI_hc,res$consensus$ARI_spec,
           res$sparseKM$ARI,
           res$sparseHclust$ARI,
           res$regular$ARI)
  time = c(res$IMPACC$time,res$MPCC$time,
           res$consensus$time,
           res$sparseKM$time,
           res$sparseHclust$time,
           res$regular$time)
  if (length(res$sparseHclust)>0){
  names(aris)[1:8]=c('IMPACC(HC)','IMPAC (Spec)','MPCC(HC)','MPCC(Spec)',
                     'Consensus (HC)','Consensus (Spec)','sparseKM','sparseHC')
  names(time)[1:5] = c('MPCC','IMPACC','Consensus','sparseKM','sparseHC')
  }else{
    names(aris)[1:7]=c('IMPACC(HC)','IMPAC (Spec)','MPCC(HC)','MPCC(Spec)',
                       'Consensus (HC)','Consensus (Spec)','sparseKM')
    names(time)[1:4] = c('MPCC','IMPACC','Consensus','sparseKM')
  }
  
  
  return(list(aris=aris,time=time))
}



impo_plot = function(impo){
  rownames(impo) =seq(1,nrow(impo))
  colnames(impo) = paste0('Feature',seq(1,ncol(impo)))
  dd = reshape2::melt(impo)
  dd$f=rep(c(rep(TRUE,25),rep(FALSE,4975)),each=530)
  wg = ggplot(data=dd, aes(x=Var1, y=value,group=Var2))+
    geom_line(color='red')+ 
    xlab('# of iteration')+
    ylab('Feature score')+
    geom_vline(xintercept = 30,size=0.5,linetype = "dashed",colour = 'green')+
    theme(axis.text=element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_text(),
          plot.title = element_text(hjust = 0.5))+
    gghighlight(f, use_direct_label = FALSE,
                unhighlighted_params = list( colour = alpha("grey", 0.4)))
  return(wg)
}
fit_cluster = function(mydata,K){
  d0 = dist(mydata,method = 'manhattan')
  fit <- hclust(d0,method = 'ward.D')
  cls = cutree(fit, k=K)
  return(cls)
}

scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- rowSds(data, center = cm)
  csd[which(csd==0)]=0.001
  (data - cm) / csd
}

library('org.Hs.eg.db')
library(edgeR)
get_top_kegg=function(feature_importance,gene_names){
  cut = mean(feature_importance)+sd(feature_importance)
  selected_features = gene_names[which(feature_importance>cut)]
  selected_features=mapIds(org.Hs.eg.db, selected_features, 'ENTREZID', 'SYMBOL')
  keg <- kegga(selected_features, species="Hs")
  return(topKEGG(keg, n=5))
}
