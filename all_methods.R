## function to run cluster algorithms in empirical study  
## regular: includes hclust, kmeans, kmedoid
source('impacc.R')
library(matrixStats)
library(sparcl)
library(Matrix)
library(mclust)
library(SC3)
library(Seurat)
library(Rtsne)
library(scater)
library(tictoc)
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- rowSds(data, center = cm)
  csd[which(csd==0)]=0.001
  (data - cm) / csd
}

all_methods = function(dat,nIter=200, m=0.1,IMPACC = T, MPACC =F,MPCC = T,
                       consensus =T,sparseKM = T,
                       sparseHC = T,
                       seurat=T,sc3=T,tsne=T,
                       early_stop=T,
                       regular =T,spec = T,verbose = F,
                       num_unchange = 5,eps = 0.001 ## for early stop
                       ){
  
  sc=dat$raw
  if (is.null(sc)){
    sc= dat$sc_cnt
  }
  sc_scale=dat$sc_cnt
  sc_standardized = scalematrix(as.matrix(sc_scale))
  sc_label = dat$sc_label
  
  if (ncol(dat$sc_cnt)<100){
    pItem = 0.5
  }else{
    pItem = 0.25
  }
  #pItem=0.25
  # pItem = max(50/ncol(biase$sc_cnt),0.25)
  # if (pItem>1){
  #   pItem = 0.9
  # }

  K = length(unique(sc_label))
  res=list()
  if (IMPACC == T){
    print('run IMPACC')
    tic.clearlog()
    tic()    
    css= IMPACC(d=sc_scale,K = K, E=3,reps=nIter,pItem=pItem,pFeature=0.1,early_stop=early_stop,verbose = verbose)
    toc(log = TRUE)

      ### early stopping 
    res$IMPACC= list(ARI_hc=adjustedRandIndex(css$labels,sc_label),
                   ARI_spec =adjustedRandIndex(sc_label,IMPACC_cluster(css$ConsensusMatrix,K,'spectral')),
                   cluster = css$labels,
                   time = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))),
                   feature_importance=css$feature_importance,
                   mat = css$ConsensusMatrix,
                   stop_point=css$nIter)
    remove(css)
  }
  if (MPACC == T){
    print('run MPACC')
    ## only adp observations 
    tic.clearlog()
    tic()    
    css= IMPACC(d=sc_scale,K = K,adaptiveFeature=FALSE, E=3,reps=nIter,pItem=pItem,pFeature=m,verbose = verbose)
    toc(log = TRUE)
    res$MPACC= list(ARI_hc=adjustedRandIndex(css$labels,sc_label),
                    ARI_spec =adjustedRandIndex(sc_label,IMPACC_cluster(css$ConsensusMatrix,K,'spectral')),
                    cluster = css$labels,
                    time = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))),
                    feature_importance=css$feature_importance,
                    mat = css$ConsensusMatrix,
                    stop_point=css$nIter)
    
    remove(css)
  }
  if (MPCC == T){
    # # random mini
    print('run MPCC')
    tic.clearlog()
    tic()    
    css = MPCC(d=sc_scale,K = K, h=0.95,reps=nIter,pItem=pItem,pFeature=m,early_stop=early_stop,verbose = verbose)
    toc(log = TRUE)
    

    res$MPCC=list(ARI_hc=adjustedRandIndex(css$labels,sc_label),
                  ARI_spec =adjustedRandIndex(sc_label,IMPACC_cluster(css$ConsensusMatrix,K,'spectral')),
                  cluster = css$labels,
                  time = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))),
                  mat = css$ConsensusMatrix,
                  stop_point=css$nIter)
    
  }
  if (consensus == T){
    print('run CSS')
    
    tic.clearlog()
    tic()
    
    css = ConsensusClusterPlus(as.matrix(sc_standardized),clusterAlg = 'hc',
                               maxK = K,pItem = pItem,
                               reps = res$MPCC$stop_point,innerLinkage="ward.D", finalLinkage="ward.D", distance="manhattan")
    ## reps = 1000
    toc(log = TRUE)
    
    res$consensus$ARI_hc = adjustedRandIndex(css[[K]]$consensusClass,sc_label)
    specc = IMPACC_cluster(css[[K]]$consensusMatrix,K,'spectral')
    res$consensus$ARI_spec = adjustedRandIndex(specc,sc_label)
    
    consensus_time = unlist(tic.log(format = TRUE))
    consensus_time = as.numeric(gsub("[^0-9.]", "",  consensus_time))
    tic.clearlog()
    res$consensus$time=consensus_time
    res$consensus$mat = css[[K]]$consensusMatrix
    remove(css)
    
  }
  
  ###### sparse kmean
  if (sparseKM == T){
    print('run sparseKM')
    
    tic()
    km.perm <- KMeansSparseCluster.permute(t(sc_scale),K=K,
                                           wbounds=seq(3,7,len=15),nperms=1)
    
    # run sparse k-means
    km.out <- KMeansSparseCluster(t(sc_scale),K=K,wbounds=km.perm$bestw)
    toc(log = TRUE)
    res[['sparseKM']]= list(cluster = km.out[[1]]$Cs,
                            ARI=adjustedRandIndex(km.out[[1]]$Cs,sc_label),
                            wcss = km.out[[1]]$wcss$wcss.perfeature,
                            bcss = km.out[[1]]$wcss$bcss.perfeature,
                            ws=km.out[[1]]$ws,
                            bestw = km.perm$bestw,
                            time = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE))))
    )
    tic.clearlog()
  }
  
  if (sparseHC == T){
    print('run sparseHC')
    
    tic()
    
    perm.out <- HierarchicalSparseCluster.permute(t(sc_scale), wbounds=c(1.5,2:6),
                                                  nperms=1)
    # Perform sparse hierarchical clustering
    sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists,
                                          wbound=perm.out$bestw, method="average")
    cl= cutree(sparsehc$hc,K)
    
    toc(log = TRUE)
    res[['sparseHclust']]=list(cluster = cl,
                               bestw=perm.out$bestw,
                               ARI=adjustedRandIndex(cl,sc_label),
                               ws = sparsehc$ws,
                               time = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE))))
    )
    tic.clearlog()
  }
  
  ## Seurat 
  if (seurat ==T){
    dd <- CreateSeuratObject(counts = as(as.matrix(sc), "sparseMatrix")  )
    all.genes <- rownames(dd)
    tic()
    dd <- FindVariableFeatures(dd, selection.method = "vst", nfeatures = 2000)
    dd <- ScaleData(dd, features = all.genes)
    
    dd <- RunPCA(dd, features = VariableFeatures(object = dd),npcs=min(ncol(sc_scale)-1,50))
    dd <- FindNeighbors(dd, dims = 1:10)
    toc(log = TRUE)
    resol = 0.1
    kk=0
    while (kk<K){
      dd2 <- FindClusters(dd, resolution = resol)
      kk = length(unique(dd2$seurat_clusters))
      resol=resol+0.05
    }

    tic()
    dd2 <- FindClusters(dd, resolution = resol-0.05)
    kk = length(unique(dd2$seurat_clusters))
    toc(log = TRUE)
    ser_clusters = dd2$seurat_clusters
    res[['Seurat']]=list(cluster = ser_clusters,
                         ARI=adjustedRandIndex(ser_clusters,sc_label),
                         time = sum(as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE))))))
    
    tic.clearlog()
  }
  ################
  ## sc3
  #####################
  
  if (sc3==T){
    tic()
    sce <- SingleCellExperiment(
      assays = list(
        counts = as.matrix(sc),
        logcounts = as.matrix(sc_scale)
      ), 
      colData = sc_label
    )
    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- sc3(sce, ks = K, biology = TRUE)
    toc(log=T)
    
    res[['SC3']]=list(cluster = sce[[2]],
                      ARI=adjustedRandIndex(sce[[2]],sc_label),
                      info = rowData(sce),
                      time = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))))
    
    tic.clearlog()
    
  }
  ################
  ## tsne
  #####################
  if (tsne == T){
  
      tic.clearlog()
      
      tic()
      tsne <- Rtsne(t(data.frame(sc_standardized)),
                    perplexity = min(floor((nrow(t(data.frame(sc_scale))) - 1) / 3),30),
                    )$Y
      toc(log = TRUE)

    time_tsne = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE))))
    tic.clearlog()
    
    tic()
    kmeans = kmeans((tsne),K)$cluster
    toc(log = TRUE)
    
    tic()
    hclust = cutree(hclust(dist((tsne),method = 'manhattan'),method='ward.D'),
                    K)
    toc(log = TRUE)
    tryCatch(
      {
      tic()
      spectral =kernlab::specc((tsne), centers=K,kernel = "rbfdot")
      toc(log = TRUE)
      tsne_spec = T
      },
     error =function(cond) {
        message("can't perform spectral")
        # Choose a return value in case of error
       tsne_spec  = F
     }
    )
    
    
    tic()
    kmedoid = cluster::pam((tsne), K)$clustering
    toc(log = TRUE)
    if (tsne_spec == T){
      cluster_Res = rbind(kmeans,hclust,spectral,kmedoid)
      ARIs = sapply(1:nrow(cluster_Res), function(i) adjustedRandIndex(cluster_Res[i,],sc_label))
      times = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE))))+time_tsne
      
      names(ARIs) =names(times) = c('tSNE+KMeans','tSNE+HC','tSNE+spectral','tSNE+KMedoid')
    }else{
      cluster_Res = rbind(kmeans,hclust,kmedoid)
      ARIs = sapply(1:nrow(cluster_Res), function(i) adjustedRandIndex(cluster_Res[i,],sc_label))
      times = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE))))+tsne_time
      names(ARIs) =names(times) = c('tSNE+KMeans','tSNE+HC','tSNE+KMedoid')
    }
    res[['tsne']] = list(cluster = cluster_Res,
                         ARI=ARIs,
                         time =times )
    tic.clearlog()
    
    
  }
  if (regular == T){
    ########## regular clusterings
    tic()
    kmeans = kmeans((t(sc_standardized)),K)$cluster
    toc(log = TRUE)
    
    tic()
    hclust = cutree(hclust(dist(t(sc_standardized),method = 'manhattan'),method='ward.D'),
                    K)
    toc(log = TRUE)
    
    tryCatch(
      {
        tic()
        spectral =kernlab::specc(t(sc_standardized), centers=K,kernel = "rbfdot")
        toc(log = TRUE)
        spec = T
      },
      error =function(cond) {
        message("can't perform spectral")
        # Choose a return value in case of error
        spec  = F
      }
    )
    tic()
    kmedoid = cluster::pam(t(sc_standardized), K)$clustering
    toc(log = TRUE)
    
    if (spec == T){
      cluster_Res = rbind(kmeans,hclust,spectral,kmedoid)
      ARIs = sapply(1:nrow(cluster_Res), function(i) adjustedRandIndex(cluster_Res[i,],sc_label))
      times = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE))))
      
      names(ARIs) =names(times) = c('KMeans','HC','spectral','KMedoid')
    }else{
      cluster_Res = rbind(kmeans,hclust,kmedoid)
      ARIs = sapply(1:nrow(cluster_Res), function(i) adjustedRandIndex(cluster_Res[i,],sc_label))
      times = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE))))
      names(ARIs) =names(times) = c('KMeans','HC','KMedoid')
    }
    res[['regular']] = list(cluster = cluster_Res,
                            ARI=ARIs,
                            time =times )
    tic.clearlog()
  }
  
  return(res)
}


