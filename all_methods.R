## function to run cluster algorithms in empirical study  
## regular: includes hclust, KMeans, KMedoid
source('IMPACC.r')
source('funmini.r')
library(ConsensusClusterPlus)
library(Rtsne)
library(SC3)
library(scater)
library(mclust)
library(Matrix)
library(Seurat)

all_methods = function(dat,nIter=300, m=0.1,n=0.25,IMPACC = T, MPACC =F,MPCC = T,
                       consensus =T,sparseKM = T,
                       sparseHC = T,regular =T,spec = T,CS3=F,seruat=F,DR_regular=F,verbose = F){
  sc_scale=dat$sc_cnt
  sc_label = dat$sc_label
  K = length(unique(sc_label))
  res=list()
  if (IMPACC == T){
    print('run IMPACC')
    tic.clearlog()
    tic()
    css= IMPACC(d=sc_scale,E=3,reps=nIter,pItem=0.25,pFeature=m,early_stop=T, verbose = verbose)
    toc(log = TRUE)
    res$IMPACC= list(ARI_hc=get_ari(css$consensus,'HClust',sc_label,K,'ward.D'),
                   ARI_spec = get_ari(css$consensus,'Spectral',sc_label,K,'ward.D'),
                   feature_importance=css$feature_importance,
                   time=as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))),
                   mat = css$consensus)
    remove(css)
  }
  if (MPCC == T){
    # # random mini
    print('run MPCC')
    tic.clearlog()
    tic()
    css = MPCC(d=sc_scale,h=0.95,reps=nIter,pItem=0.25,pFeature=m,verbose = verbose)
    toc(log = TRUE)
    stop_point=css$stop_point
    res$MPCC=list(ARI_hc=get_ari(css$consensus,'HClust',sc_label,K,'ward.D'),
                    ARI_spec = get_ari(css$consensus,'Spectral',sc_label,K,'ward.D'),
                    feature_importance=css$feature_importance,
                    time=as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))),
                    mat = css$consensus)
  }
  if (consensus == T){
    print('run CSS')
    
    tic.clearlog()
    tic()
    css = ConsensusClusterPlus(as.matrix(sc_scale),clusterAlg = 'hc',
                               maxK = K,pItem = 0.25,
                               reps = stop_point,innerLinkage="ward.D", finalLinkage="ward.D", distance="manhattan")
    ## reps = 1000
    toc(log = TRUE)
    res$consensus=list(ARI_hc=adjustedRandIndex(css[[K]]$consensusClass,sc_label),
                  ARI_spec = adjustedRandIndex(get_cluster(css[[K]]$consensusMatrix,K,'Spectral'),sc_label),
                  time=as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))),
                  mat = css[[K]]$consensusMatrix)
    remove(css)
  }
  ###### sparse KMean
  if (sparseKM == T){
    print('run sparseKM')
    tic.clearlog()
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
  if (regular == T){
    ########## regular clusterings
    tic()
    KMeans = kmeans(t(sc_scale),K)$cluster
    toc(log = TRUE)
    
    tic()
    hclust = cutree(hclust(dist(t(sc_scale),method = 'manhattan'),method='ward.D'),
                    K)
    toc(log = TRUE)
    
    if (spec == T){
      tic()
      spectral =kernlab::specc(t(sc_scale), centers=K,kernel = "rbfdot")
      toc(log = TRUE)
    }
    tic()
    KMedoid = cluster::pam(t(sc_scale), K)$clustering
    toc(log = TRUE)
    
    if (spec == T){
      cluster_Res = rbind(KMeans,hclust,spectral,KMedoid)
      ARIs = sapply(1:nrow(cluster_Res), function(i) adjustedRandIndex(cluster_Res[i,],sc_label))
      times = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE))))
      
      names(ARIs) =names(times) = c('KMeans','HClust','Spectral','KMedoid')
    }else{
      cluster_Res = rbind(KMeans,hclust,KMedoid)
      ARIs = sapply(1:nrow(cluster_Res), function(i) adjustedRandIndex(cluster_Res[i,],sc_label))
      times = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE))))
      names(ARIs) =names(times) = c('KMeans','HClust','KMedoid')
    }
    res[['regular']] = list(cluster = cluster_Res,
                            ARI=ARIs,
                            time =times )
    tic.clearlog()
  }
  
  if (DR_regular == T){
    tic.clearlog()
    
    tic()
    tsne <- Rtsne(t(data.frame(sc_scale)))$Y
    toc(log = TRUE)
    time_dr = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE))))
    tic.clearlog()
    
    tic()
    KMeans = kmeans((tsne),K)$cluster
    toc(log = TRUE)
    
    tic()
    hclust = cutree(hclust(dist((tsne),method = 'manhattan'),method='ward.D'),
                    K)
    toc(log = TRUE)
    
    tic()
    spectral =kernlab::specc((tsne), centers=K,kernel = "rbfdot")
    toc(log = TRUE)
    
    tic()
    KMedoid = cluster::pam((tsne), K)$clustering
    toc(log = TRUE)
    
    cluster_Res = rbind(KMeans,hclust,spectral,KMedoid)
    ARIs = sapply(1:nrow(cluster_Res), function(i) adjustedRandIndex(cluster_Res[i,],sc_label))
    times = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE))))+time_dr
    names(ARIs) =names(times) = c('tSNE+KMeans','tSNE+hclust','tSNE+spectral','tSNE+KMedoid')
  
  res[['DR_regular']] = list(cluster = cluster_Res,
                          ARI=ARIs,
                          time =times )
  tic.clearlog()
  }
  
  if (seruat == T){
    tic.clearlog()
    
    dd <- CreateSeuratObject(counts = as(sc_scale, "sparseMatrix")  )
    all.genes <- rownames(dd)
    tic()
    dd <- FindVariableFeatures(dd, selection.method = "vst", nfeatures = 2000)
    dd <- ScaleData(dd, features = all.genes)
    dd <- RunPCA(dd, features = VariableFeatures(object = dd))
    dd <- FindNeighbors(dd, dims = 1:10)
    dd <- FindClusters(dd, resolution = 0.3)
    toc(log = TRUE)
    ser_clusters = dd$seurat_clusters
    res[['Seurat']]=list(cluster = ser_clusters,
                          ARI=adjustedRandIndex(ser_clusters,sc_label),
                               time = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))))
                               
  }
  if (CS3 == T){
    sce <- SingleCellExperiment(
      assays = list(
        counts = as.matrix(sc_scale),
        
        logcounts = as.matrix(sc_scale)
        # logcounts = log2(as.matrix(sc) + 1)
      ), 
      colData = sc_label
    )
    rowData(sce)$feature_symbol <- rownames(sce)
    tic()
    sce <- sc3(sce, ks = K, biology = TRUE)
    toc(log = TRUE)
    res[['CS3']]=list(cluster = sce[[2]],
                         ARI=adjustedRandIndex(sce[[2]],sc_label),
                         time = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))))
    
  }
  return(res)
}
