## use rank ANOVA to measure importance 
library(matrixStats)
library(cluster) 
library(nnet)
library(Matrix)
library(mclust)
library(SC3)
library(Seurat)
library(Rtsne)
library(scater)
library(tictoc)
args = commandArgs(trailingOnly=TRUE)

scalematrix <- function(data) {
    cm <- rowMeans(data)
    csd <- rowSds(data, center = cm)
    csd[which(csd==0)]=0.001
    (data - cm) / csd
}
source('impacc.R')
source('impacc_variants.R')
dd = args[1]
print(dd) 
dat=readRDS(paste0('data/',dd,'.rds'))
    sc_scale = dat$sc_cnt
    sc_label = dat$sc_label
    K=length(unique(sc_label))
for (j in c(1:10)){    
	print('rank')    
    res = list()
    tic.clearlog()
    tic()    
    css= IMPACC_rankANOVA(d=sc_scale,K=K,reps=500)
    toc(log = TRUE)

    res$IMPACC_rankANOVA= list(ARI_hc=adjustedRandIndex(css$labels,sc_label),
             ARI_spec =adjustedRandIndex(sc_label,IMPACC_cluster(css$ConsensusMatrix,K,'spectral')),
             cluster = css$labels,
             time = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))),
             feature_importance=css$feature_importance,
             mat = css$ConsensusMatrix,
             stop_point=css$nIter)
        saveRDS(res,paste0('results/variants_',dd,'_',j,'.rds'))
print('mn')  


        tic.clearlog()
        tic()    
        css= IMPACC_multinomial(d=sc_scale,K=K,reps=500)
        toc(log = TRUE)
        
        res$IMPACC_multinomial= list(ARI_hc=adjustedRandIndex(css$labels,sc_label),
                              ARI_spec =adjustedRandIndex(sc_label,IMPACC_cluster(css$ConsensusMatrix,K,'spectral')),
                              cluster = css$labels,
                              time = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))),
                              feature_importance=css$feature_importance,
                              mat = css$ConsensusMatrix,
                              stop_point=css$nIter)
        saveRDS(res,paste0('results/variants_',dd,'_',j,'.rds'))

        tic.clearlog()
        tic()    
        css= IMPACC_silh(d=sc_scale,K=K,reps=500)
        toc(log = TRUE)
        
        res$IMPACC_silh= list(ARI_hc=adjustedRandIndex(css$labels,sc_label),
                            ARI_spec =adjustedRandIndex(sc_label,IMPACC_cluster(css$ConsensusMatrix,K,'spectral')),
                            cluster = css$labels,
                            time = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))),
                            feature_importance=css$feature_importance,
                            mat = css$ConsensusMatrix,
                            stop_point=css$nIter)
        saveRDS(res,paste0('results/variants_',dd,'_',j,'.rds'))
        print('oracle')  
        tic.clearlog()
        tic()    
        css= IMPACC_oracle(d=sc_scale,K=K,reps=500)
        toc(log = TRUE)
        
        res$IMPACC_oracle= list(ARI_hc=adjustedRandIndex(css$labels,sc_label),
                              ARI_spec =adjustedRandIndex(sc_label,IMPACC_cluster(css$ConsensusMatrix,K,'spectral')),
                              cluster = css$labels,
                              time = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))),
                              feature_importance=css$feature_importance,
                              mat = css$ConsensusMatrix,
                              stop_point=css$nIter)
        saveRDS(res,paste0('results/variants_',dd,'_',j,'.rds'))
        tic.clearlog()
        
    tic()
    css= IMPACC(d=sc_scale,K=K,reps=500)
    toc(log = TRUE)

    res$IMPACC= list(ARI_hc=adjustedRandIndex(css$labels,sc_label),
             ARI_spec =adjustedRandIndex(sc_label,IMPACC_cluster(css$ConsensusMatrix,K,'spectral')),
             cluster = css$labels,
             time = as.numeric(gsub("[^0-9.]", "",  unlist(tic.log(format = TRUE)))),
             feature_importance=css$feature_importance,
             mat = css$ConsensusMatrix,
             stop_point=css$nIter)
        saveRDS(res,paste0('results/variants_',dd,'_',j,'.rds'))	
}
