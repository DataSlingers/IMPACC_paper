source('impacc_record.R')
source('load_packages.R')

#for (j in c(1:10)){
for (dd in c('yan','goolam','biase','hiseq','darmanis','coil')){
 
	print(dd)
      dat=readRDS(paste0('data/',dd,'.rds'))
               
                
                sc_scale = as.matrix(dat$sc_cnt)
                sc_label = dat$sc_label
                K=length(unique(sc_label))

	    
        if (ncol(sc_scale)<100){
            n=0.5
        }else{
            n=0.25
        }	
	
	res = list()
    hs = c(0.7,0.8,0.9,0.95,0.99)

    for (hi in hs){
    print(hi)
        tic.clearlog()
        tic()    
        css= IMPACC_record(d=sc_scale,K =K, reps=200,h =hi,pItem = n,early_stop = F)
        toc(log = TRUE)


        labels = lapply(css$record, function(i) IMPACC_cluster(i,K))
        
        res$h[[as.character(hi)]] = sapply(labels, function(i) adjustedRandIndex(i,sc_label))
    }
        saveRDS(res,paste0('results/hyper_',dd,'.rds'))

    ps = c(0.01,0.05,0.1,0.15,0.2)
    
    for (pp in ps){
	  print(pp)
        tic.clearlog()
        tic()    
        css= IMPACC_record(d=sc_scale,K =K, reps=200,pp=pp,pItem = n,early_stop = F)
        toc(log = TRUE)
                labels = lapply(css$record, function(i) IMPACC_cluster(i,K))

        res$ps[[as.character(pp)]] = sapply(labels, function(i) adjustedRandIndex(i,sc_label))
 
    }
    saveRDS(res,paste0('results/hyper_',dd,'.rds'))
    alphas = seq(0.1,0.9,by=0.1)
    for (a in alphas){
        tic.clearlog()
        tic()    
                css= IMPACC_record(d=sc_scale,K =K, reps=200,alpha_F =a,pItem = n,early_stop = F)

        toc(log = TRUE)

	labels = lapply(css$record, function(i) IMPACC_cluster(i,K))

        res$alpha_f[[as.character(a)]] =  sapply(labels, function(i) adjustedRandIndex(i,sc_label))
    }
    saveRDS(res,paste0('results/hyper_',dd,'.rds'))
    alphas = seq(0.1,0.9,by=0.1)
    for (a in alphas){
        tic.clearlog()
        tic()    
        css= IMPACC_record(d=sc_scale,K =K, reps=200,alpha_I = a,pItem = n,early_stop = F)
        toc(log = TRUE)
        labels = lapply(css$record, function(i) IMPACC_cluster(i,K))
        res$alpha_i[[as.character(a)]] =  sapply(labels, function(i) adjustedRandIndex(i,sc_label))
    }
    saveRDS(res,paste0('results/hyper_',dd,'.rds'))
	qIs = c(0.7,0.8,0.9,0.95,0.99)
    for (qq in qIs){
    print(qq)
    	    tic.clearlog()
        tic()    
	        css= IMPACC_record(d=sc_scale,K =K, reps=200,qI = qq,pItem = n,early_stop = F)

        toc(log = TRUE)
        labels = lapply(css$record, function(i) IMPACC_cluster(i,K))
        res$qI[[as.character(qq)]] = sapply(labels, function(i) adjustedRandIndex(i,sc_label))
    }
    saveRDS(res,paste0('results/hyper_',dd,'.rds'))

    qFs = c(0.9,0.95,0.99,1,1.5,2)
    for (qq in qFs){
        tic.clearlog()
        tic()    
        css= IMPACC_record(d=sc_scale,K =K, reps=200,qF= qq,pItem = n,early_stop = F)
        toc(log = TRUE)
             labels = lapply(css$record, function(i) IMPACC_cluster(i,K))
        res$qF[[as.character(qq)]] = sapply(labels, function(i) adjustedRandIndex(i,sc_label))
    }
    saveRDS(res,paste0('results/hyper_',dd,'.rds'))

    models = c('EE+prob','EE','prob')
    for (m in models){
        print(m)
        tic.clearlog()
        tic()
        css= IMPACC_record(d=sc_scale,K =K, reps=200,sample=m,pItem = n,early_stop = F)
        toc(log = TRUE)
        labels = lapply(css$record, function(i) IMPACC_cluster(i,K))
        
        res$model[[m]] = sapply(labels, function(i) adjustedRandIndex(i,sc_label))
        saveRDS(res,paste0('results/hyper_',dd,'.rds'))
        
    }
    
}

