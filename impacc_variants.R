library(cluster) 
pv_anova_rank = function(X,y){
    X=as.matrix(X)
    n=length(y)
    n_y=levels(y)
    R = t(apply(X, 1, rank))
    TIES = apply(X, 1, table)
    TIEE = sapply(TIES, function(x) sum(x^3-x))
    STATISTIC = rowsums(sapply(n_y, function(w) rowsums(R[,y==w])^2/sum(y==w)))
    STATISTICS <- ((12 * STATISTIC / (n * (n + 1)) - 3 * (n + 1)) /
                      (1 - TIEE/ (n^3 - n)))
    PARAMETER <-  nlevels(y) - 1L
    PVAL <- sapply(STATISTICS, function(s) pchisq(s, PARAMETER, lower.tail = FALSE))
    return(PVAL)
}      
 


pv_multinom = function(X,y){
    multinom_score = function(x,y){
        fit=multinom(y ~ x)
        pred = predict(fit, x, "probs")
        if (nlevels(y)==2){
            pred = cbind(1-pred,pred)
            predd = pred[cbind(1:nrow(pred), y)]
        }else{
            predd = pred[cbind(1:nrow(pred), y)]
        }
        return(mean(predd))
    }
    return(sapply(c(1:nrow(X)), function(i) multinom_score(as.matrix(X)[i,],y)))
}


IMPACC_rankANOVA = function(d=NULL,
                            K=NULL,
                            adaptiveFeature = TRUE,
                            reps=100, ## number of iterations
                            pItem=0.25, ## minipatch size for observations
                            pFeature=0.1, ## minipatch size for features
                            innerLinkage="ward.D", ## internal HC linkage
                            distance="manhattan",## internal HC distance
                            h=0.95, ## cut internal tree at h quantile
                            E= 3, ## number of epochs in burn-in stage
                            qI=0.95, ## observations with greater than (qI) quantile of weights are classified to the set of high uncertain observations in adaptive sampling
                            qF=1, ## features with weights greater than (mean+qF*std) are classified to the set of high importance genes in adaptive sampling
                            alpha_I=0.5, ## learning rate for updating observation weights
                            alpha_F=0.5, ## learning rate for updating feature weights
                            pp=0.05, ## feature support threshold, a feature will be add to feature support if its p-value in ANOVA is smaller than pp
                            finalAlgorithm='hclust',
                            finalLinkage='ward.D',
                            early_stop=TRUE, ## whether perform early stop criteria
                            num_unchange = 5,
                            eps = 0.00001,
                            seed=NULL,
                            verbose=TRUE,
                            corUse="everything"){
    ###########
    ## The IMPACC function perform IMPACC (Interpretable MiniPatch Adaptive Consensus Clustering) with adaptive sampling scheme on both observations and features.
    ## It returns (1)consensus: final N by N consensus matrix; (2):feature_importance: feature importance scores
    ############
    # adaptiveFeature = TRUE
    # reps=100  ## number of iterations
    # pItem=0.25  ## minipatch size for observations
    # pFeature=0.1  ## minipatch size for features
    # innerLinkage="ward.D"  ## internal HC linkage
    # distance="manhattan" ## internal HC distance
    # h=0.95  ## cut internal tree at h quantile
    # E= 3  ## number of epochs in burn-in stage
    # qI=0.95  ## observations with greater than (qI) quantile of weights are classified to the set of high uncertain observations in adaptive sampling
    # qF=1  ## features with weights greater than (mean+qF*std) are classified to the set of high importance genes in adaptive sampling
    # alpha_I=0.5  ## learning rate for updating observation weights
    # alpha_F=0.5  ## learning rate for updating feature weights
    # pp=0.05  ## feature support threshold, a feature will be add to feature support if its p-value in ANOVA is smaller than pp
    # finalAlgorithm='hclust'
    # finalLinkage='ward.D'
    # early_stop=TRUE ## whether perform early stop criteria
    # seed=NULL
    # verbose=TRUE
    # corUse="everything"
    ######################
    ## check validity of data input
    check(d)
    d=data.frame(scalematrix(as.matrix(d)))
    ###############
    if(is.null(seed)==TRUE){
        seed=timeSeed = as.numeric(Sys.time())
    }
    set.seed(seed)
    
    
    update_IMPACC = function(){
        if(verbose){
            message(paste("  i =",i))
        }
        ### update observation weights
        confusion = rowMeans(CoAsso*(1-CoAsso))
        ww = (i+n_burnin)/subsample_o*confusion
        wi <<- alpha_I*wi+(1-alpha_I)*(ww/sum(ww))
        if (adaptiveFeature==TRUE){
            wi_p<<- alpha_F*wi_p+(1-alpha_F)*feature_score
        }else{
            wi_p<<-NULL
        }
        ## subsample
        sample_x = sample_MiniPatch(d, pItem, pFeature,pi_item=pi_item[i],pi_feature=pi_feature[i],
                                    weightsItem = wi, weightsFeature = wi_p,qI=qI,qF=qF)
        subsample_o<<- subsample_o+colnames(d)%in%colnames(d)[sample_x$subcols]
        subsample_f<<- subsample_f+rownames(d)%in%rownames(d)[sample_x$subrows]
        ############################
        ## clustering
        ########################
        this_assignment = cluster_algo(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage)
        
        #######################################
        ######### for feature importance
        #######################################
        if (adaptiveFeature==TRUE){
            
            pvalu_rank = pv_anova_rank(as.matrix(sample_x$submat),as.factor(this_assignment))
            pv_rank = quantile(na.omit(pvalu_rank),pp)
            pvalue = pvalu_rank<=pv_rank
            
            feature_support[rownames(sample_x$submat)[which(pvalue)]]<<- 1+ feature_support[rownames(sample_x$submat)[which(pvalue)]]
            feature_score<<- feature_support/subsample_f
        }else{
            feature_support = feature_score <<- NULL
        }
        #######################################
        
        ##########################################
        ##mCount stores number of times a sample pair was sampled together.
        mCount<<- connectivityMatrix(rep(1,length(sample_x$subcols)),
                                     mCount,
                                     sample_x$subcols)
        
        ##ml stores number of times a sample pair was sampled together.
        ml<<- connectivityMatrix(this_assignment,
                                 ml,
                                 sample_x$subcols)
        
        this_coA = ml / mCount
        this_coA[mCount==0]=0
        CoAsso  <<- this_coA
        
    }
    ## initialize parameters
    n=ncol(d)
    if (adaptiveFeature==TRUE){
        pi_item= pi_feature <- seq(0.5,1,by=0.5/reps)
    }else{
        pi_item <- seq(0.5,1,by=0.5/reps)
        pi_feature <- rep(1,reps)
    }
    
    ############################
    ## calculate # of iterations needed for burn-in stage
    ############################
    if (100%%(pItem*100)==0){
        p_sample = rep(pItem,1/pItem)
    }else{
        p_sample = c(rep(pItem,floor(1/pItem)),1%%pItem)
    }
    if (100%%(pFeature*100)==0){
        p_sample_row = rep(pFeature,1/pFeature)
    }else{
        p_sample_row = c(rep(pFeature,floor(1/pFeature)),1%%pFeature)
    }
    lcm = max(length(p_sample_row),length(p_sample))
    num_partition = E*lcm ## number of partition needed for burn-in stage
    conf_record = array(0,dim = c(reps+num_partition,n))
    
    ######################################
    ## BURN IN STAGE
    ##############################################
    message('Burn-in stage')
    sample_Burn = sampleBurin(d,num_partition,p_sample,p_sample_row) ## 'sampleBurin' is in funmini.R
    n_burnin=length(sample_Burn)
    subsample_o=rep(0,ncol(d))
    subsample_f=rep(0,nrow(d))
    if (adaptiveFeature==TRUE){
        feature_support = rep(0,nrow(d)) ## feature support
        names(feature_support)= rownames(d)
    }
    
    CoAsso=mCount=ml=matrix(c(0),ncol=n,nrow=n)
    
    for (i in 1:length(sample_Burn)){
        if(verbose){
            message(paste("  i =",i))
        }
        
        sample_x = sample_Burn[[i]]
        this_assignment = cluster_algo(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage)
        ## update records on minipatch sampling
        subsample_o = subsample_o+colnames(d)%in%colnames(d)[sample_x$subcols]
        subsample_f = subsample_f+rownames(d)%in%rownames(d)[sample_x$subrows]
        if (adaptiveFeature==TRUE){
            ##########################################
            ######## for feature importance
            ###########################################
            pvalu_rank = pv_anova_rank(as.matrix(sample_x$submat),as.factor(this_assignment))
            pv_rank = quantile(na.omit(pvalu_rank),pp)
            pvalue = pvalu_rank<=pv_rank
            
            feature_support[rownames(sample_x$submat)[which(pvalue)]]=1+ feature_support[rownames(sample_x$submat)[which(pvalue)]]
            feature_score=feature_support/(subsample_f)
            feature_score[subsample_f==0]=0
        }else{
            feature_support = feature_score<-NULL
        }
        ###################################
        mCount <- connectivityMatrix(rep(1,length(sample_x$subcols)),
                                     mCount,
                                     sample_x$subcols)
        
        ml <- connectivityMatrix(this_assignment,
                                 ml,
                                 sample_x$subcols)
        CoAsso = ml / mCount
        CoAsso[mCount==0]=0
        
    }
    
    wi=rep(1/n,n)
    wi_p=rep(1/nrow(d),nrow(d))
    
    # records = list('CoAsso'=CoAsso, 'mCount'=mCount,'ml'=ml,
    #                'wi'=wi,'wi_p'=wi_p,'subsample_o'=subsample_o,'subsample_f'=subsample_f,
    #                'feature_support'=feature_support,'feature_score'=feature_score)
    
    ############################################
    ## adaptive stage
    #####################################################
    ##message("adaptive stage")
    
    if (early_stop == F){
        message('No early stopping')
        for (i in seq(1,reps)){
            update_IMPACC()
        }
    }else{
        conf_q=NULL
        continue=T
        i=1
        while (continue==T & i<reps){
            update_IMPACC()
            
            conf_q=c(conf_q,quantile(rowMeans(CoAsso*(1-CoAsso)),0.9)) ## take 90% quantile of the confusion
            if (i>num_unchange){
                ## find change of quantile
                cm = sapply(c(i:(i-(num_unchange-1))), function(x) abs(conf_q[x]-conf_q[x-1]))
                if (max(cm)<eps){
                    continue=F
                    message(paste0('Stop at iteration ',i+num_partition))
                    
                }else{
                    continue=T
                }
            }
            i=i+1
        }
    }
    
    #heatmap(CoAsso)
    labels = IMPACC_cluster(CoAsso,K=K,finalAlgorithm=finalAlgorithm,finalLinkage=finalLinkage)
    
    return(list(ConsensusMatrix = CoAsso,
                labels=as.character(labels),
                feature_importance=feature_score,
                nIter = i))
}


#############################
## multinomial regression 
############################
IMPACC_multinomial = function(d=NULL,
                              K=NULL,
                              adaptiveFeature = TRUE,
                              reps=100, ## number of iterations
                              pItem=0.25, ## minipatch size for observations
                              pFeature=0.1, ## minipatch size for features
                              innerLinkage="ward.D", ## internal HC linkage
                              distance="manhattan",## internal HC distance
                              h=0.95, ## cut internal tree at h quantile
                              E= 3, ## number of epochs in burn-in stage
                              qI=0.95, ## observations with greater than (qI) quantile of weights are classified to the set of high uncertain observations in adaptive sampling
                              qF=1, ## features with weights greater than (mean+qF*std) are classified to the set of high importance genes in adaptive sampling
                              alpha_I=0.5, ## learning rate for updating observation weights
                              alpha_F=0.5, ## learning rate for updating feature weights
                              pp=0.05, ## feature support threshold, a feature will be add to feature support if its p-value in ANOVA is smaller than pp
                              finalAlgorithm='hclust',
                              finalLinkage='ward.D',
                              early_stop=TRUE, ## whether perform early stop criteria
                              num_unchange = 5,
                              eps = 0.00001,
                              seed=NULL,
                              verbose=TRUE,
                              corUse="everything"){
    ###########
    ## The IMPACC function perform IMPACC (Interpretable MiniPatch Adaptive Consensus Clustering) with adaptive sampling scheme on both observations and features.
    ## It returns (1)consensus: final N by N consensus matrix; (2):feature_importance: feature importance scores
    
    ######################
    ## check validity of data input
    check(d)
    d=data.frame(scalematrix(as.matrix(d)))
    ###############
    if(is.null(seed)==TRUE){
        seed=timeSeed = as.numeric(Sys.time())
    }
    set.seed(seed)
    
    
    update_IMPACC = function(){
        if(verbose){
            message(paste("  i =",i))
        }
        ### update observation weights
        confusion = rowMeans(CoAsso*(1-CoAsso))
        ww = (i+n_burnin)/subsample_o*confusion
        wi <<- alpha_I*wi+(1-alpha_I)*(ww/sum(ww))
        if (adaptiveFeature==TRUE){
            wi_p<<- alpha_F*wi_p+(1-alpha_F)*feature_score
        }else{
            wi_p<<-NULL
        }
        ## subsample
        sample_x = sample_MiniPatch(d, pItem, pFeature,pi_item=pi_item[i],pi_feature=pi_feature[i],
                                    weightsItem = wi, weightsFeature = wi_p,qI=qI,qF=qF)
        subsample_o<<- subsample_o+colnames(d)%in%colnames(d)[sample_x$subcols]
        subsample_f<<- subsample_f+rownames(d)%in%rownames(d)[sample_x$subrows]
        ############################
        ## clustering
        ########################
        this_assignment = cluster_algo(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage)
        
        #######################################
        ######### for feature importance
        #######################################
        if (adaptiveFeature==TRUE){
	prob_all = pv_multinom(as.matrix(sample_x$submat),as.factor(this_assignment))
            cutoff = quantile(na.omit(prob_all),1-pp)
            pvalue = prob_all>=cutoff
            
            feature_support[rownames(sample_x$submat)[which(pvalue)]]<<- 1+ feature_support[rownames(sample_x$submat)[which(pvalue)]]
            feature_score<<- feature_support/subsample_f
        }else{
            feature_support = feature_score <<- NULL
        }
        #######################################
        
        ##########################################
        ##mCount stores number of times a sample pair was sampled together.
        mCount<<- connectivityMatrix(rep(1,length(sample_x$subcols)),
                                     mCount,
                                     sample_x$subcols)
        
        ##ml stores number of times a sample pair was sampled together.
        ml<<- connectivityMatrix(this_assignment,
                                 ml,
                                 sample_x$subcols)
        
        this_coA = ml / mCount
        this_coA[mCount==0]=0
        CoAsso  <<- this_coA
        
    }
    ## initialize parameters
    n=ncol(d)
    if (adaptiveFeature==TRUE){
        pi_item= pi_feature <- seq(0.5,1,by=0.5/reps)
    }else{
        pi_item <- seq(0.5,1,by=0.5/reps)
        pi_feature <- rep(1,reps)
    }
    
    ############################
    ## calculate # of iterations needed for burn-in stage
    ############################
    if (100%%(pItem*100)==0){
        p_sample = rep(pItem,1/pItem)
    }else{
        p_sample = c(rep(pItem,floor(1/pItem)),1%%pItem)
    }
    if (100%%(pFeature*100)==0){
        p_sample_row = rep(pFeature,1/pFeature)
    }else{
        p_sample_row = c(rep(pFeature,floor(1/pFeature)),1%%pFeature)
    }
    lcm = max(length(p_sample_row),length(p_sample))
    num_partition = E*lcm ## number of partition needed for burn-in stage
    conf_record = array(0,dim = c(reps+num_partition,n))
    
    ######################################
    ## BURN IN STAGE
    ##############################################
    message('Burn-in stage')
    sample_Burn = sampleBurin(d,num_partition,p_sample,p_sample_row) ## 'sampleBurin' is in funmini.R
    n_burnin=length(sample_Burn)
    subsample_o=rep(0,ncol(d))
    subsample_f=rep(0,nrow(d))
    if (adaptiveFeature==TRUE){
        feature_support = rep(0,nrow(d)) ## feature support
        names(feature_support)= rownames(d)
    }
    
    CoAsso=mCount=ml=matrix(c(0),ncol=n,nrow=n)
    
    for (i in 1:length(sample_Burn)){
        if(verbose){
            message(paste("  i =",i))
        }
        
        sample_x = sample_Burn[[i]]
        this_assignment = cluster_algo(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage)
        ## update records on minipatch sampling
        subsample_o = subsample_o+colnames(d)%in%colnames(d)[sample_x$subcols]
        subsample_f = subsample_f+rownames(d)%in%rownames(d)[sample_x$subrows]
        if (adaptiveFeature==TRUE){
            ##########################################
            ######## for feature importance
            ###########################################
        prob_all = pv_multinom(as.matrix(sample_x$submat),as.factor(this_assignment))

            cutoff = quantile(na.omit(prob_all),1-pp)
            pvalue = prob_all>=cutoff
            
            feature_support[rownames(sample_x$submat)[which(pvalue)]]=1+ feature_support[rownames(sample_x$submat)[which(pvalue)]]
            feature_score=feature_support/(subsample_f)
            feature_score[subsample_f==0]=0
        }else{
            feature_support = feature_score<-NULL
        }
        ###################################
        mCount <- connectivityMatrix(rep(1,length(sample_x$subcols)),
                                     mCount,
                                     sample_x$subcols)
        
        ml <- connectivityMatrix(this_assignment,
                                 ml,
                                 sample_x$subcols)
        CoAsso = ml / mCount
        CoAsso[mCount==0]=0
        
    }
    
    wi=rep(1/n,n)
    wi_p=rep(1/nrow(d),nrow(d))
    
    # records = list('CoAsso'=CoAsso, 'mCount'=mCount,'ml'=ml,
    #                'wi'=wi,'wi_p'=wi_p,'subsample_o'=subsample_o,'subsample_f'=subsample_f,
    #                'feature_support'=feature_support,'feature_score'=feature_score)
    
    ############################################
    ## adaptive stage
    #####################################################
    ##message("adaptive stage")
    
    if (early_stop == F){
        message('No early stopping')
        for (i in seq(1,reps)){
            update_IMPACC()
        }
    }else{
        conf_q=NULL
        continue=T
        i=1
        while (continue==T & i<reps){
            update_IMPACC()
            
            conf_q=c(conf_q,quantile(rowMeans(CoAsso*(1-CoAsso)),0.9)) ## take 90% quantile of the confusion
            if (i>num_unchange){
                ## find change of quantile
                cm = sapply(c(i:(i-(num_unchange-1))), function(x) abs(conf_q[x]-conf_q[x-1]))
                if (max(cm)<eps){
                    continue=F
                    message(paste0('Stop at iteration ',i+num_partition))
                    
                }else{
                    continue=T
                }
            }
            i=i+1
        }
    }
    
    #heatmap(CoAsso)
    labels = IMPACC_cluster(CoAsso,K=K,finalAlgorithm=finalAlgorithm,finalLinkage=finalLinkage)
    
    return(list(ConsensusMatrix = CoAsso,
                labels=as.character(labels),
                feature_importance=feature_score,
                nIter = i))
}


## cut tree by siloutte score 
IMPACC_silh = function(d=NULL,
                       K=NULL,
                       adaptiveFeature = TRUE,
                       reps=100, ## number of iterations
                       pItem=0.25, ## minipatch size for observations
                       pFeature=0.1, ## minipatch size for features
                       innerLinkage="ward.D", ## internal HC linkage
                       distance="manhattan",## internal HC distance
                       h=0.95, ## cut internal tree at h quantile
                       E= 3, ## number of epochs in burn-in stage
                       qI=0.95, ## observations with greater than (qI) quantile of weights are classified to the set of high uncertain observations in adaptive sampling
                       qF=1, ## features with weights greater than (mean+qF*std) are classified to the set of high importance genes in adaptive sampling
                       alpha_I=0.5, ## learning rate for updating observation weights
                       alpha_F=0.5, ## learning rate for updating feature weights
                       pp=0.05, ## feature support threshold, a feature will be add to feature support if its p-value in ANOVA is smaller than pp
                       finalAlgorithm='hclust',
                       finalLinkage='ward.D',
                       early_stop=TRUE, ## whether perform early stop criteria
                       num_unchange = 5,
                       eps = 0.00001,
                       seed=NULL,
                       verbose=TRUE,
                       corUse="everything"){
    
    ######################
    ## check validity of data input
    check(d)
    d=data.frame(scalematrix(as.matrix(d)))
    ###############
    if(is.null(seed)==TRUE){
        seed=timeSeed = as.numeric(Sys.time())
    }
    set.seed(seed)
    
    
    update_IMPACC = function(){
        if(verbose){
            message(paste("  i =",i))
        }
        ### update observation weights
        confusion = rowMeans(CoAsso*(1-CoAsso))
        ww = (i+n_burnin)/subsample_o*confusion
        wi <<- alpha_I*wi+(1-alpha_I)*(ww/sum(ww))
        if (adaptiveFeature==TRUE){
            wi_p<<- alpha_F*wi_p+(1-alpha_F)*feature_score
        }else{
            wi_p<<-NULL
        }
        ## subsample
        sample_x = sample_MiniPatch(d, pItem, pFeature,pi_item=pi_item[i],pi_feature=pi_feature[i],
                                    weightsItem = wi, weightsFeature = wi_p,qI=qI,qF=qF)
        subsample_o<<- subsample_o+colnames(d)%in%colnames(d)[sample_x$subcols]
        subsample_f<<- subsample_f+rownames(d)%in%rownames(d)[sample_x$subrows]
        ############################
        ## clustering
        ########################
        this_assignment = cluster_algo_silh(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage)
        
        #######################################
        ######### for feature importance
        #######################################
        if (adaptiveFeature==TRUE){
            pvalu = pv_anova(as.matrix(sample_x$submat),as.factor(this_assignment))
            pv = quantile(na.omit(pvalu),pp)
            pvalue = pvalu<=pv
            feature_support[rownames(sample_x$submat)[which(pvalue)]]<<- 1+ feature_support[rownames(sample_x$submat)[which(pvalue)]]
            feature_score<<- feature_support/subsample_f
        }else{
            feature_support = feature_score <<- NULL
        }
        #######################################
        
        ##########################################
        ##mCount stores number of times a sample pair was sampled together.
        mCount<<- connectivityMatrix(rep(1,length(sample_x$subcols)),
                                     mCount,
                                     sample_x$subcols)
        
        ##ml stores number of times a sample pair was sampled together.
        ml<<- connectivityMatrix(this_assignment,
                                 ml,
                                 sample_x$subcols)
        
        this_coA = ml / mCount
        this_coA[mCount==0]=0
        CoAsso  <<- this_coA
        
    }
    ## initialize parameters
    n=ncol(d)
    if (adaptiveFeature==TRUE){
        pi_item= pi_feature <- seq(0.5,1,by=0.5/reps)
    }else{
        pi_item <- seq(0.5,1,by=0.5/reps)
        pi_feature <- rep(1,reps)
    }
    
    ############################
    ## calculate # of iterations needed for burn-in stage
    ############################
    if (100%%(pItem*100)==0){
        p_sample = rep(pItem,1/pItem)
    }else{
        p_sample = c(rep(pItem,floor(1/pItem)),1%%pItem)
    }
    if (100%%(pFeature*100)==0){
        p_sample_row = rep(pFeature,1/pFeature)
    }else{
        p_sample_row = c(rep(pFeature,floor(1/pFeature)),1%%pFeature)
    }
    lcm = max(length(p_sample_row),length(p_sample))
    num_partition = E*lcm ## number of partition needed for burn-in stage
    conf_record = array(0,dim = c(reps+num_partition,n))
    
    ######################################
    ## BURN IN STAGE
    ##############################################
    message('Burn-in stage')
    sample_Burn = sampleBurin(d,num_partition,p_sample,p_sample_row) ## 'sampleBurin' is in funmini.R
    n_burnin=length(sample_Burn)
    subsample_o=rep(0,ncol(d))
    subsample_f=rep(0,nrow(d))
    if (adaptiveFeature==TRUE){
        feature_support = rep(0,nrow(d)) ## feature support
        names(feature_support)= rownames(d)
    }
    
    CoAsso=mCount=ml=matrix(c(0),ncol=n,nrow=n)
    
    for (i in 1:length(sample_Burn)){
        if(verbose){
            message(paste("  i =",i))
        }
        
        sample_x = sample_Burn[[i]]
        this_assignment = cluster_algo_silh(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage)
        ## update records on minipatch sampling
        subsample_o = subsample_o+colnames(d)%in%colnames(d)[sample_x$subcols]
        subsample_f = subsample_f+rownames(d)%in%rownames(d)[sample_x$subrows]
        if (adaptiveFeature==TRUE){
            ##########################################
            ######## for feature importance
            ###########################################
            pvalu = pv_anova(as.matrix(sample_x$submat),as.factor(this_assignment))
            pv = quantile(na.omit(pvalu),pp)
            pvalue = pvalu<=pv
            feature_support[rownames(sample_x$submat)[which(pvalue)]]=1+ feature_support[rownames(sample_x$submat)[which(pvalue)]]
            feature_score=feature_support/(subsample_f)
            feature_score[subsample_f==0]=0
        }else{
            feature_support = feature_score<-NULL
        }
        ###################################
        mCount <- connectivityMatrix(rep(1,length(sample_x$subcols)),
                                     mCount,
                                     sample_x$subcols)
        
        ml <- connectivityMatrix(this_assignment,
                                 ml,
                                 sample_x$subcols)
        CoAsso = ml / mCount
        CoAsso[mCount==0]=0
        
    }
    
    wi=rep(1/n,n)
    wi_p=rep(1/nrow(d),nrow(d))
    
    # records = list('CoAsso'=CoAsso, 'mCount'=mCount,'ml'=ml,
    #                'wi'=wi,'wi_p'=wi_p,'subsample_o'=subsample_o,'subsample_f'=subsample_f,
    #                'feature_support'=feature_support,'feature_score'=feature_score)
    
    ############################################
    ## adaptive stage
    #####################################################
    ##message("adaptive stage")
    
    if (early_stop == F){
        message('No early stopping')
        for (i in seq(1,reps)){
            update_IMPACC()
        }
    }else{
        conf_q=NULL
        continue=T
        i=1
        while (continue==T & i<reps){
            update_IMPACC()
            
            conf_q=c(conf_q,quantile(rowMeans(CoAsso*(1-CoAsso)),0.9)) ## take 90% quantile of the confusion
            if (i>num_unchange){
                ## find change of quantile
                cm = sapply(c(i:(i-(num_unchange-1))), function(x) abs(conf_q[x]-conf_q[x-1]))
                if (max(cm)<eps){
                    continue=F
                    message(paste0('Stop at iteration ',i+num_partition))
                    
                }else{
                    continue=T
                }
            }
            i=i+1
        }
    }
    
    #heatmap(CoAsso)
    labels = IMPACC_cluster(CoAsso,K=K,finalAlgorithm=finalAlgorithm,finalLinkage=finalLinkage)
    
    return(list(ConsensusMatrix = CoAsso,
                labels=as.character(labels),
                feature_importance=feature_score,
                nIter = i))
}


### SPECIFY CLUSTERING MODELS
cluster_algo_silh = function(submat,h,distance,innerLinkage){
    acceptable.distance <- c( "euclidean", "maximum", "manhattan", "canberra", "binary","minkowski",
                              "pearson", "spearman" )
    ## check distance function
    if ( inherits( distance, "character" ) ) {
        if ( ! distance %in%  acceptable.distance  &  ( class(try(get(distance),silent=T))!="function") ){
            stop("unsupported distance.")}
        
        
        if(distance=="pearson" | distance=="spearman"){
            dis <- as.dist( 1-cor(submat,method=distance,use=corUse ))
        }else if( class(try(get(distance),silent=T))=="function"){
            dis <- get(distance)(t(submat))
        }else{
            dis = dist(t(submat),method = distance)
        }
    }else{stop("unsupported distance specified.")}
    
    
    this_cluster = hclust(dis, method=innerLinkage)
    silh=rep(0,9)
    for (k in c(2:10)){
        ct = cutree(this_cluster,k)
        ss <- silhouette(ct,dis)
        silh[k-1] = mean(ss[,3])
    }
    kk =c(2:10)[which.max(silh)]
    this_assignment = cutree(this_cluster,kk)
    return(this_assignment=this_assignment)
}
cluster_algo_oracle = function(submat,h,distance,innerLinkage,K){
    acceptable.distance <- c( "euclidean", "maximum", "manhattan", "canberra", "binary","minkowski",
                              "pearson", "spearman" )
    ## check distance function
    if ( inherits( distance, "character" ) ) {
        if ( ! distance %in%  acceptable.distance  &  ( class(try(get(distance),silent=T))!="function") ){
            stop("unsupported distance.")}
        
        
        if(distance=="pearson" | distance=="spearman"){
            dis <- as.dist( 1-cor(submat,method=distance,use=corUse ))
        }else if( class(try(get(distance),silent=T))=="function"){
            dis <- get(distance)(t(submat))
        }else{
            dis = dist(t(submat),method = distance)
        }
    }else{stop("unsupported distance specified.")}
    
    
    this_cluster = hclust(dis, method=innerLinkage)
    this_assignment = cutree(this_cluster,K)
    return(this_assignment=this_assignment)
}

IMPACC_oracle = function(d=NULL,
                         K=NULL,
                         adaptiveFeature = TRUE,
                         reps=100, ## number of iterations
                         pItem=0.25, ## minipatch size for observations
                         pFeature=0.1, ## minipatch size for features
                         innerLinkage="ward.D", ## internal HC linkage
                         distance="manhattan",## internal HC distance
                         h=0.95, ## cut internal tree at h quantile
                         E= 3, ## number of epochs in burn-in stage
                         qI=0.95, ## observations with greater than (qI) quantile of weights are classified to the set of high uncertain observations in adaptive sampling
                         qF=1, ## features with weights greater than (mean+qF*std) are classified to the set of high importance genes in adaptive sampling
                         alpha_I=0.5, ## learning rate for updating observation weights
                         alpha_F=0.5, ## learning rate for updating feature weights
                         pp=0.05, ## feature support threshold, a feature will be add to feature support if its p-value in ANOVA is smaller than pp
                         finalAlgorithm='hclust',
                         finalLinkage='ward.D',
                         early_stop=TRUE, ## whether perform early stop criteria
                         num_unchange = 5,
                         eps = 0.00001,
                         seed=NULL,
                         verbose=TRUE,
                         corUse="everything"){
    
    ######################
    ## check validity of data input
    check(d)
    d=data.frame(scalematrix(as.matrix(d)))
    ###############
    if(is.null(seed)==TRUE){
        seed=timeSeed = as.numeric(Sys.time())
    }
    set.seed(seed)
    
    
    update_IMPACC = function(){
        if(verbose){
            message(paste("  i =",i))
        }
        ### update observation weights
        confusion = rowMeans(CoAsso*(1-CoAsso))
        ww = (i+n_burnin)/subsample_o*confusion
        wi <<- alpha_I*wi+(1-alpha_I)*(ww/sum(ww))
        if (adaptiveFeature==TRUE){
            wi_p<<- alpha_F*wi_p+(1-alpha_F)*feature_score
        }else{
            wi_p<<-NULL
        }
        ## subsample
        sample_x = sample_MiniPatch(d, pItem, pFeature,pi_item=pi_item[i],pi_feature=pi_feature[i],
                                    weightsItem = wi, weightsFeature = wi_p,qI=qI,qF=qF)
        subsample_o<<- subsample_o+colnames(d)%in%colnames(d)[sample_x$subcols]
        subsample_f<<- subsample_f+rownames(d)%in%rownames(d)[sample_x$subrows]
        ############################
        ## clustering
        ########################
        this_assignment = cluster_algo_oracle(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage,K=K)
        
        #######################################
        ######### for feature importance
        #######################################
        if (adaptiveFeature==TRUE){
            pvalu = pv_anova(as.matrix(sample_x$submat),as.factor(this_assignment))
            pv = quantile(na.omit(pvalu),pp)
            pvalue = pvalu<=pv
            feature_support[rownames(sample_x$submat)[which(pvalue)]]<<- 1+ feature_support[rownames(sample_x$submat)[which(pvalue)]]
            feature_score<<- feature_support/subsample_f
        }else{
            feature_support = feature_score <<- NULL
        }
        #######################################
        
        ##########################################
        ##mCount stores number of times a sample pair was sampled together.
        mCount<<- connectivityMatrix(rep(1,length(sample_x$subcols)),
                                     mCount,
                                     sample_x$subcols)
        
        ##ml stores number of times a sample pair was sampled together.
        ml<<- connectivityMatrix(this_assignment,
                                 ml,
                                 sample_x$subcols)
        
        this_coA = ml / mCount
        this_coA[mCount==0]=0
        CoAsso  <<- this_coA
        
    }
    ## initialize parameters
    n=ncol(d)
    if (adaptiveFeature==TRUE){
        pi_item= pi_feature <- seq(0.5,1,by=0.5/reps)
    }else{
        pi_item <- seq(0.5,1,by=0.5/reps)
        pi_feature <- rep(1,reps)
    }
    
    ############################
    ## calculate # of iterations needed for burn-in stage
    ############################
    if (100%%(pItem*100)==0){
        p_sample = rep(pItem,1/pItem)
    }else{
        p_sample = c(rep(pItem,floor(1/pItem)),1%%pItem)
    }
    if (100%%(pFeature*100)==0){
        p_sample_row = rep(pFeature,1/pFeature)
    }else{
        p_sample_row = c(rep(pFeature,floor(1/pFeature)),1%%pFeature)
    }
    lcm = max(length(p_sample_row),length(p_sample))
    num_partition = E*lcm ## number of partition needed for burn-in stage
    conf_record = array(0,dim = c(reps+num_partition,n))
    
    ######################################
    ## BURN IN STAGE
    ##############################################
    message('Burn-in stage')
    sample_Burn = sampleBurin(d,num_partition,p_sample,p_sample_row) ## 'sampleBurin' is in funmini.R
    n_burnin=length(sample_Burn)
    subsample_o=rep(0,ncol(d))
    subsample_f=rep(0,nrow(d))
    if (adaptiveFeature==TRUE){
        feature_support = rep(0,nrow(d)) ## feature support
        names(feature_support)= rownames(d)
    }
    
    CoAsso=mCount=ml=matrix(c(0),ncol=n,nrow=n)
    
    for (i in 1:length(sample_Burn)){
        if(verbose){
            message(paste("  i =",i))
        }
        
        sample_x = sample_Burn[[i]]
        this_assignment = cluster_algo_oracle(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage,K=K)
        ## update records on minipatch sampling
        subsample_o = subsample_o+colnames(d)%in%colnames(d)[sample_x$subcols]
        subsample_f = subsample_f+rownames(d)%in%rownames(d)[sample_x$subrows]
        if (adaptiveFeature==TRUE){
            ##########################################
            ######## for feature importance
            ###########################################
            pvalu = pv_anova(as.matrix(sample_x$submat),as.factor(this_assignment))
            pv = quantile(na.omit(pvalu),pp)
            pvalue = pvalu<=pv
            feature_support[rownames(sample_x$submat)[which(pvalue)]]=1+ feature_support[rownames(sample_x$submat)[which(pvalue)]]
            feature_score=feature_support/(subsample_f)
            feature_score[subsample_f==0]=0
        }else{
            feature_support = feature_score<-NULL
        }
        ###################################
        mCount <- connectivityMatrix(rep(1,length(sample_x$subcols)),
                                     mCount,
                                     sample_x$subcols)
        
        ml <- connectivityMatrix(this_assignment,
                                 ml,
                                 sample_x$subcols)
        CoAsso = ml / mCount
        CoAsso[mCount==0]=0
        
    }
    
    wi=rep(1/n,n)
    wi_p=rep(1/nrow(d),nrow(d))
    
    # records = list('CoAsso'=CoAsso, 'mCount'=mCount,'ml'=ml,
    #                'wi'=wi,'wi_p'=wi_p,'subsample_o'=subsample_o,'subsample_f'=subsample_f,
    #                'feature_support'=feature_support,'feature_score'=feature_score)
    
    ############################################
    ## adaptive stage
    #####################################################
    ##message("adaptive stage")
    
    if (early_stop == F){
        message('No early stopping')
        for (i in seq(1,reps)){
            update_IMPACC()
        }
    }else{
        conf_q=NULL
        continue=T
        i=1
        while (continue==T & i<reps){
            update_IMPACC()
            
            conf_q=c(conf_q,quantile(rowMeans(CoAsso*(1-CoAsso)),0.9)) ## take 90% quantile of the confusion
            if (i>num_unchange){
                ## find change of quantile
                cm = sapply(c(i:(i-(num_unchange-1))), function(x) abs(conf_q[x]-conf_q[x-1]))
                if (max(cm)<eps){
                    continue=F
                    message(paste0('Stop at iteration ',i+num_partition))
                    
                }else{
                    continue=T
                }
            }
            i=i+1
        }
    }
    
    #heatmap(CoAsso)
    labels = IMPACC_cluster(CoAsso,K=K,finalAlgorithm=finalAlgorithm,finalLinkage=finalLinkage)
    
    return(list(ConsensusMatrix = CoAsso,
                labels=as.character(labels),
                feature_importance=feature_score,
                nIter = i))
}



