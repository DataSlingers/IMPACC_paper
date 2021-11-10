library(tictoc)
source('funmini.r')
######### random minipatch 
MPCC = function(d=NULL,
                h=0.95, ## cut internal tree at h quantile 
                reps=100, ## number of iterations 
                pItem=0.25, ## minipatch size for observations 
                pFeature=0.1, ## minipatch size for features 
                innerLinkage="ward.D", ## internal HC linkage 
                distance="manhattan",## internal HC distance
                verbose=T,
                early_stop=T, ## whether perform early stop criteria 
                corUse="everything"){
  ###########
  ## The MPCC function perform MPCC (MiniPatch Consensus Clustering) with uniform sampling scheme on both observations and features. 
  ## It returns final N by N consensus matrix 
  ############
  #######################################
  #######function to update consensus matrix 
  #######################################
  
  update_MPCC = function(){
    if(verbose){
      message(paste("  i =",i))
    }
    #################
    ## create MP
    ############# 'sampleCols_RD' is in funmini.R
    sample_x = sampleCols_RD(d, pItem, pFeature) 
    this_assignment = cluster_algo(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage)
    ##mCount stores number of times a sample pair was sampled together.
    mCount <<- connectivityMatrix(rep(1,length(sample_x$subcols)),
                                  mCount,
                                  sample_x$subcols) 
    
    ##ML stores number of times a sample pair was CLUSTERED together.
    ml <<- connectivityMatrix(this_assignment,
                              ml,
                              sample_x$subcols)
    this_Co = ml / mCount
    this_Co[mCount==0]=0
    CoAsso<<-this_Co
  }
  
  ##############################
  ############ initialize
  n=ncol(d)
  CoAsso=mCount=ml=matrix(c(0),ncol=n,nrow=n)
  
  if (early_stop == F){
    ###############################################################
    ### if not perform early stopping, run given number of iterations 
    ################################################################
    for (i in 1:(reps)){
      ## update consensus matrix 
      update_MPCC()
    }
  }else{
    conf_q=NULL ## record quantile of confusions  
    continue=T
    i=1
    while (continue==T & i<reps){
      update_MPCC() ## update css matrix
      conf_q=c(conf_q,quantile(rowMeans(CoAsso*(1-CoAsso)),0.9)) ## record 90% quantile of the confusion
      ########################
      ### stop the iteration if change of 90% quantile is less than 0.00001 for past 5 iterations 
      #######################
      if (i>5){ 
        ## find changes of quantile 
        cm = sapply(c(i:(i-4)), function(x) abs(conf_q[x]-conf_q[x-1]))
        if (max(cm)<0.00001){
          continue=F
          message(paste0('Stop at iteration ',i))
        }else{
          continue=T
        }
      }
      i=i+1
    }
  }
  message("end fraction")
  return(list(consensus = CoAsso,stop_point=i))
}





IMPACC = function(d=NULL,
                  h=0.95, ## cut internal tree at h quantile 
                  reps=100, ## number of iterations 
                  pItem=0.25, ## minipatch size for observations 
                  pFeature=0.1, ## minipatch size for features 
                  innerLinkage="ward.D", ## internal HC linkage 
                  distance="manhattan",## internal HC distance
                  verbose=T,
                  early_stop=T, ## whether perform early stop criteria 
                  E= 3, ## number of epochs in burn-in stage
                  qI=0.95, ## observations with greater than (qI) quantile of weights are classified to the set of high uncertain observations in adaptive sampling
                  qF=1, ## features with weights greater than (mean+qF*std) are classified to the set of high importance genes in adaptive sampling
                  alpha_I=0.5, ## learning rate for updating observation weights  
                  alpha_F=0.5, ## learning rate for updating feature weights  
                  pp=0.05, ## feature support threshold, a feature will be add to feature support if its p-value in ANOVA is smaller than pp
                  corUse="everything"){
  ###########
  ## The IMPACC function perform IMPACC (Interpretable MiniPatch Adaptive Consensus Clustering) with adaptive sampling scheme on both observations and features. 
  ## It returns (1)consensus: final N by N consensus matrix; (2):feature_importance: feature importance scores
  ############
  
  update_IMPACC = function(){  
    if(verbose){
      message(paste("  i =",i))
    }
    ### update observation weights
    confusion = rowMeans(CoAsso*(1-CoAsso))
    ww = (i+length(sample_Burn))/subsample_o*confusion
    wi <<- alpha_I*wi+(1-alpha_I)*(ww/sum(ww))
    wi_p<<- alpha_F*wi_p+(1-alpha_F)*feature_score
    
    ## subsample 
    sample_x = sample_EEprob_EEprob(d, pItem, pFeature,pi_item=pi_item[i],pi_feature=pi_feature[i],
                                    weightsItem = wi, weightsFeature = wi_p,qI=qI,qF=qF)
    subsample_o<<- subsample_o+colnames(d)%in%colnames(d)[sample_x$subcols]
    subsample_f<<- subsample_f+rownames(d)%in%rownames(d)[sample_x$subrows]
    ############################
    ## clustering 
    ######################## 'cluster_algo' is in funmini.R
    this_assignment = cluster_algo(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage)
    
    #######################################
    ######### for feature importance
    #######################################
    ### 'pv_anova' is in funmini.R
    
    pvalu = pv_anova(as.matrix(sample_x$submat),as.factor(this_assignment))
    pv = quantile(na.omit(pvalu),pp)
    pvalue = pvalu<=pv
    feature_support[rownames(sample_x$submat)[which(pvalue)]]<<- 1+ feature_support[rownames(sample_x$submat)[which(pvalue)]]
    feature_score<<- feature_support/subsample_f
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
  pi_item= pi_feature = seq(0.5,1,by=0.5/reps)
  
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
  sample_Burn = sampleBurnin(d,num_partition,p_sample,p_sample_row) ## 'sampleBurnin' is in funmini.R
  subsample_o=rep(0,ncol(d))
  subsample_f=rep(0,nrow(d))
  feature_support = rep(0,nrow(d)) ## feature support 
  names(feature_support)= rownames(d)
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
    ##########################################
    ######## for feature importance
    ###########################################
    pvalu = pv_anova(as.matrix(sample_x$submat),as.factor(this_assignment))
    pv = quantile(na.omit(pvalu),pp)
    pvalue = pvalu<=pv
    feature_support[rownames(sample_x$submat)[which(pvalue)]]=1+ feature_support[rownames(sample_x$submat)[which(pvalue)]]
    feature_score=feature_support/(subsample_f)
    feature_score[subsample_f==0]=0
    ###################################
    ## 'connectivityMatrix' is in funmini.R
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

  ############################################
  ## adaptive stage
  #####################################################
  ##message("adaptive stage")
  
  if (early_stop == F){
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
      if (i>5){
        ## find change of quantile 
        cm = sapply(c(i:(i-4)), function(x) abs(conf_q[x]-conf_q[x-1]))
        if (max((cm))<0.00001){
          continue=F
          message(paste0('Stop at iteration ',i+num_partition))
        }else{
          continue=T
        }
      }
      i=i+1
    }
  }
  
  return(list(consensus = CoAsso,
              feature_importance=feature_score))
}



update_IMPACC2 = function(i, this_coA,verbose=T){  
  if(verbose){
    message(paste("  i =",i))
  }
  ### update observation weights
  confusion = rowMeans(this_coA*(1-this_coA))
  ww = (i+length(sample_Burn))/subsample_o*confusion
  wi<<-alpha_I*wi+(1-alpha_I)*(ww/sum(ww))
  wi_p<<-alpha_F*wi_p+(1-alpha_F)*feature_score
  
  ## subsample 
  sample_x = sample_EEprob_EEprob(d, pItem, pFeature,pi_item=pi_item[i],pi_feature=pi_feature[i],
                                  weightsItem = wi, weightsFeature = wi_p,qI=qI,qF=qF)
  subsample_o<<-subsample_o+colnames(d)%in%colnames(d)[sample_x$subcols]
  subsample_f<<-subsample_f+rownames(d)%in%rownames(d)[sample_x$subrows]
  ############################
  ## clustering 
  ######################## 'cluster_algo' is in funmini.R
  this_assignment = cluster_algo(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage)
  
  #######################################
  ######### for feature importance
  #######################################
  ### 'pv_anova' is in funmini.R
  
  pvalu = pv_anova(as.matrix(sample_x$submat),as.factor(this_assignment))
  pv = quantile(na.omit(pvalu),pp)
  pvalue = pvalu<=pv
  feature_support[rownames(sample_x$submat)[which(pvalue)]]<<-1+ feature_support[rownames(sample_x$submat)[which(pvalue)]]
  feature_score <<-feature_support/subsample_f
  #######################################
  
  ##########################################
  ##mCount stores number of times a sample pair was sampled together.
  mCount<<-connectivityMatrix(rep(1,length(sample_x$subcols)),
                              mCount,
                              sample_x$subcols) 
  
  ##ml stores number of times a sample pair was sampled together.
  ml<<-connectivityMatrix(this_assignment,
                          ml,
                          sample_x$subcols)
  
  CoAsso<<-ml / mCount
  CoAsso[mCount==0]=0
}

#############################################################
############# adaptive OBSERVATION sampling scheme for IMAPCC
##############################################################
sample_EEprob_EEprob <- function(d,
                                 pSamp=NULL,
                                 pRow=NULL,
                                 weightsItem=NULL, ## vector 
                                 weightsFeature=NULL,
                                 pi_item=0.5,
                                 pi_feature=0.5,
                                 qI=0.95,
                                 qF=0.95 ){
  
  
  
  ## returns a list with the sample columns, as well as the sub-matrix & sample features (if necessary)
  space <-  ncol(d)
  sampleN <- floor(space*pSamp)
  upper = which(weightsItem>=quantile(weightsItem,qI))
  sampleN1 = ceiling(min(sampleN,pi_item*length(upper)))
  sampCols1 <- sort(sample(upper, sampleN1, replace = FALSE,prob = weightsItem[upper]))
  sampCols2 <- sort(sample((1:space)[-upper], sampleN-sampleN1, replace = FALSE))
  sampCols = sort(c(sampCols1,sampCols2))
  
  
  this_sample <- sampRows <- NA
  
  ## sample rows 
  space = nrow(d)
  sampleN = floor(space*pRow)
  if (qF<1&qF>0){
    upper = which(weightsFeature>=quantile(weightsFeature,qF))
  }else{
    upper = which(weightsFeature>=mean(weightsFeature)+qF*sd(weightsFeature))
  }
  sampleN1 = ceiling(min(sampleN,pi_feature*length(upper)))
  # if (sampleN1>1){
  sampRows1 <- sort(sample(upper, sampleN1, replace = FALSE,prob = weightsFeature[upper]))
  sampRows2 <- sort(sample((1:space)[-upper], sampleN-sampleN1, replace = FALSE))
  sampRows = sort(c(sampRows1,sampRows2))
  # }else{
  #   sampRows <- sort(sample((1:space), sampleN, replace = FALSE))
  # }
  
  
  this_sample <- d[sampRows,sampCols]
  # dimnames(this_sample) <- NULL
  return( list(submat=this_sample,
               subrows=sampRows,
               subcols=sampCols ))
}


#########################################
####### obtain final clustering results 

get_cluster = function(ml,K,finalalgorithm,linkage){
  return(tryCatch({
    if(finalalgorithm=='hclust'){
      hc=hclust(as.dist(1-ml),method=linkage)
      ct = cutree(hc,K)
    }else if (finalalgorithm=='kmedoid'){
      ct = pam(as.dist(1-ml), K)$clustering
    }else if (finalalgorithm=='spectral'){
      ct <- SNFtool::spectralClustering(ml, K)
      
      # ct <- tryCatch(SNFtool::spectralClustering(ml, K), error=function(err) 0)
    }else if (finalalgorithm == 'convex'){
      cap=NULL
      mds=cmdscale(1-ml)
      cap = CARP(mds)
      ct= cutree(cap$dendrogram,K)
      
    }else if (finalalgorithm == 'AP'){
      apres <- apcluster(ml)
      aggres <- aggExCluster(ml, apres)
      res = cutree(aggres, K)
      ct = rep(0,ncol(ml))
      for (kk in 1:K){
        ct[res[[kk]]]=kk
      }
    }
    # names(ct) = colnames(d)
    return(ct)
  }, error=function(e) NA))
}

get_ari = function(CoAsso,finalalgorithm,sc_label,KK,linkage){
  cluster_result=get_cluster(CoAsso,KK,finalalgorithm,linkage)
  return(tryCatch(adjustedRandIndex(cluster_result,sc_label), error=function(err) 0))
}

#################################