library(mclust)
library(tictoc)
sampleBurnin <- function(d,num_partition,
                        p_sample,
                        p_sample_row){
  list_col = lapply(1:ceiling(num_partition/length(p_sample)), function(a) {
    g = sample(cut(seq(ncol(d)),ncol(d)*cumsum(c(0,p_sample))))
    partCol = split(1:ncol(d), g)
    names(partCol)=NULL
    partCol
  })
  list_col = unlist(list_col, recursive = FALSE)
  
  list_row = lapply(1:ceiling(num_partition/length(p_sample_row)), function(a) {
    g = sample(cut(seq(nrow(d)),nrow(d)*cumsum(c(0,p_sample_row))))
    partRow = split(1:nrow(d), g)
    names(partRow)=NULL
    partRow
  })
  list_row = unlist(list_row, recursive = FALSE)
  res = lapply(1:num_partition, function(i){
    list(submat=d[list_row[[i]],list_col[[i]]],
         subrows=list_row[[i]],
         subcols=list_col[[i]])
  })
  return(res)
}

########### update matrix that record joint clustering info (M^{h}),(N × N) connectivity matrix corresponding to dataset D(h) 
connectivityMatrix <- function(clusterAssignments, m, sampleKey){
  ##input: named vector of cluster assignments, matrix to add connectivities
  ##output: connectivity matrix

  names(clusterAssignments) <- sampleKey 
  #list samples by clusterId
  cls <- lapply(unique(clusterAssignments), function(i) as.numeric(names(clusterAssignments[clusterAssignments%in%i])))  
  for ( i in 1:length(cls)){
    nelts <- 1:ncol(m)
    cl <- as.numeric( nelts %in% cls[[i]] ) ## produces a binary vector
    ## cl = 1 if this obs is in cluster i 
    updt <- outer( cl, cl ) 
    #product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster
    ## (i,j) = 1 if (i,j) are in same cluster
    m <- m + updt
  }
  return(m)
}

cluster_algo = function(submat,h,distance,innerLinkage){
  dis = dist(t(submat),method = distance)
  this_cluster = hclust(dis, method=innerLinkage)
  this_cluster$height <- round(this_cluster$height, 6) 
  this_assignment = cutree(this_cluster,h=quantile(this_cluster$height,h))
  return(this_assignment=this_assignment)
}


sampleCols_RD <- function(d,
                          pSamp=NULL,
                          pRow=NULL,
                          weightsItem=NULL, ## vector 
                          weightsFeature=NULL){
  ## returns a list with the sample columns, as well as the sub-matrix & sample features (if necessary)
  space <-  ncol(d)
  sampleN <- floor(space*pSamp)
  sampCols <- sort(sample(space, sampleN, replace = FALSE, prob = weightsItem))
  ## sample rows 
  space = nrow(d)
  sampleN = floor(space*pRow)
  sampRows = sort( sample(space, sampleN, replace = FALSE, prob = weightsFeature) )
  this_sample <- d[sampRows,sampCols]
  # dimnames(this_sample) <- NULL
  return( list(submat=this_sample,
               subrows=sampRows,
               subcols=sampCols ))
}


rowsums = function(xx){
  if (is.null(dim(xx))){
    xx
  }else{
    rowSums(xx)
  }
}
rowmeans = function(xx){
  if (is.null(dim(xx))){
    xx
  }else{
    rowMeans(xx)
  }
}
pv_anova=function(X,y){
  # 
  # X = as.matrix(sample_x$submat)
  # y=as.factor(this_assignment)
  # 
  n_y=levels(y)
  ss_resi = rowsums(sapply(n_y, function(w) rowsums((X[,y==w]-rowmeans(X[,y==w]))^2)))
  ss_explained =  rowSums((X-rowMeans(X))^2)-ss_resi
  df1=length(n_y)-1
  df2 = ncol(X)-length(n_y)
  FF = (ss_explained/df1)/(ss_resi/df2 )
  return(pf(FF, df1, df2, lower.tail = FALSE))
}
