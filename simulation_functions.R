generateGaussianData <- function(n, p,center, sigma){
  data = rmvnorm(n, mean = center,sigma = sigma)
  data = data.frame(data)
  names(data) = paste0('p',seq(1:p))
  return(data)
}

###############################
### autoregressive simulation
###############################
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

simulate_ar1 = function(n,p,num_cluster,signal,nr_other_vars){
  dataset=NULL
  for (k in (1:num_cluster)){  
    mu = rep(signal[[k]],p/length(signal[[k]]))
    if (length(mu)<p){
      mu=c(mu,signal[[k]][1])
    }
    mu=c(mu,rep(0,nr_other_vars))
    
    sigma = ar1_cor(nr_other_vars+p, 0.5)
    
    
    mat = (generateGaussianData(n[k],p,mu, sigma))
    colnames(mat) = c(paste0('p',seq(1:p)),paste0('X',seq(1:nr_other_vars)))
    dataset[[k]] = mat
  }
  
  df <- (t(do.call("rbind", dataset)))
  label = rep(factor(seq(1,num_cluster)),times = n)
  ########## add noise 
  l = list()
  l$sc_cnt =data.frame(df)
  l$sc_label = label
  return(l)
}


###############################
### sparse block-diagonal simulation
###############################
simulate_sparse = function(n,p,num_cluster,signal,nr_other_vars){
  dataset=NULL
  for (k in (1:num_cluster)){  
    mu = rep(signal[[k]],p/length(signal[[k]]))
    if (length(mu)<p){
      mu=c(mu,signal[[k]][1])
    }
    
    sigma=diag(1,p)
    n_block = p/5
    for (a in c(1:(n_block-1))){
      l=a*5-4
      r=a*5
      sigma[l:r,l:r] = 0.5
    }
    
    diag(sigma)=1
    
    d1 = (generateGaussianData(n[k],p,mu, sigma))
    
    
    sigma_noise=diag(1,nr_other_vars)
    n_block = nr_other_vars/5
    for (a in c(1:(n_block-1))){
      l=a*5-4
      r=a*5
      sigma_noise[l:r,l:r] = 0.5
    }
    
    diag(sigma_noise)=1
    
    mat = generateGaussianData(n[k],nr_other_vars,rep(0,nr_other_vars), sigma_noise)
    # mat <- matrix(stats::rnorm(nr_other_vars*n[k]),n[k],nr_other_vars)
    colnames(mat) = paste0('X',seq(1:ncol(mat)))
    dataset[[k]] = cbind(d1,mat)
  }
  
  df <- (t(do.call("rbind", dataset)))
  label = rep(factor(seq(1,num_cluster)),times = n)
  ########## add noise 
  l = list()
  l$sc_cnt =data.frame(df)
  l$sc_label = label
  return(l)
}

###############################
### weak sparse block-diagonal simulation
###############################
simulate_weaksparse = function(n,p,num_cluster,signal,nr_other_vars){
  dataset=NULL
  mu_noise = rnorm(nr_other_vars,0,1)
  for (k in (1:num_cluster)){  
    mu = rep(signal[[k]],p/length(signal[[k]]))
    if (length(mu)<p){
      mu=c(mu,signal[[k]][1])
    }
    
    sigma=diag(1,p)
    n_weak = p/5
    for (a in c(1:(n_weak-1))){
      l=a*5-4
      r=a*5
      sigma[l:r,l:r] = 0.5
    }
    
    diag(sigma)=1
    
    d1 = (generateGaussianData(n[k],p,mu, sigma))
    
    sigma_noise=diag(1,nr_other_vars)
    n_weak = nr_other_vars/5
    for (a in c(1:(n_weak-1))){
      l=a*5-4
      r=a*5
      sigma_noise[l:r,l:r] = 0.5
    }
    
    diag(sigma_noise)=1
    
    
    ### weak sparse , noise with mu = 0.1 
    mat = generateGaussianData(n[k],nr_other_vars,mu_noise, sigma_noise)
    # mat <- matrix(stats::rnorm(nr_other_vars*n[k]),n[k],nr_other_vars)
    colnames(mat) = paste0('X',seq(1:ncol(mat)))
    dataset[[k]] = cbind(d1,mat)
  }
  
  df <- (t(do.call("rbind", dataset)))
  label = rep(factor(seq(1,num_cluster)),times = n)
  ########## add noise 
  l = list()
  l$sc_cnt =data.frame(df)
  l$sc_label = label
  return(l)
}

###############################
### no sparse block-diagonal simulation
###############################
simulate_nosparse = function(n,p,num_cluster,signal,nr_other_vars){
  dataset=NULL
  for (k in (1:num_cluster)){  
    mu = rep(signal[[k]],p/length(signal[[k]]))
    if (length(mu)<p){
      mu=c(mu,signal[[k]][1])
    }
    mu = mu+rnorm(p,0,0.1)
    sigma=diag(1,p)
    n_block = p/5
    for (a in c(1:(n_block-1))){
      l=a*5-4
      r=a*5
      sigma[l:r,l:r] = 0.5
    }
    
    diag(sigma)=1
    
    mat = (generateGaussianData(n[k],p,mu, sigma))
    # mat <- matrix(stats::rnorm(nr_other_vars*n[k]),n[k],nr_other_vars)
    colnames(mat) = paste0('p',seq(1:ncol(mat)))
    dataset[[k]] = mat
  }
  
  df <- (t(do.call("rbind", dataset)))
  label = rep(factor(seq(1,num_cluster)),times = n)
  ########## add noise 
  l = list()
  l$sc_cnt =data.frame(df)
  l$sc_label = label
  return(l)
}