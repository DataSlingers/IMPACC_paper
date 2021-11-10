setwd("~/Documents/Rice/research/mini_clustering/IMPACC/code")
source('simulation_function.r')
source('all_methods.r')
####################################################
## simulate sparse data sets with SNR = 1,..,8, 
# number of signal features =25  
########################################################
M = 5000
N = 500
p=25
num_clu=4
SNR = c(1,2,3,4,5,6,7,8)
muu  = SNR/sqrt(p)
clu_fea_sel=list()

for (j in 1:5){ 
  print(j)
  sim_sparse=list()
  clu_fea_sel[[j]]=list()
  # sim_sparse=readRDS(paste0('data/sim_block_',j,'.rds'))
  
  for(snr in c(1:8)){
    set.seed(snr*98+11*j)
    print(snr)
    #################################################
    ########### simulate sparse data set
    ##################################################
    sim_sparse[[as.character(SNR[snr])]] = simulate_sparse(n = c(20,80,120,280),
                                                       p=p,
                                                       signal =list(c(muu[snr],muu[snr]),
                                                                    c(muu[snr],-muu[snr]),
                                                                    c(-muu[snr],muu[snr]),
                                                                    c(-muu[snr],-muu[snr])),
                                                       num_cluster=num_clu,
                                                       nr_other_vars=M-p)


    tic.clearlog()
    ###########################################
    ############## clustering results
    ############################################
    clu_fea_sel[[j]][[as.character(SNR[snr])]]=all_methods(sim_sparse[[as.character(SNR[snr])]],spec = F,sparseKM = F,
                                                           sparseHC = F,regular =F)
  }
}

## save to results folder
saveRDS(clu_fea_sel,'results/res_sim_sparse.rds')
