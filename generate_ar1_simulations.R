args = commandArgs(trailingOnly=TRUE)

source('simulation_function.R')
source('load_packages.R')
source('all_methods.R')
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
j = as.integer(args[1])
#snr = as.integer(args[2])

sim_sparse=list()


for(snr in c(1:8)){
    set.seed(snr*98+11*j)
    print(snr)
 #  if (length(sim_sparse[[as.character(SNR[snr])]])==0){
    #################################################
    ########### simulate sparse data set
    ##################################################
     sim_sparse[[as.character(snr)]] = simulate_sparse(n = c(20,80,120,280),
                                                        p=p,
                                                        signal =list(c(muu[snr],muu[snr]),
                                                                     c(muu[snr],-muu[snr]),
                                                                     c(-muu[snr],muu[snr]),
                                                                    c(-muu[snr],-muu[snr])),
                                                        num_cluster=num_clu,
                                                        nr_other_vars=M-p)
     
	saveRDS(sim_sparse,paste0('data/sim_block_',j,'.rds'))


    tic.clearlog()
    ###########################################
    ############## clustering results
    ############################################
   clu_fea_sel[[as.character(SNR[snr])]]=all_methods(sim_sparse[[as.character(SNR[snr])]],spec = F)
saveRDS(clu_fea_sel,paste0('results/res_sim_ar_',j,'.rds'))
  

}

## save to results folder
#saveRDS(clu_fea_sel,paste0('results/res_sim_ar_',j,'.rds'))


