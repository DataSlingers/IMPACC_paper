args = commandArgs(trailingOnly=TRUE)

source('load_packages.R')
source('all_methods.R')
####################################################
## simulate sparse data sets with SNR = 1,..,8, 
# number of signal features =25  
########################################################

library(splatter)
library(scater)
goolam  = readRDS('data/goolam.rds')
j = as.integer(args[1])

set.seed(123*j+5)
real <- goolam$raw[, sample(1:ncol(goolam$raw), 100)]
real <- real[rowSums(goolam$raw) > 0, ]
params <- splatEstimate(as.matrix(real))


sim_real=list()


for(snr in c(1:8)){
	seed = snr*98+11*j
    set.seed(seed)
    print(snr)
    #################################################
    ########### simulate sparse data set
    ##################################################
    sim.groups = splatSimulateGroups(
			params=params,	
			group.prob = c(0.15, 0.25,0.6),
                            #method = "groups",
                            nGenes = 10000,
                            batchCells = 500,
                            de.prob         = 0.03, ##controls the probability that a gene will be selected to be differentially expressed
                            de.facLoc       = snr/10,
                            de.facScale     = 0.3,
                        #    dropout.present = FALSE,
                            verbose = FALSE,
			    seed = seed) 
	sim = list()
    sim$raw = data.frame(counts(sim.groups))
	sim$sc_cnt = normalizeCounts(sim$raw)
	sim$sc_label = sim.groups$Group
	sim$info = rowData(sim.groups)

    
    sim_real[[as.character(snr)]] = sim 
   
   
    saveRDS(sim_real,paste0('data/sim_goolam_',j,'.rds'))
    tic.clearlog()
    ###########################################
    ############## clustering results
    ############################################
    clu_fea_sel[[as.character(SNR[snr])]]=all_methods(sim_real[[as.character(SNR[snr])]],spec = F)
	saveRDS(clu_fea_sel,paste0('results/res_sim_splatter_',j,'.rds'))
}
