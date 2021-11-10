##########################################
## data proprocessing
####################################

######################
###### PANCAN tumor cell
########################
########read data and labels
dat = read.csv('data/TCGA-PANCAN-HiSeq/data.csv')
rownames(dat)=dat[,1]
dat = dat[,2:ncol(dat)]
l = read.csv('data/TCGA-PANCAN-HiSeq/labels.csv')
l = l[[2]]
x1 = t(dat)
# Filter out cells with library size greater than 95-th percentile.
x2 <- x1[, which(colSums(x1) <= quantile(colSums(x1),0.95))] 
label = l[which(colSums(x1) <= quantile(colSums(x1),0.95))]
# Remove genes with mean expression less than 25-th percentile.
x3 <- x2[rowMeans(x2) >= quantile(rowMeans(x2),0.25), ]
# Remove genes with less than 15-th percentile non-zero cells.
x4 <- x3[rowSums(x3 != 0) >= quantile(rowSums(x3 != 0),0.15), ]
hiseq=list()
hiseq$sc_cnt = x4
hiseq$sc_label = label

##################################
################brain cells 
#################################
load('data/darmanis.rdata') ## brain cells 
sc = darmanis$sc_cnt
## conduct log2 transform to expression
sc_scale= log2(sc+1)
sc_label = darmanis$sc_label
K = length(unique(sc_label))

## combine labels of replicates 
new_label=sc_label
new_label[c(which(sc_label=='fetal_replicating'),which(sc_label=='fetal_quiescent'))]='fetal'
new_sc=sc_scale


brain=list()
brain$sc_cnt=new_sc
brain$sc_label=new_label


#############################
### neoplastic scRNA-seq
#############################

library(data.table)
library(Matrix)
metadata = fread('data/E-GEOD-84465-normalised-files/ExpDesign-E-GEOD-84465.tsv',header=T)
labels = as.matrix(metadata[,c(1,16)])## 16 , 24
dat = readMM('data/E-GEOD-84465-normalised-files/E-GEOD-84465.aggregated_filtered_normalised_counts.mtx')
colnames(dat) = as.matrix(fread('data/E-GEOD-84465-normalised-files/E-GEOD-84465.aggregated_filtered_normalised_counts.mtx_cols',header = F))
rows= as.matrix(fread('data/E-GEOD-84465-normalised-files/E-GEOD-84465.aggregated_filtered_normalised_counts.mtx_rows',header = F))
rownames(dat)=rows[,1]
dat=log2(1+as.matrix(dat))
label= labels[match(colnames(dat),labels[,1]),2]

geod = list()
geod$sc_cnt = dat
geod$sc_label = label
#####################################
#### coil image data 
#######################################

library(R.matlab)
a = readMat("data/COIL20/COIL_1.mat")
b = readMat("data/COIL20/COIL_2.mat")

dat = cbind(a[[1]],a[[2]])
rownames(dat) = paste0('p',1:nrow(dat))
colnames(dat)=paste0('x',1:ncol(dat))
label = rbind(a[[3]],a[[4]])[,1]
coil=list()
coil$sc_cnt = dat
coil$sc_label = label


#####################
## model fitting, train on real data 
######################
## get ARI, computation time
source('all_methods.r')
res_pancan=all_methods(hiseq,nIter=200,sparseHC=F) ##sparseHC is too slow 
res_brain=all_methods(brain,nIter=200)

res_brain=all_methods(brain,nIter=200,CS3=T,seruat=T,DR_regular=T)

res_neo=all_methods(geod,nIter=200,sparseHC=F)
res_coil=all_methods(coil,nIter=200)

#####################################################
########## validation on real data, generate comparison
######################################################
####### clustering accuracy and running time
###############################################
source('result_functions.r')
get_ari_time(res_brain)
get_ari_time(res_pancan)
get_ari_time(res_neo)
get_ari_time(res_coil)

###################################################
######################## pathway analysis 
####################################################
get_top_kegg(res_brain$IMPACC$feature_importance,rownames(brain$sc_cnt))
get_top_kegg(res_brain$sparseKM$ws,rownames(brain$sc_cnt))
get_top_kegg(res_brain$sparseHclust$ws,rownames(brain$sc_cnt))


get_top_kegg(res_pancan$IMPACC$feature_importance,rownames(brain$sc_cnt))
get_top_kegg(res_pancan$sparseKM$ws,rownames(brain$sc_cnt))


get_top_kegg(res_neo$IMPACC$feature_importance,rownames(brain$sc_cnt))
get_top_kegg(res_neo$sparseKM$ws,rownames(brain$sc_cnt))




