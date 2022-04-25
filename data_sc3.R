######################
###### PANCAN tumor cell
########################
########read data and labels
dat = read.csv('data/TCGA-PANCAN-HiSeq/data.csv')
rownames(dat)=dat[,1]
dat = dat[,2:ncol(dat)]
l = read.csv('data/TCGA-PANCAN-HiSeq/labels.csv')
label = l[[2]]
hiseq=list()
hiseq$raw = t(dat)
hiseq$sc_cnt = normalized_dat=data.frame(log2(1+t(dat)))
hiseq$sc_label = label

hiseq=readRDS('data/hiseq.rds')
#############
##### bias https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57249
###############
library(data.table)
library(Matrix)
label = c(rep('Zygote',9),rep('Two-cell Embryo',20),
           rep('Four-cell Embryo',20)) ## from GEO website
dat = fread('data/GSE57249_fpkm.txt',header=T)
gene = dat$ID
dat = data.frame(dat)[,2:50]
rownames(dat)=gene
normalized_dat=data.frame(log2(1+dat))

biase = list()
biase$raw = dat
biase$sc_cnt = normalized_dat
biase$sc_label = label
saveRDS(biase,'biase.rds')

###########################
##### Goolam
#############################
metadata=fread('https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3321/E-MTAB-3321.sdrf.txt',header=T,fill=TRUE)
dat = data.frame(fread('data/Goolam_et_al_2015_count_table.tsv',header=T,fill=TRUE))
coln = colnames(dat)[1:124]
rownames(dat) = dat[,1]
dat = dat[,2:125]
colnames(dat)=coln
normalized_dat=data.frame(log2(dat+1))
##quality control: remove cells with low expression?
label = metadata$`Characteristics[developmental stage]`
label = label[match(sub("X", "",coln),metadata$`Source Name`)]

goolam = list()
goolam$raw = dat
goolam$sc_cnt = normalized_dat
goolam$sc_label = label
saveRDS(goolam,'goolam.rds')



############### yan https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36552
library(stringr)
dat = NULL
setwd("data/GSE36552")
files <- sort(list.files(pattern = "*.txt.gz$"))
for (f in files){
    print(f)
    # f = "GSM1278043_fibroblast_21_BxC_expression.txt.gz"
    # cell_name = regmatches(f, regexec("GSM1112490_\\s*(.*?)\\s*_expression.txt.gz",f))[[1]][2]
    x=read.table(gzfile(f),sep="\t", row.names=NULL,fill=T)
    if (is.null(dat)){
        dat = x[2:nrow(x), c(1,2)]
    }else{
        # x = data.frame(x[,c(1,4)])
        
        dat = merge(dat, x[2:nrow(x), c(1,2)], by = 'V1',all=TRUE)
        
    }
}
rownames(dat) =make.names(dat[,1], unique=TRUE) 
dat=dat[,2:ncol(dat)]
dim(dat)
cell_name = sapply(strsplit(files,'_'),"[[",1)
colnames(dat) = cell_name
dat[is.na(dat)]= 0
setwd("~/Documents/Documents - Alice/Rice/research/mini_clustering/PLOS review")

library(data.table)
metadata = fread('data/GSE36552_series_matrix.txt',header=T,fill=TRUE)
label = (sub("#.*", "", metadata[52,2:ncol(metadata)])) 

dat = dat[,label!='hESC passage']
label = label[label!='hESC passage']


yan = list()
yan$raw = dat
yan$sc_cnt = log2(1+dat)
yan$sc_label = label

saveRDS(yan,'yan.rds')

#######
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

################################
####### run clustering methods
###############################
source('all_methods.r')
res_pancan=all_methods(hiseq,nIter=500,sparseHC=FALSE,SC3=FALSE) ##sparseHC is too slow 
res_biase=all_methods(biase,nIter=500) ##sparseHC is too slow 
res_goolam=all_methods(goolam,nIter=500) ##sparseHC is too slow 
res_yan=all_methods(yan,nIter=500) ##sparseHC is too slow 
res_coil=all_methods(coil,nIter=500,sparseHC=FALSE,SC3=FALSE)

