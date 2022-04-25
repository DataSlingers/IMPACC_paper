###
source('load_packages.R')
source('result_functions.R')

snr = c(1,2,3,4,5,6,7,8)

files <- Sys.glob('results/res_sim_ar_*')
res_sim=vector(mode = "list", length =10)
for (snr in c(1:8)){
  for (j in c(1:10)){
    res_sim[[j]][[snr]] = readRDS(paste0('results/res_sim_ar_',snr,'_',j,'.rds'))
  }
}

# files <- Sys.glob('results/res_sim_block*')
# files
# 
# res_sim=list()
# for (j in c(1:length(files))){
#     res_sim[[j]] = readRDS(paste0('results/res_sim_block_',j,'.rds'))
# }



### Plot ARI, F1, computation time results 
res = vector(mode = "list", length = 4)
names(res)=c('F1_oracle','F1_dd','aris','time')

snr = c(1:8)

for (j in 1:length(res_sim)){
  mini=res_sim[[j]]
  cm=list()

  ## oracle , know number of featrues 
  F1_o = array(0,dim = c(8,3))
  F1_d = array(0,dim = c(8,3))
  time = array(0,dim = c(8,3))
  # true=c(rep(1,100),rep(0,9900))
  true=c(rep(1,25),rep(0,4975))
  
  for (v in 1:8){
    ### IMPACC
    ## find stop point 
    fea = mini[[v]]$IMPACC$feature_importance
    cut = mean(fea)+sd(fea) ##>1 all ok
    # cut = kde_based_pi_thr(fea)
    prediction_o = prediction_d = rep(0,length(true))
    prediction_o[sort(fea, index.return=TRUE, decreasing=TRUE)$ix[1:25]]=1
    prediction_d[fea>cut]=1
    sum(prediction_d)
    F1_o[v,1] = get_F1(prediction_o,true)
    F1_d[v,1] = get_F1(prediction_d,true)
    time[v,1] = mini[[v]]$IMPACC$time
    
    
    ## sparsekM
    ws = mini[[v]][['sparseKM']]$ws
    prediction_o = prediction_d = rep(0,length(true))
    if ((sum(ws)>0)<25){
      prediction_o[ws>0]=1
    }else{
      prediction_o[sort(ws, index.return=TRUE, decreasing=TRUE)$ix[1:25]]=1
    }
    prediction_d[ws>0]=1
    F1_o[v,2] = get_F1(prediction_o,true)
    F1_d[v,2] = get_F1(prediction_d,true)
    time[v,2] = mini[[v]]$sparseKM$time
    
    ## sparseHC
    ws =mini[[v]][['sparseHclust']]$ws 
    prediction_o = prediction_d = rep(0,length(true))
    if ((sum(ws)>0)<25){
      prediction_o[ws>0]=1
    }else{
      prediction_o[sort(ws, index.return=TRUE, decreasing=TRUE)$ix[1:25]]=1
    }
    prediction_d[ws>0]=1

    F1_o[v,3] = get_F1(prediction_o,true)
    F1_d[v,3] = get_F1(prediction_d,true)
    time[v,3] = mini[[v]]$sparseHclust$time
  }
  time = cbind(time,
               (unlist(sapply(mini, function(x) x$MPCC$time))),
               (unlist(sapply(mini, function(x) x$consensus$time))),
               t(unlist(sapply(mini, function(x) x$regular$time))))
  time=time[,c(1,4:5,2:3,6:8)]
  aris = array(0,dim=c(8,11))
  aris[,1] = (unlist(sapply(mini, function(x) x$IMPACC$ARI_hc)))
  aris[,2] = (unlist(sapply(mini, function(x) x$IMPACC$ARI_spec)))
  aris[,3] = (unlist(sapply(mini, function(x) x$MPCC$ARI_hc)))
  aris[,4] = (unlist(sapply(mini, function(x) x$MPCC$ARI_spec)))
  aris[,5] = (unlist(sapply(mini, function(x) x$consensus$ARI_hc)))
  aris[,6] = (unlist(sapply(mini, function(x) x$consensus$ARI_spec)))
  
  aris[,7] = unlist(sapply(mini, function(x) x$sparseKM$ARI)) 
  aris[,8] = unlist(sapply(mini, function(x) x$sparseHclust$ARI)) 
  
  aris[,9:11] = t(unlist(sapply(mini, function(x) x$regular$ARI)))
  
  colnames(time) = c('IMPACC','MPCC','CSS','sparseKM','sparseHC','KMeans','HClust','KMedoid')
  colnames(aris) =c('IMPACC(HC)','IMPACC(Spec)','MPCC(HC)','MPCC(Spec)','CSS(HC)','CSS(Spec)','sparseKM','sparseHC','KMeans','HClust','KMedoid') 
  # colnames(time) = colnames(aris)= c('sparseKM','sparseHC','IMPACC','MPCC','Consensus','KMeans','HClust','KMedoid')
  rownames(F1_o)=rownames(F1_d)= rownames(time) = rownames(aris)= seq(1:8)
  F1_o[is.na(F1_o)]=0
  F1_d[is.na(F1_d)]=0
  colnames(F1_d) =colnames(F1_o) = c('IMPACC','sparseKM','sparseHC' )

  res$F1_dd[[j]] = F1_d 
  res$F1_oracle[[j]] = F1_o
  res$aris[[j]] = aris
  res$time[[j]] = log(time)
}


##### summarize the results
summ = get_summary(res)
F11_oracle = data.frame(summ[[1]])
F11_dd = data.frame(summ[[2]])
arii = data.frame(summ[[3]])
timee = data.frame(summ[[4]])

# arii$data = c(rep('F1 Data Driven',nrow(summ[[2]])),rep('ARI',nrow(summ[[3]])),rep('Computational Time (log2)',nrow(summ[[4]])))
colnames(arii)[1:2]=colnames(F11_oracle)[1:2]=colnames(F11_dd)[1:2]=colnames(timee)[1:2]=c('SNR','Method')
arii$SNR = as.numeric(as.character(arii$SNR))
arii$value[arii$value<0]=0

#### make plots 
cols = c('red1',"red3",
         'purple1','magenta2',
         'steelblue1','steelblue3',
         'seagreen4','seagreen3',
         'yellow3','peachpuff3','sandybrown')

p_ari_all <- ggplot(arii, aes(x=SNR, y=value,group=Method,col = Method,fill = Method))+
  geom_point(alpha=0.2,size = 0.8)+
  geom_smooth(alpha = 0.2,se=F)+
  theme_bw()+
  scale_y_continuous(limit=c(0,1))+
  ggtitle("ARI")+
  xlab('SNR')+
  ylab('ARI')+
  theme(plot.title = element_text(size=20,hjust = 0.5),
        axis.text=element_text(size=10*2),legend.position = 'bottom',
        legend.box.margin=margin(0),
        legend.text=element_text(size=14),axis.title = element_text(size = 20),
        legend.title=element_blank(),
        strip.text = element_text(size = 18))+
  scale_color_manual(values =cols)+
  scale_fill_manual(values=cols)
p_ari_all 
ggsave('simu_all_auto.png',p_ari_all,height = 7,width = 7)

## include only (HC) results
del = c(which(arii$Method=='IMPACC(Spec)'),
        which(arii$Method=='MPCC(Spec)'),
        which(arii$Method=='CSS(Spec)'))

arii2 = arii[-del,]  
arii2$Method=timee$Method
cols_t = cols[c(1,3,5,7:11)]

p_ari <- ggplot(arii2, aes(x=SNR, y=value,group=Method,col = Method,fill = Method))+
  geom_point(alpha=0.2,size = 0.8)+
  geom_smooth(alpha = 0.2,se=F)+
  theme_bw()+
  scale_y_continuous(limit=c(0,1))+
  xlab('SNR')+
  ylab('ARI')+
  ggtitle("ARI")+
  theme(plot.title = element_text(size=20,hjust = 0.5),
        axis.text=element_text(size=10*2),legend.position = 'bottom',
        legend.box.margin=margin(0),
        legend.text=element_text(size=18),axis.title = element_text(size = 20),
        legend.title=element_blank(),
        strip.text = element_text(size = 18))+
  scale_color_manual(values =cols_t)+
  scale_fill_manual(values=cols_t)

p_ari
## time 
p_time <- ggplot(timee, aes(x=SNR, y=value,group=Method,col = Method,fill = Method))+
  geom_point(alpha=0.2,size = 0.8)+
  geom_smooth(alpha = 0.2,se=F)+
  theme_bw()+
  xlab('SNR')+
  ylab('Computational Time (log)')+
  ggtitle("Computational Time")+
  theme(plot.title = element_text(size=20,hjust = 0.5),
        axis.text=element_text(size=10*2),legend.position = 'bottom',
        legend.text=element_text(size=18),axis.title = element_text(size = 20),
        legend.box.margin=margin(15),
        legend.title=element_blank(),
        strip.text = element_text(size = 18))+
  scale_color_manual(values =cols_t)+
  scale_fill_manual(values=cols_t)+
  scale_y_continuous(breaks = log(c(1,5,30,60,120,300,600)),
                     labels = c('1 s', '5 s', '30 s', '1 min','2 min', '5 min','10 min'))
p_time
## F1 
dim(F11_oracle)
F11 = rbind(F11_oracle,F11_dd)
F11$ltype = c(rep('Oracle',240),rep('Data Driven',240))
F11$method = paste0(F11$Method,' (',F11$ltype,')')
cols_f1 = c("red1","red3","seagreen4","seagreen2","seagreen3","seagreen1")

p_F1 <- ggplot(F11, aes(x=SNR, y=value,group=method,col = method,fill = method,linetype = ltype))+
  geom_point(alpha=0.2,size = 0.8)+
  geom_smooth(alpha = 0.2,se=F)+
  theme_bw()+
 #scale_y_continuous(limit=c(0,1))+
  coord_cartesian(ylim=c(0, 1))+
  xlab('SNR')+
  ylab('F1')+
  ggtitle("F1 score")+
  theme(plot.title = element_text(size=20,hjust = 0.5),
        axis.text=element_text(size=10*2),legend.position = 'bottom',
        legend.text=element_text(size=18),axis.title = element_text(size = 20),
        legend.box.margin=margin(30),
        legend.title=element_blank(),
        strip.text = element_text(size = 18))+
  scale_color_manual(values =cols_f1)+
  scale_fill_manual(values=cols_f1)
fitt =  ggplot_build(p_F1)$data[[2]] 

fitt$y[fitt$y<0]=0
fitt$y[fitt$y>1]=1
fitt$method = fitt$colour
fitt$ltype=fitt$linetype
fitt[fitt=='red1']='IMPACC (Data Driven)'
fitt[fitt=='seagreen4']='sparseKM (Data Driven)'
fitt[fitt=='seagreen3']='sparseHC (Data Driven)'
fitt[fitt=='red3']='IMPACC (Oracle)'
fitt[fitt=='seagreen2']='sparseKM (Oracle)'
fitt[fitt=='seagreen1']='sparseHC (Oracle)'
fitt[fitt=='solid']='Data Driven'
fitt[fitt=='22']='Oracle'

p_F1 <- ggplot(F11, aes(x=SNR, y=value,group=method,col = method,fill = method,linetype = ltype))+
  geom_point(alpha=0.2,size = 0.8)+
  geom_line(aes(x,y), data=fitt,size=1)+
  theme_bw()+
  #scale_y_continuous(limit=c(0,1))+
  coord_cartesian(ylim=c(0, 1))+
  xlab('SNR')+
  ylab('F1')+
  ggtitle("F1 score")+
  theme(plot.title = element_text(size=20,hjust = 0.5),
        axis.text=element_text(size=10*2),legend.position = 'bottom',
        legend.text=element_text(size=13),axis.title = element_text(size = 20),
        legend.box.margin=margin(30),
        legend.title=element_blank(),
        strip.text = element_text(size = 18))+
  scale_color_manual(values =cols_f1)+
  scale_fill_manual(values=cols_f1)
p_F1 <- p_F1 + guides(colour = guide_legend(nrow=3,byrow=TRUE,override.aes = list(linetype = c('solid','22','solid','22','solid','22'))))
p_F1 <- p_F1 + scale_linetype(guide = FALSE)


# grid.arrange(p_ari, p_time, p_F1,ncol =3)
gt <- arrangeGrob(p_ari, p_time, p_F1,                          # box plot and scatter plot
                  ncol = 3)
p <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C"), size = 18,
                  x = c(0, 0.34, 0.67), 
                 # x = c(0, 0.25,0.5, 0.75), 
                  y = c(1, 1, 1)) # Add labels
p
ggsave('simu_auto.png',p,height = 7,width = 21)
