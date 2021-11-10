###
source('load_packages.r')
source('result_functions.r')

snr = c(1,2,3,4,5,6,7,8)
res_sim = readRDS('results/res_sim_sparse.rds')

### Plot ARI, F1, computation time results 
res = vector(mode = "list", length = 3)
names(res)=c('F1','aris','time')


## oracle, know number of featrues 
F1 = array(0,dim = c(length(snr),3))
true=c(rep(1,25),rep(0,4975))
for (i in c(1:length(res_sim))){
  for (v in 1:length(snr)){
    ## sparsekM
    bcss = res_sim[[i]][[v]][['sparseKM']]$ws
    cut = mean(bcss)+5*sd(bcss)
    prediction = rep(0,length(true))
    prediction[bcss>cut]=1
    F1[v,2] = get_F1(prediction,true)
    ## sparseHC
    ws =res_sim[[i]][[v]][['sparseHclust']]$ws 
    cut = mean(ws)+5*sd(ws)
    prediction = rep(0,length(true))
    prediction[ws>cut]=1
    F1[v,3] = get_F1(prediction,true)
    ### IMPACC
    ## find stop point 
    fi =res_sim[[i]][[v]]$IMPACC$feature_importance 
    cut = mean(fi)+5*sd(fi)
    prediction = rep(0,length(true))
    prediction[which(fi>cut)]=1
    F1[v,1] = get_F1(prediction,true)
  }
  
  time = cbind(
    (unlist(sapply(res_sim[[i]], function(x) x$IMPACC$time))),
    (unlist(sapply(res_sim[[i]], function(x) x$MPCC$time))),
    (unlist(sapply(res_sim[[i]], function(x) x$sparseKM$time))),
    (unlist(sapply(res_sim[[i]], function(x) x$sparseHclust$time))),
               (unlist(sapply(res_sim[[i]], function(x) x$consensus$time))),
               t(unlist(sapply(res_sim[[i]][1:8], function(x) x$regular$time))))

  aris = cbind(
    (unlist(sapply(res_sim[[i]], function(x) x$IMPACC$ARI_hc))),
    (unlist(sapply(res_sim[[i]], function(x) x$IMPACC$ARI_spec))),
    
    (unlist(sapply(res_sim[[i]], function(x) x$MPCC$ARI_hc))),
    (unlist(sapply(res_sim[[i]], function(x) x$MPCC$ARI_spec))),
    (unlist(sapply(res_sim[[i]], function(x) x$consensus$ARI_hc))),
    (unlist(sapply(res_sim[[i]], function(x) x$consensus$ARI_spec))),
    unlist(sapply(res_sim[[i]], function(x) x$sparseKM$ARI)),
    unlist(sapply(res_sim[[i]], function(x) x$sparseHclust$ARI)),
    t(unlist(sapply(res_sim[[i]], function(x) x$regular$ARI))))  

  colnames(time) = c('IMPACC','MPCC','CSS','sparseKM','sparseHC','KMeans','HClust','KMedoid')
  colnames(aris) =c('IMPACC(HC)','IMPACC(Spec)','MPCC(HC)','MPCC(Spec)','CSS(HC)','CSS(Spec)','sparseKM','sparseHC','KMeans','HClust','KMedoid') 
  rownames(time) = rownames(aris)= seq(1:8)
  F1[is.na(F1)]=0
  colnames(F1) = c('IMPACC' ,'sparseKM','sparseHC')
  rownames(F1)=rownames(aris)=rownames(time) =snr
  
  res$F1[[i]] = F1
  res$aris[[i]]= aris
  res$time[[i]] = log(time)
}


##### summarize the results
summ = get_summary(res)

F11 = data.frame(summ[[1]])
arii = data.frame(summ[[2]])
timee = data.frame(summ[[3]])

# arii$data = c(rep('F1 Data Driven',nrow(summ[[2]])),rep('ARI',nrow(summ[[3]])),rep('Computational Time (log2)',nrow(summ[[4]])))
colnames(arii)[1:2]=colnames(F11)[1:2]=colnames(timee)[1:2]=c('SNR','Method')
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
# ggsave('simu_all.png',p_ari_all,height = 7,width = 7)

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
  theme(plot.title = element_text(size=20,hjust = 0.5),
        axis.text=element_text(size=10*2),legend.position = 'bottom',
        legend.box.margin=margin(0),
        legend.text=element_text(size=18),axis.title = element_text(size = 20),
        legend.title=element_blank(),
        strip.text = element_text(size = 18))+
  scale_color_manual(values =cols_t)+
  scale_fill_manual(values=cols_t)


## time 
p_time <- ggplot(timee, aes(x=SNR, y=value,group=Method,col = Method,fill = Method))+
  geom_point(alpha=0.2,size = 0.8)+
  geom_smooth(alpha = 0.2,se=F)+
  theme_bw()+
  xlab('SNR')+
  ylab('Computational Time (log)')+
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

## F1 
p_F1 <- ggplot(F11, aes(x=SNR, y=value,group=Method,col = Method,fill = Method))+
  geom_point(alpha=0.2,size = 0.8)+
  geom_smooth(alpha = 0.2,se=F)+
  theme_bw()+
  scale_y_continuous(limit=c(0,1))+
  xlab('SNR')+
  ylab('F1')+
  theme(plot.title = element_text(size=20,hjust = 0.5),
        axis.text=element_text(size=10*2),legend.position = 'bottom',
        legend.text=element_text(size=18),axis.title = element_text(size = 20),
        legend.box.margin=margin(30),
        legend.title=element_blank(),
        strip.text = element_text(size = 18))+
  scale_color_manual(values =cols_t[c(1,4:5)])+
  scale_fill_manual(values=cols_t[c(1,4:5)])


gt <- arrangeGrob(p_ari, p_time, p_F1,                             # box plot and scatter plot
                  ncol = 3)
p <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C"), size = 18,
                  x = c(0, 0.34, 0.67), y = c(1, 1, 1)) # Add labels
p

# ggsave('simu.png',p,height = 7,width = 21)
