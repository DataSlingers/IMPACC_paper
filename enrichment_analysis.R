library(limma)
library(ggplot2)
library(forcats)
library(org.Hs.eg.db)
library(clusterProfiler)
############## maker genes 
##########
### yan
##########https://www.nature.com/articles/nsmb.2660
yan=readRDS('data/yan.rds')
clus_yan=readRDS('results/res_yan.rds')
gene = rownames(yan$sc_cnt)
######## ENRICHMENT go analysis 
impo=clus_yan[[1]]$IMPACC$feature_importance
sorted_gene_IMPACC = tolower(gene[sort(impo, index.return=TRUE, decreasing=TRUE)$ix])
sorted_impo = sort(impo, decreasing=TRUE)
cutoff = mean(sorted_impo)+sd(sorted_impo)
gg_impacc = names(which(sorted_impo>cutoff))
write.csv(gg_impacc,'de_impacc.csv')

de_cs3 =  clus_yan[[1]]$SC3$info$sc3_7_de_padj
gg_sc3 = gene[which(de_cs3<0.05)]
gg_km = gene[which(clus_yan[[1]]$sparseKM$ws>0)]

ego <- enrichGO(gene = toupper(gg), OrgDb='org.Hs.eg.db', 
                keyType= 'SYMBOL',ont = "ALL", minGSSize = 1, pvalueCutoff = 0.05, qvalueCutoff = 1) 
ego_sc3 <- enrichGO(gene = toupper(gg_sc3), OrgDb='org.Hs.eg.db', 
                    keyType= 'SYMBOL',ont = "ALL", minGSSize = 1, pvalueCutoff = 0.05, qvalueCutoff = 1) 
ego_km <- enrichGO(gene = toupper(gg_km), OrgDb='org.Hs.eg.db', 
                   keyType= 'SYMBOL',ont = "ALL", minGSSize = 1, pvalueCutoff = 0.05, qvalueCutoff = 1) 

## fold enrichment is defined as the ratio of the frequency of input genes annotated in a term to the frequency of all genes annotated to that term, and it is easy to calculate by dividing geneRatio by BgRatio.
egoo <- mutate(ego, foldEnrich = as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("/\\d+", "", BgRatio)))

imp_go = cbind(ego$Description,round(ego$p.adjust,5),ego$GeneRatio,ego$BgRatio,ego$Count)
sc3_go = cbind(ego_sc3$Description,round(ego_sc3$p.adjust,5),ego_sc3$GeneRatio,ego_sc3$BgRatio,ego_sc3$Count)
km_go = cbind(ego_km$Description,round(ego_km$p.adjust,5),ego_km$GeneRatio,ego_km$BgRatio,ego_km$Count)
colnames(imp_go)=colnames(sc3_go)=colnames(km_go)=c('Pathway','p-value','Gene Ratio','Bg Ratio','Count')
kable(km_go, "latex",  booktabs = T) 
write.csv(sc3_go,'GO_SC3.csv')


go=ggplot(egoo, showCategory = 10,
          aes(foldEnrich, fct_reorder(Description, foldEnrich))) +
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=p.adjust, size = Count)) +
    scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                          trans = "log10",
                          guide=guide_colorbar(reverse=TRUE, order=1)) +
    scale_size_continuous(range=c(2, 10)) +
    xlab("Fold enrichment") +
    ylab('Pathway')+
    ggtitle("IMPACC GO Enrichment Analysis on Yan")
ggsave('yan_go.png',go,height = 6,width = 9)



####################################
## plot heatmap of marker genes
####################################
library(pheatmap)
library(cowplot)
sorted_gene_IMPACC = gene[sort(impo, index.return=TRUE, decreasing=TRUE)$ix]
sel  =yan$sc_cnt[sorted_gene_IMPACC[1:50],]

hc = hclust(as.dist(1-clus_yan$IMPACC$mat),method='ward.D')
anno = data.frame(yan$sc_label)
rownames(anno)=colnames(yan$sc_cnt)
anno_colors <- list(Cluster = c("darkorange2", "purple1", "steelblue1",
                                "seagreen4","yellow3",'peachpuff3','magenta2'))
names(anno_colors$Cluster) = unique(yan$sc_label)
names(anno)='Oracle Clusters'
 
## cluster the genes 
kmeans = kmeans(scalematrix(sel),length(unique(yan$sc_label)))$cluster
 

yan_heat = yan_heat[sort(kmeans,index.return=TRUE)$ix,]

p = pheatmap(yan_heat,
             cluster_cols = hc,
             cluster_rows = F,
             cutree_cols = 7,
             gaps_row = c(10, 20),
             annotation_col = anno,
             annotation_names_col = FALSE,
             annotation_colors = anno_colors,
             show_colnames = F,
             legend_labels = F,
             silent = T)

plot_grid(p$gtable, nrow = 1)
ggsave('yan_heatmap.png',plot_grid(p$gtable, nrow = 1),height = 9,width = 9)
