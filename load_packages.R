# Package names
packages <- c("ggplot2", 'data.table',
              "matrixStats", 
              "mclust", "tictoc", "mvtnorm", "reshape2", "dplyr", "sparcl",
              'cowplot','gghighlight','gridExtra','ggpubr',
              'ConsensusClusterPlus',
              'Rtsne',
              'SC3',
              'scater',
              'mclust',
              'Matrix',
              'Seurat')

########### bio packages

package_bio = c('ConsensusClusterPlus','SC3',
                'scater','Seurat')
installed_packages <- package_bio %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(package_bio[!installed_packages])
}

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

install.packages()
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


package_bio = c('ConsensusClusterPlus')
installed_packages <- package_bio %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(package_bio[!installed_packages])
}

BiocManager::install('ConsensusClusterPlus')