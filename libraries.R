#if (!require("BiocManager", quietly = TRUE))
 # install.packages("Homo.sapiens")
#BiocManager::install("annotables");
#source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R");
#installAnRichment();
#BiocManager::install("Bioconductor/BiocFileCache")


library(tidyverse)
library(dplyr)
library(DESeq2)
library(edgeR)
library(GenomicFeatures)
library(biomaRt)
library(data.table)
library(plotly)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(cluster)
#library(genefilter)
library(umap)
#library(impute)
library(ggforce)
library(limma)
library(org.Hs.eg.db)
library(knitr)
library(reshape2)
library(gplots)
library(plyr)
library(statmod)
#library(MOFA2)
library(factoextra)
library(pheatmap)
#library(apeglm)
#library(EnhancedVolcano)
library(umap)
library(vcfR)
library(adegenet)
library(glmmTMB)
library(trainomeHelper)
library(trainomeMetaData)
library(DHARMa)
library(marginaleffects)
library(clusterProfiler)
library(GOfuncR)
library(AnnotationHub)
library(ensembldb)
library(annotables)

library(Homo.sapiens)
library(topGO)
library(KEGGREST)
library(DGCA)
library(Hmisc)
library(igraph)


