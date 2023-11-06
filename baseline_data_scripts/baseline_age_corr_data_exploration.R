# Exploration fo the differentially expressed lncRNA by age and correlated 
# protein coding data of from baseline data

source("./libraries.R")
source("./R/Trainome_functions.R")



#Load the correlation data
#contains the diffe expressed lncRNAs by age, and the protein-coding genes they corelate with


age_base_df <- readRDS("./data/baseline_data/Correlation_results_of_significance_for_age.RDS") 
  

age_abs_df <- age_base_df%>%
  dplyr::filter(correlation_coefficient >= 0.5 | correlation_coefficient < -0.5)#%>%
  ggplot(aes(y = p_value, x = correlation_coefficient)) + geom_point()


#store the coexpressed pairs in a g object

#g <- graph.data.frame(age_abs_df[, 1:2], directed = F)

#calculate the degree and betweeness objects respectively


# degree <- igraph::degree(g)
# 
# betweenness <- igraph::betweenness(g)
# 
# Node_nw_st <- data.frame(degree, betweenness)
# 
# 
# Rank_stat <- rowMeans(cbind(rank(Node_nw_st[,1]), rank(Node_nw_st[,2])))
# Node_nw_st <- cbind(Node_nw_st, Rank_stat)
# 
# plot(g)
# 
# #Which transcript has the highest degree
# which.max(degree)

#Save as txt file for cytoscape
# write.table(Node_nw_st, file = "./data/baseline_data/nodes_baseline_age.txt", sep = "\t",
#             col.names = NA, quote = F)

#save the filtered correlation table as well
# write.table(age_abs_df, file = "./data/baseline_data/filtered_correlation_age_absolutes.txt",
#             sep= "\t", row.names = F, quote = F)

#get the ensemble gene IDs as gsea doesnt seem to work with transcript names
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)


#extract transcript biotypes and transcript names
results_end_1 <- getBM(attributes = c("ensembl_gene_id",
                                      "external_gene_name", 
                                      "ensembl_gene_id_version",
                                      "external_transcript_name",
                                      "transcript_biotype"), 
                       values = age_abs_df$protein_coding_gene, mart = ensembl )


#length(unique(results_end_1$external_gene_name))

#merge and extract the protein coding transcripts

prots_of_int <- merge(age_abs_df, results_end_1, by.x = "protein_coding_gene", by.y = "external_transcript_name")

length(unique(prots_of_int$external_gene_name))
plot(prots_of_int$p_value)

plot(prots_of_int$correlation_coefficient)




#browseVignettes(AnnotationHub())

#explore the grch38 table loaded by the annotable library
ego_df <- enrichGO(gene = prots_of_int$ensembl_gene_id,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = T)

## Output results from GO analysis to a table
cluster_summary <- data.frame(ego_df)

dotplot(ego_df, showCategory = 20,
        
        font.size = 8) 


# pairwise_termsim(ego_df)
# emapplot(ego_df)
