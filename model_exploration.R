source("./libraries.R")
source("./R/Trainome_functions.R")

#Load the model
model <- readRDS("./data/lncRNA_model.RDS")


#transform the model into a dataframe
model_df <- data.frame( bind_rows(model$model_summarises) %>%
                             mutate(target = rep(names(model$model_summarises), each = 9))) %>%
  #extract the adjusted p value
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2), 
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) #%>%
# ~ 1.4 fold change









#subset to include only those from or below p.value of 0.05 and significant threshold
#model_df_filtered <- filter(model_df,  adj.p <= 0.05 & fcthreshold == "s")%>%
  #remove the rows containing intercept
 # subset(!coef == "(Intercept)")


#get the ensemble gene IDs as gsea doesnt seem to work with transcript names
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
#extract transcript biotypes and transcript names
results_end_1 <- getBM(attributes = c("ensembl_gene_id",
                                      "ensembl_transcript_id","external_gene_name", "external_transcript_name"), 
                       values = model_df$target, mart = ensembl )

targets_updated <-  merge(model_df, results_end_1, by.x = "target", by.y = "external_transcript_name")%>%
  subset(!coef == "(Intercept)")



#save the updated model information
saveRDS(targets_updated, "./data/model_updated.RDS")


#subset to include only those from or below p.value of 0.05 and significant threshold
model_df_filtered <- filter(targets_updated,  adj.p <= 0.05 & fcthreshold == "s")
  #remove the rows containing intercept
  

#save the filtered model

#saveRDS(model_df_filtered, "./data/filtered_model.RDS")



#select the upregulated genes 
upreg_genes <- model_df_filtered %>%
  filter(Estimate > 0) %>%
  select(ensembl_gene_id)

#select the down_regulated genes
downreg_genes <- model_df_filtered %>%
  filter(Estimate < 0) %>%
  select(ensembl_gene_id)

#ont= Bp stands for biological function
#gsea_up <- enrichGO(targets_updated$ensembl_gene_id, OrgDb = org.Hs.eg.db, 
#                keyType = "ENSEMBL", ont = "BP",
#                universe = upreg_genes$ensembl_gene_id)
gsea_down <- enrichGO(targets_updated$ensembl_gene_id, OrgDb = org.Hs.eg.db, 
                      keyType = "ENSEMBL", ont = "BP",
                      universe = downreg_genes$ensembl_gene_id)

gsea_up <- simplify(gsea_up, by = )

dsea_up_df <- data.frame(gsea_up)



background <- bitr(targets_updated$ensembl_gene_id, fromType = "ENSEMBL",
                   toType =  c( "ENTREZID", "SYMBOL"),
                   OrgDb = org.Hs.eg.db)


entrez_id_up <- bitr(upreg_genes[, 1], fromType = "ENSEMBL",
                     toType = c( "ENTREZID", "SYMBOL"),
                     OrgDb = org.Hs.eg.db)

#Biological process
bp_up <- enrichGO(entrez_id_up[, 2], OrgDb = 'org.Hs.eg.db', 
                  ont = "BP", 
                  universe = background[,2], 
                  readable = TRUE)

#cellular component
cc_up <- enrichGO(entrez_id_up[, 2], OrgDb = 'org.Hs.eg.db', 
                  ont = "CC", 
                  universe = background[,2], 
                  readable = TRUE)


#molecular function

mf_up <- enrichGO(entrez_id_up[, 2], OrgDb = 'org.Hs.eg.db', 
                  ont = "MF", 
                  universe = background[,2], 
                  readable = TRUE)

bp_up2 <- simplify(bp_up)

cc_up2 <- simplify(cc_up, cutoff = 0.7, by = "p.adjust", select_fun = min)
mf_up2 <- simplify(mf_up, cutoff = 0.7, by = "p.adjust", select_fun = min)


simpl.df.up <- data.frame(bp_up) %>%
  mutate(ont = "BP") %>%
  rbind(data.frame(cc_up) %>%
          mutate(ont = "CC")) %>%
  rbind(data.frame(mf_up) %>%
          mutate(ont = "MF")) %>%
  mutate(change = "up")








