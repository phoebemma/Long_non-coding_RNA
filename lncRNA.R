source("./libraries.R")
source("./R/Trainome_functions.R")
#install.packages('dplyr', type = 'source')

#Load the lncRNA data and the metadata


lncRNA <- readr::read_csv("./data/lncRNAs.csv") #%>%
 


# mutate_at(2:619, ~ as.integer(round(., 0))) %>%
 # print()

metadata <- readr::read_csv("./data/all_metadata.csv") %>%
  #drop the "MidExc rows 
  subset(time != "MidExc") %>%
  mutate(time = factor(time, levels = c("PreExc",    "PostExc")), 
         age = factor(age, levels = c("young", "old")), 
         sex = factor(sex, levels = c("male", "female"))) %>%
  
  print()
  
  

# Extract sample_ids that are common between lncRNAand metadata
lncRNA_intersect <- (intersect(colnames(lncRNA), metadata$sample_id))
  
  
#subsett the lncRNA data to only include the intersects
lncRNA_data <- lncRNA %>%
  subset( select = c("transcript_name", lncRNA_intersect))%>%
  #Round counts to one significant figure
  mutate_at(2:356, ~ as.integer(round(., 0))) %>%
  print()
colnames(lncRNA_data)
  

#extract genes with zero counts in all rows 
#non_zero_counts <- rowSums(lncRNA_data[, -1]) > 0

#percentage of non_zero counts 
#sum(non_zero_counts, na.rm = T)/length(non_zero_counts)

#lncRNA_data <- lncRNA_data[non_zero_counts,]


#select only metadata that intersect
meta_df <- metadata %>%
  filter(sample_id %in% c(lncRNA_intersect)) %>%
  print()

 


#How many samples of each study are represented in the data
table(meta_df$study)
table(meta_df$sex)
table(meta_df$time)

colnames(lncRNA_data)
#Extract all except the first column
transcript <- lncRNA_data %>%
  pivot_longer(names_to = "sample_id",
               values_to = "y",
               cols = X102PostExcVLL14:X98.subj40sample7) %>%
  #mutate(counts = as.integer(round(counts,0))) %>%
  print()




args<- list(formula = y ~  time  + age:time + (1|participant),
                       family = glmmTMB::nbinom2())
                     





results <- seq_wrapper(fitting_fun = glmmTMB::glmmTMB,
                       arguments = args,
                       data = lncRNA_data,
                       metadata = meta_df,
                       samplename = "sample_id",
                       summary_fun = sum_fun,
                       eval_fun = eval_mod,
                       exported = list(),
                       #additional_vars = NULL,
                       #subset = NULL,
                       cores = ncores)

summary(results) 

#save the result
#saveRDS(results, file = "./data/lncRNA_model.RDS")


results <- readRDS("./data/lncRNA_model.RDS")

#comparisons(results, variables = list(time = c("PreExc", "PostExc")),vcov = F,
            #newdata = datagrid(sex = "female")
           # )




model_eval_df <- bind_rows(results$model_evaluations) %>%
  mutate(target = names(results$model_evaluations))



model_sum_df <- bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 4))%>%
  subset(!coef == "(Intercept)") %>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"))




#merge the model summaries and model evaluation into one dataframe
model_results <- merge(model_sum_df, model_eval_df, by = "target")




#check those that pass model summary and evaluation characteristics
model_filtered <- model_results %>% filter(Pr...z.. <= 0.05 & fcthreshold == "s" & pval.disp >= 0.5 & pval.unif >= 0.5)%>%
  print()



 #saveRDS(model_filtered, "./data/filtered_all_transcripts.RDS")

plot( table(model_filtered$coef))


unique(model_filtered$target)




















                        

comparisons(results$model_fits[[1]], vcov = FALSE,
            newdata = "marginalmeans",
            variables = list(age = c("young", "old")), 
            by = "time",
            type = "response")


###Visualizing the datasets




bind_rows(results$model_evaluations) %>%
  mutate(target = names(results$model_evaluations)) %>%
  #print()
  
  #dplyr::filter(pval.unif < 0.05)%>%
  
  
  ggplot(aes(pval.unif)) + geom_histogram()
print()



#A pval.zinfl value < 1 means that the observed data has less zeros than expected,
#a value > 1 means that it has more zeros than expected
# 




# Instead of rejecting the null hypothesis based on a single p-value less than 5% significance level,
# this method compares an empirical probability distribution of p-values to a uniform distribution.
# Hence, this method entails a larger sample and naturally produces relatively reliable results.
# 

## Problem is, this is a simulated result

# bind_rows(results$model_summarises) %>%
#   mutate(target = rep(names(results$model_summarises), each = 9)) %>%
#   filter(coef == "timeMidExc") %>%
#   ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
#   ggtitle("Pvalues time(MidExc)")+
#   theme(axis.text = element_text(size = 15), text = element_text(size = 15),
#         plot.title = element_text(hjust = 0.5))


 


bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 4)) %>%
  filter(coef == "timePostExc") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("Pvalues time in full dataset")+
  theme(axis.text = element_text(size = 15), text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))





bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 4)) %>%
  filter(coef == "timePostExc:ageold") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("Pvalues interaction time with age (timePostExc:ageold)")+
  theme(axis.text = element_text(size = 15), text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))





bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 4)) %>%
  filter(coef == "timePreExc:ageold") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("Pvalues interaction time with age (timePreExc:ageold)")+
  theme(axis.text = element_text(size = 15), text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))






# 
# x <- results%>%
#   dplyr::select(gene, method, coef, estimate, p.val)


#Plot the values for time and age
all_coefs_df <- data.frame( bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 4))) %>%
  #filter(coef == "timePostExc:ageold")) %>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2), 
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
  #remove the rows containing intercept
  subset(!coef == "(Intercept)")%>%
  # ~ 1.4 fold change
  #subset to include only those from or below p.value of 0.05 and significant threshold
  filter( Pr...z.. <= 0.05 & fcthreshold == "s")


#check the distribution of the coefficients
table(all_coefs_df$coef)

#Extract the different coefficients into different dataframse

time_df <- all_coefs_df %>%
  filter(coef == "timePostExc")

#extract the interaction between preexercise and age
preExc_age <- all_coefs_df %>%
  filter(coef== "timePreExc:ageold")


postExc_age <- all_coefs_df %>%
  filter(coef == "timePostExc:ageold")



#save them in a subfolder of the data folder

#saveRDS(time_df, file = "./data/significant_coefs-from_model/significant_lncrnas_time.RDS")

#saveRDS(preExc_age, file = "./data/significant_coefs-from_model/significant_lncrnas_PreExc_age.RDS")

#saveRDS(postExc_age, file = "./data/significant_coefs-from_model/significant_lncrnas_PostExc_age.RDS")









#get the ensemble gene IDs as gsea doesnt seem to work with transcript names
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
#extract transcript biotypes and transcript names
results_end_1 <- getBM(attributes = c("ensembl_gene_id",
                                      "ensembl_transcript_id", "external_transcript_name"), 
                       values = time_vs_age$target, mart = ensembl )

targets_updated <-  merge(time_vs_age, results_end_1, by.x = "target", by.y = "external_transcript_name")





#subset to include only those from or below p.value of 0.05 and significant threshold
time_age_filtered <- filter(targets_updated,  adj.p <= 0.05 & fcthreshold == "s")%>%
  #remove the rows containing intercept
  subset(!coef == "(Intercept)")


#select the upregulated genes 
upreg_genes <- time_age_filtered %>%
  filter(Estimate > 0) %>%
  select(ensembl_gene_id)

#select the down_regulated genes
downreg_genes <- time_age_filtered %>%
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















ggplot(aes( Pr...z..)) + geom_histogram(bins = 80) +
   ggtitle(" lncRNA Pvalues time: age")




  # subset to include only those with values = or less than 0.05
time_vs_age_filtered <- subset(time_vs_age, Pr...z.. <= 0.05 & fcthreshold == "s")


#By time alone
 time_alone <- data.frame(bind_rows(results$model_summarises) %>%
   mutate(target = rep(names(results$model_summarises), each = 9)) %>%
   filter(coef == "timePostExc")) %>%
   subset(Pr...z.. <= 0.05)
 
 
   ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
   ggtitle("lncRNA Pvalues timePostExc")
 
 
 bind_rows(results$model_summarises) %>%
   mutate(target = rep(names(results$model_summarises), each = 9)) %>%
   filter(coef == "timeMidExc") %>%
   ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
   ggtitle("lncRNA Pvalues timeMidExc")
 
 #by age alone
 age_alone <- data.frame(bind_rows(results$model_summarises) %>%
   mutate(target = rep(names(results$model_summarises), each = 9)) %>%
   filter(coef == "ageold")) %>%
   subset(Pr...z.. <= 0.05)
   
   ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
   ggtitle("lncRNA Pvalues ageold")
 
 #by sex
 bind_rows(results$model_summarises) %>%
   mutate(target = rep(names(results$model_summarises), each = 9)) %>%
   filter(coef == "sexfemale") %>%
   ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
   ggtitle("lncRNA Pvalues sexfemale")
 
 bind_rows(results$model_summarises) %>%
   mutate(target = rep(names(results$model_summarises), each = 9)) %>%
   filter(coef == "timePostExc:sexfemale") %>%
   ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
   ggtitle("lncRNA Pvalues timePostExc:sexfemale")
