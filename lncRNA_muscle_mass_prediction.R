source("./libraries.R")
source("./R/Trainome_functions.R")

#This dataset cntains only pre-exercise samples and the percentage change in muscle mass post exercise
#meta_df <- readr::read_csv("./data/all_metadata_with_pct_change.csv") %>%
  #dplyr::select("participant", "sample_id", "study",
#                "sex", "age", "condition", "pct_change", "mm_rank")%>%
#  mutate(sex = factor(sex, levels = c("male", "female"))) %>%
  
#  print()


contratrain_meta <- readr::read_csv("./muscle_mass_rank_metadata/contratrain_baseline_with_mm.csv") %>%
  dplyr::select(participant, sample_id, study,age, time, sex, pct_change, mm_rank) 
#Rename the preexercise to PreExc to match the other metadataframes
contratrain_meta$time[contratrain_meta$time == "t1"] <- "PreExc"

hist(contratrain_meta$mm_rank, ylab = "Number of samples",col = "lightblue", border = "black",
     main = "Distribution of  change in muscle mass ranking in the Contratrain data",
     xlab = "Muscle mass ranking in Contratrain data",
     cex = 1.5, cex.lab = 1.5, cex.sub = 1.5, cex.main = 1.5)
hist(contratrain_meta$pct_change, ylab = "Number of samples",col = "lightblue", border = "black",
     main = "Distribution of percentage change in the Contratrain data",
     xlab = "Percentage change in Contratrain data",
     cex = 1.5, cex.lab = 1.5, cex.sub = 1.5, cex.main = 1.5)


copd_meta <- readr::read_csv("./muscle_mass_rank_metadata/copd_baseline_with_mm.csv") %>%
  dplyr::select(participant, sample_id, study, age, time, sex, pct_change, mm_rank) 
#The "X" in front of the sample_id is missing. To match that on the transcript data we add it
copd_meta$sample_id <- sub("^", "X", copd_meta$sample_id)


hist(copd_meta$pct_change,ylab = "Number of samples", col = "lightblue", border = "black",
     main = "Distribution of percentage change in the COPD data",
     xlab = "Percentage change in COPD data" ,
     cex = 1.5, cex.lab = 1.5, cex.sub = 1.5, cex.main = 1.5)
hist(copd_meta$mm_rank, ylab = "Number of samples",col = "lightblue", border = "black",
     main = "Distribution of  change in muscle mass ranking in the COPD data",
     xlab = "Muscle mass ranking in COPD data",
     cex = 1.5, cex.lab = 1.5, cex.sub = 1.5, cex.main = 1.5)


vol_meta <- readr::read_csv("./muscle_mass_rank_metadata/volume_baseline_with_mm.csv") %>%
  dplyr::select(participant, sample_id, study, age, time, sex, pct_change, mm_rank) 

#Rename the preexercise to PreExc to match the other metadataframes
vol_meta$time[vol_meta$time == "w0"] <- "PreExc"


hist(vol_meta$pct_change,  ylab = "Number of samples",  col = "lightblue", border = "black",
     main = "Distribution of percentage change in the volume data",
     xlab = "Percentage change in Volume data",
     cex = 1.5, cex.lab = 1.5, cex.sub = 1.5, cex.main = 1.5)
hist(vol_meta$mm_rank, ylab = "Number of samples", col = "lightblue", border = "black",
     main = "Distribution of  change in muscle mass ranking in the volume data",
     xlab = "Muscle mass ranking in Volume data",
     cex = 1.5, cex.lab = 1.5, cex.sub = 1.5, cex.main = 1.5)



metadata <- rbind(copd_meta, contratrain_meta, vol_meta)%>%
  drop_na() %>%
  #create an age classification
  mutate(age_class = case_when(age <= 40 ~ "young",
                               age > 40 ~ "old"))

#meta_df$sample_id[!duplicated(meta_df$sample_id), ]

#Load the lncRNA data
lncRNA  <- readr::read_csv("./data/lncRNAs.csv")

# Extract sample_ids that are common between lncRNAand metadata
lncRNA_intersect <- (intersect(colnames(lncRNA), metadata$sample_id))


#subsett the lncRNA data to only include the intersects
lncRNA_data <- lncRNA %>%
  subset( select = c("transcript_name", lncRNA_intersect))%>%
  #Round counts to one significant figure
  mutate_at(2:146, ~ as.integer(round(., 0))) %>%
  print()
colnames(lncRNA_data)
#write_csv(lncRNA_data, "./data/PreExc_lncRNA.csv")

meta_df <- metadata %>%
  filter(sample_id %in% c(lncRNA_intersect)) %>%
  drop_na()
  print()

#How many samples of each study are represented in the data
table(meta_df$study)
table(meta_df$sex)
table(meta_df$age_class)

hist(meta_df$pct_change, col = "lightblue", border = "black", ylab = "Number of samples",
     main = "Distribution of percentage change across samples",
     xlab = "Percentage change across samples",
     cex = 1.5, cex.lab = 1.5, cex.sub = 1.5, cex.main = 1.5)
hist(meta_df$mm_rank, col = "lightblue", border = "black", ylab = "Number of samples",
     main = "Distribution of ranks across samples",
     xlab = "Ranks across samples (0-1)",
     cex = 1.5, cex.lab = 1.5, cex.sub = 1.5, cex.main = 1.5)

#plot scatterplot of the changes in muscle mass across studies
ggplot(meta_df, aes(x= study, y= mm_rank)) +
  #geom_line()+
  geom_point()

set.seed(12345)
args<- list(formula = y ~  age  + age_class + sex + age_class:sex + (1|participant),
            family = glmmTMB::nbinom2())

#args2 <- list(formula = y ~  age  + mm_rank + age:mm_rank + (1|participant),
#              family = glmmTMB::nbinom2())

#args3 <- list(formula = y ~  pct_change + (1|participant),
#              family = glmmTMB::nbinom2())
ncores <- parallel::detectCores()
results <- seq_wrapper(fitting_fun = glmmTMB::glmmTMB,
                       arguments = args,
                       data = lncRNA_data,
                       metadata = meta_df,
                       samplename = "sample_id",
                       summary_fun = sum_fun,
                       eval_fun = eval_mod,
                       exported = list(),
                       additional_vars = NULL,
                       subset = NULL,
                       cores = ncores)
# results2 <- seq_wrapper(fitting_fun = glmmTMB::glmmTMB,
#                        arguments = args2,
#                        data = lncRNA_data,
#                        metadata = meta_df,
#                        samplename = "sample_id",
#                        summary_fun = sum_fun,
#                        eval_fun = eval_mod,
#                        exported = list(),
#                        additional_vars = NULL,
#                        subset = NULL,
#                        cores = ncores)
# 
# results3 <- seq_wrapper(fitting_fun = glmmTMB::glmmTMB,
#                         arguments = args3,
#                         data = lncRNA_data,
#                         metadata = meta_df,
#                         samplename = "sample_id",
#                         summary_fun = sum_fun,
#                         eval_fun = eval_mod,
#                         exported = list(),
#                         additional_vars = NULL,
#                         subset = NULL,
#                         cores = ncores)

#saveRDS(results, file = "./data/muscle_mass_model_updated.RDS")
results$model_summarises

results$model_evaluations


comparisons(results$model_fits[[1]], vcov = F,
            newdata = "marginalmeans",
            #focal variable is age_class. I am interested in the cpmparison between old and young
            variables = list(age_class = c("old", "young")), 
            #average over sex
            by = "sex",
            type = "response", 
            #hypothesis = "b1 = b2"
            )

# 
 # bind_rows(comparisons(results$model_fits[], vcov = F,
 #                       newdata = "marginalmeans",
 #                       variables = list(age_class = c("old", "young")),
 #                       ))






bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 5)) %>%
  filter(coef == "age") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues age in numbers") +
   theme(axis.text = element_text(size = 15), text = element_text(size = 15))  


bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 5)) %>%
  filter(coef == "sexmale") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues for sex") +
  theme(axis.text = element_text(size = 15), text = element_text(size = 15)) 






bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 5)) %>%
  filter(coef == "age_classyoung") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues age as a class")


bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 5)) %>%
  filter(coef == "age_classyoung:sexmale") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues impact age as a class on sex ")+
  theme(axis.text = element_text(size = 15), text = element_text(size = 15)) 


# 
# bind_rows(results2$model_summarises) %>%
#   mutate(target = rep(names(results2$model_summarises), each = 4)) %>%
#   filter(coef == "mm_rank") %>%
#   ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
#   ggtitle("baseline lncRNA Pvalues muscle ??? mass rank") +
#   theme(axis.text = element_text(size = 15), text = element_text(size = 15)) 
# 
# 
# bind_rows(results2$model_summarises) %>%
#   mutate(target = rep(names(results2$model_summarises), each = 4)) %>%
#   filter(coef == "age") %>%
#   ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
#   ggtitle("baseline lncRNA Pvalues age in numbers model 2") +
#   theme(axis.text = element_text(size = 15), text = element_text(size = 15))
# 
# 
# 
# 
# bind_rows(results3$model_summarises) %>%
#   mutate(target = rep(names(results3$model_summarises), each = 2)) %>%
#   filter(coef == "pct_change") %>%
#   ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
#   ggtitle("baseline lncRNA Pvalues pct_alone") +
#   theme(axis.text = element_text(size = 15), text = element_text(size = 15))  
#Doing the fdr step
model_df <- data.frame( bind_rows(results$model_summarises) %>%
                          mutate(target = rep(names(results$model_summarises), each = 5))) %>%
  #extract the adjusted p value
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2), 
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
  subset(!coef == "(Intercept)") 
  
# mm_muscle <- model_df %>%
#   filter(fcthreshold == "mm_rank")


#filter only those considered significant based on p and threshold values
model_filt <- model_df %>%
  filter(Pr...z.. <= 0.05 & fcthreshold == "s")



#plotting based on the different coefficients
filt_age_class <- model_filt %>%
  filter(coef == "age_classyoung")
ggplot(filt_age_class, aes(Pr...z..)) +
  geom_histogram() +
  ggtitle("P values of statistically significant values of differentially expressed lncRNAs by age")
  



filt_age_num <- model_filt %>%
  filter(coef == "age")
ggplot(filt_age_num, aes(Pr...z..)) +
  geom_histogram() +
  ggtitle("P values of statistically significant values of differentially expressed lncRNAs by age in numbers")


table(model_filt$coef)
#none were considered significant for age in numbers



filt_age_sex <- model_filt %>%
  filter(coef == "age_classyoung:sexmale")
ggplot(filt_age_sex, aes(Pr...z..)) +
  geom_histogram() +
  ggtitle("P values of statistically significant values of differentially expressed lncRNAs by interaction age:sex")



filt_sex <- model_filt %>%
  filter(coef == "sexmale")
ggplot(filt_sex, aes(Pr...z..)) +
  geom_histogram() +
  ggtitle("P values of statistically significant values of differentially expressed lncRNAs by sex")







#write_csv(mm_muscle_filt1, "./lncRNA_relevant_for_mm.csv")
#list_ <- mm_muscle_filt1$target

#save as csv file
#capture.output(list_, file = "./muscle_transcript_list.csv")

#c <- readr::read_csv("./muscle_transcript_list.csv")
  
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
#extract transcript biotypes and transcript names from filt-age-class
results_end_1 <- getBM(attributes = c("ensembl_gene_id",
                                      "ensembl_transcript_id","external_gene_name", "external_transcript_name"), 
                       values = filt_age_class$target, mart = ensembl )

targets_updated <-  merge(filt_age_class, results_end_1, by.x = "target", by.y = "external_transcript_name")


### Correlation analyses

#load the protein_coding transcripts
protein_coding <- readr::read_csv("./data/protein_coding_transcripts.csv")


#extract the transcripts in targets_updated that exist in the lncRNA_DATA

age_class_intersect <- (intersect(filt_age_class$target, lncRNA_data$transcript_name))

age_lncRNA_data <- lncRNA_data %>%
  subset(select =c(transcript_name %in% age_class_intersect))


# Extract sample_ids that are common between lncRNAand protein_coding
corr_intersect <- (intersect(colnames(lncRNA_data), colnames(protein_coding)))

#extract the intersect from the protein coding data
protein_df <- protein_coding %>%
  subset( select = c("transcript_name", corr_intersect))

#There appears to be a duplicate column, remove it 
protein_df <- protein_df[, !duplicated(colnames(protein_df))]%>%
  #Round counts to one significant figure
  mutate_at(2:146, ~ as.integer(round(., 0))) 


res <- cor(t(lncRNA_data[, -1]), t(protein_df[, -1]))

res.head()

#list_ <- targets_updated$ensembl_transcript_id

#save as csv file
#capture.output(list_, file = "./muscle_transcript_list.txt")

#c <- readr::read_csv("./muscle_transcript_list.txt")



#Gene ontology using topGo and KEGGREST
#as described here https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html
 #geneList <- targets_updated$Pr...z..
 #names(geneList) <- targets_updated$ensembl_transcript_id
 # Create topGOData object
 #GOdata <- new("topGOdata",
 #              ontology = "BP",
 #              allGenes = geneList,
 #              geneSelectionFun = function(x)x,
 #              annot = annFUN.org, mapping = "org.Hs.eg.db")

 #listDatabases()
 #pull all pathways for humans
# pathways.list <- keggList("pathway", "hsa") 
# head(pathways.list)
# 
# # Pull all genes for each pathway
# pathway.codes <- sub("path:", "", names(pathways.list)) 
# genes.by.pathway <- sapply(pathway.codes,
#                            function(pwid){
#                              pw <- keggGet(pwid)
#                              if (is.null(pw[[1]]$GENE)) return(NA)
#                              pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
#                              pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
#                              return(pw2)
#                            }
# )
# tail(genes.by.pathway)
# 
# 
# #use the targets_updated file to create a gene list
# geneList <- targets_updated$Pr...z..
# names(geneList) <- targets_updated$external_gene_name
# head(geneList)
# 
# 
# #Apply Wilcoxon rank-sum test to each pathway, testing if "in" p-values are smaller than "out" p-values:
# 
# # Wilcoxon test for each pathway
# pVals.by.pathway <- t(sapply(names(genes.by.pathway),
#                              function(pathway) {
#                                pathway.genes <- genes.by.pathway[[pathway]]
#                                list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
#                                list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
#                                scores.in.pathway <- geneList[list.genes.in.pathway]
#                                scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
#                                if (length(scores.in.pathway) > 0){
#                                  p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
#                                } else{
#                                  p.value <- NA
#                                }
#                                return(c(p.value = p.value, Annotated = length(list.genes.in.pathway)))
#                              }
# ))
# 
# 
# # Assemble output table
# outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
# outdat$pathway.name <- pathways.list[outdat$pathway.code]
# outdat$p.value <- pVals.by.pathway[,"p.value"]
# outdat$Annotated <- pVals.by.pathway[,"Annotated"]
# outdat <- outdat[order(outdat$p.value),]
# head(outdat)
# 
# 
# browseVignettes("DGCA")
# diffcor00