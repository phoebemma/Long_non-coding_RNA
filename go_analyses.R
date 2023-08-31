source("./libraries.R")
source("./R/Trainome_functions.R")

#Gene ontology enrichment analyses as shown on https://rpubs.com/Akhil_Velluva/GOfuncR'

#load the filtered model
model_df_filtered <- readRDS("./data/filtered_model.RDS")


#load the full model from which the filtered model was extracted
model_df <- readRDS("./data/model_updated.RDS")


unique(model_df_filtered$coef)


filtered_list <- model_df_filtered %>%
  dplyr::select(target, external_gene_name, coef)

unfiltered_list <- model_df %>%
  dplyr::select(target, external_gene_name, coef)


go_list <- unfiltered_list %>%
  mutate(is_candidate = ifelse(unfiltered_list$external_gene_name %in% filtered_list$external_gene_name, 1, 0))
  #model_df %>%
  #extract the significant genes as one, and none significant as zero
  #mutate(is_candidate = if_else(fcthreshold == "s", 1, 0)) %>%
  #select(c(ensembl_gene_id, is_candidate))

#extract the individual dataframes based on their coefficients
age <- go_list %>%
  subset(coef == "ageold")%>%
  dplyr::select(external_gene_name, is_candidate)

postExc_age <- go_list %>%
  subset(coef == "timePostExc:ageold")



go_enriched <- go_enrich(age)

#Get enrichment results
results <- go_enriched$results


#categorize enriched results by p value

Over_represented <- results[results$raw_p_overrep <= 0.05, ]
under_represented <- results[results$raw_p_underrep <= 0.05, ]


#Genes used in enrichment
genes_used <- go_enriched$genes
