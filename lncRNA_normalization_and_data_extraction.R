source("./libraries.R")
source("./R/Trainome_functions.R")



#copd_2 <- extract_rsem_isoform_counts("./data/COPD_rsem_isoforms/COPD_rsem_isoforms/")


#rename by removing everything after the underscore'
#this would match the samplenames in the metadata
#colnames(copd_2) <- gsub("_.*", "", colnames(copd_2) )

#colnames(copd_2)[1] = "transcript_id"
#drop the duplicated column names. This is because COPD samples are in duplicates
#copd_2 <- copd_2[, !duplicated(colnames(copd_2))]

#write_csv(copd_2, "./data/copd_isoform_counts.csv")



copd_counts <-  readr::read_csv("./data/copd_isoform_counts.csv") 
#write_csv(copd_counts, "./data/copd_isoform_counts.csv")
contratrain_counts <- readr::read_csv("./data/contratrain_isoform_counts.csv") 

volume_counts <- readr::read_csv("./data/volume_isoform_counts.csv") 



metadata <- readr::read_csv("./data/all_metadata.csv")

#volume_ <- extract_rsem_isoform_counts("./data/Volume_rsem_isoforms/Volume_rsem_isoforms/")
#write_csv(volume_, "./data/volume_isoform_counts.csv")

#merge all the datasets based on trasncript_id
full_counts <- copd_counts %>%
  inner_join(contratrain_counts, join_by("transcript_id")) %>%
  inner_join(volume_counts, join_by("transcript_id")) %>%
  separate(transcript_id, c("transcript_id", "transcript_name"), sep = "_", extra = "merge") %>%
  #drop transcript_id
  dplyr::select(-(transcript_name))


#save dataset
#write_csv(full_counts, "./data/all_isoforms_unnormalized.csv")




#Use the ensemble database to get the annotation of the transcripts
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
#extract transcript biotypes and transcript names
results_end_1 <- getBM(attributes = c("transcript_biotype", "ensembl_transcript_id_version"), values = full_counts$transcript_id, mart = ensembl )
all_cpm_with_annotation <- merge(full_counts, results_end_1, by.x = "transcript_id", by.y = "ensembl_transcript_id_version")




#list all the biotypes in the dataset
unique(all_cpm_with_annotation$transcript_biotype)


#Extracting from the dataframe only the genes with biotype "lncRNA"
lncRNA <- filter(all_cpm_with_annotation, transcript_biotype == "lncRNA")


#check if only lncRNAs are contained in the dataframe
unique(lncRNA$transcript_biotype)

#subset the dataframe to exclude the biotype COLUMN
lncRNA <- subset(lncRNA, select = -(transcript_biotype))

#save dataset
#write_csv(lncRNA, "./data/full_lncRNAs.csv")
#saveRDS(lncRNA, "./data/full_lncRNAs.RDS")

#remove rows with all zeros
#lncRNA[-1] <- lncRNA[-1][rowSums(lncRNA[-1]) > 0]

genes <- lncRNA$transcript_id


