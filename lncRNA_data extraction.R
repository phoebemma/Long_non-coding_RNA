source("./libraries.R")
source("./R/Trainome_functions.R")

library(trainomeHelper)
library(glmmTMB)

all_cpm <- data.table(readRDS("./data/cpm_values_all.RDS"), keep.rownames = T)%>% 
  #split the first column containing the gene 
  separate(rn, c("transcript_id", "transcript_name"), sep = "_", extra = "merge") %>%
  #drop transcript_id
  select (-(transcript_id))

#loadinf ensembl attributes as I would be using them to annotate the genes in the dataset.
#I am picking the gene_id, gene_name and gene biotypes and matching them to the listed genes in my dataset
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
#extract transcript biotypes and transcript names
results_end_1 <- getBM(attributes = c("transcript_biotype", "external_transcript_name"), values = all_cpm$transcript_name, mart = ensembl )
all_cpm_with_annotation <- merge(all_cpm, results_end_1, by.x = "transcript_name", by.y = "external_transcript_name")

#list all the biotypes in the dataset
unique(all_cpm_with_annotation$transcript_biotype)


#Extracting from the dataframe only the genes with biotype "lncRNA"
lncRNA <- filter(all_cpm_with_annotation, transcript_biotype == "lncRNA")


#check if only lncRNAs are contained in the dataframe
unique(lncRNA$transcript_biotype)

#subset the dataframe to exclude the biotype COLUMN
lncRNA <- subset(lncRNA, select = -(transcript_biotype))

#write_csv(lncRNA, "./data/lncRNAs.csv")

# Filter zero count rows
#lncRNA <- lncRNA[rowSums(lncRNA[-1,]) > 0]

#Load the differet metadata

#Contratrain metadata
ct_metadata <- readr::read_csv("./data/contratrain_metadata.csv") %>%
  select(sample_id,participant, sex, condition, time, study, age)

#COPD metadata
copd_metadata <- readr::read_csv("./data/copd_metadata.csv") %>%
  select(sample_id,participant, sex, condition, time, study, age)


#Volume meatadata

volu_metadata <- readr::read_csv("./data/volume_metadata.csv") %>%
  select(sample_id,participant, sex, condition, time, study, age)

metadata <- rbind(ct_metadata, volu_metadata, copd_metadata)

#write_csv(metadata, "./data/all_metadata.csv")




metadata <- readr::read_csv("./data/all_metadata.csv")




colnames(lncRNA)

#Extract all except the first column
transcript <- lncRNA %>%
  pivot_longer(names_to = "sample_id",
               values_to = "counts",
               cols = X1023WVLL11:X99.subj40sample3) %>%
  mutate(counts = as.integer(round(counts,0))) %>%
  print()










results <- seq_wrapper(fitting_fun = glmmTMB::glmmTMB,
                      arguments = counts ~ time + time:condition + sex + condition + study + (1|participant),
                      data = lncRNA,
                       metadata = metadata,
                     samplename = "sample_id",
                      cores = "max")










