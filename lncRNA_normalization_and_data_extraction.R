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
