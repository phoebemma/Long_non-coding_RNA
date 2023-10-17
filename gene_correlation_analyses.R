source("./libraries.R")
source("./R/Trainome_functions.R")

all_cpm <- data.table(readRDS("./data/cpm_values_all.RDS"), keep.rownames = T)%>% 
  #split the first column containing the gene 
  separate(rn, c("transcript_id", "transcript_name"), sep = "_", extra = "merge") %>%
  #drop transcript_id
  dplyr::select (-(transcript_id))

#Load the protein coding transcripts

protein_coding <- readr::read_csv("./data/protein_coding_transcripts.csv")

#load the filtered lncRNA data of significance

lncRNA_data <- readr::read_csv("./data/PreExc_lncRNA.csv")

#since this is just the preexercise data, we have to extract the intersect of both
#that is column names present in both dataframes

# Extract sample_ids that are common between lncRNAand protein_coding
corr_intersect <- (intersect(colnames(lncRNA_data), colnames(protein_coding)))

#extract the intersect from the protein coding data
protein_df <- protein_coding %>%
  subset( select = c("transcript_name", corr_intersect))
 
#There appears to be a duplicate column, remove it 
protein_df <- protein_df[, !duplicated(colnames(protein_df))]%>%
  #Round counts to one significant figure
  mutate_at(2:146, ~ as.integer(round(., 0))) 


#Load the dataframe of interest. That is those considered significant 
#as regards muscle mass change at baseline

mm_muscle_filt1 <- readr::read_csv("./lncRNA_relevant_for_mm.csv")
#extract the transcripts here that are in the lnc data
#lncRNA_intersect <- intersect(lncRNA_data$transcript_name, mm_muscle_filt1$target)

lnc_of_interest <- merge(lncRNA_data, mm_muscle_filt1, by.x = "transcript_name", by.y = "target", all.y = F) %>%
  #remove all the extra columns from x
  dplyr::select (-c(coef, Estimate,Std..Error, z.value, Pr...z..,
                    adj.p, log2fc, fcthreshold      ))
colnames(lnc_of_interest)


# Extract just the numeric data into a matrix with named rows by gene
rownames(lnc_of_interest) <- lnc_of_interest$transcript_name
geneExp_matrix <- as.matrix(lnc_of_interest[2:146])
head(geneExp_matrix)
heatmap(geneExp_matrix,
        Rowv=NA, Colv=NA)


##check correlation of the genes of interest to the protein_coding_ones
#This will only work when transformed
ress <- rcorr(as.matrix(t(protein_df[,-1]), t(lnc_of_interest[,-1])))

ress$r

# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values

head(ress$P)
