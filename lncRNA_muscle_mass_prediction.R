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

hist(contratrain_meta$mm_rank, ylab = "Number of samples",
     main = "Distribution of  change in muscle mass ranking in the Contratrain data",
     xlab = "Muscle mass ranking in Contratrain data")
hist(contratrain_meta$pct_change, ylab = "Number of samples",
     main = "Distribution of percentage change in the Contratrain data",
     xlab = "Percentage change in Contratrain data")


copd_meta <- readr::read_csv("./muscle_mass_rank_metadata/copd_baseline_with_mm.csv") %>%
  dplyr::select(participant, sample_id, study, age, time, sex, pct_change, mm_rank) 
#The "X" in front of the sample_id is missing. To match that on the transcript data we add it
copd_meta$sample_id <- sub("^", "X", copd_meta$sample_id)


hist(copd_meta$pct_change,ylab = "Number of samples",
     main = "Distribution of percentage change in the COPD data",
     xlab = "Percentage change in COPD data" )
hist(copd_meta$mm_rank, ylab = "Number of samples",
     main = "Distribution of  change in muscle mass ranking in the COPD data",
     xlab = "Muscle mass ranking in COPD data")


vol_meta <- readr::read_csv("./muscle_mass_rank_metadata/volume_baseline_with_mm.csv") %>%
  dplyr::select(participant, sample_id, study, age, time, sex, pct_change, mm_rank) 

#Rename the preexercise to PreExc to match the other metadataframes
vol_meta$time[vol_meta$time == "w0"] <- "PreExc"


hist(vol_meta$pct_change,  ylab = "Number of samples",
     main = "Distribution of percentage change in the volume data",
     xlab = "Percentage change in Volume data")
hist(vol_meta$mm_rank, ylab = "Number of samples",
     main = "Distribution of  change in muscle mass ranking in the volume data",
     xlab = "Muscle mass ranking in Volume data")



metadata <- rbind(copd_meta, contratrain_meta, vol_meta)

#meta_df$sample_id[!duplicated(meta_df$sample_id), ]

#Load the lncRNA data
lncRNA  <- readr::read_csv("./data/lncRNAs.csv")

# Extract sample_ids that are common between lncRNAand metadata
lncRNA_intersect <- (intersect(colnames(lncRNA), metadata$sample_id))


#subsett the lncRNA data to only include the intersects
lncRNA_data <- lncRNA %>%
  subset( select = c("transcript_name", lncRNA_intersect))%>%
  #Round counts to one significant figure
  mutate_at(2:147, ~ as.integer(round(., 0))) %>%
  print()
colnames(lncRNA_data)

meta_df <- metadata %>%
  filter(sample_id %in% c(lncRNA_intersect)) %>%
  print()

#How many samples of each study are represented in the data
table(meta_df$study)
table(meta_df$sex)
table(meta_df$time)

hist(meta_df$pct_change, col = "blue", border = "white", ylab = "Number of samples",
     main = "Distribution of percentage change across samples",
     xlab = "Percentage change across samples")
hist(meta_df$mm_rank, col = "blue", border = "white", ylab = "Number of samples",
     main = "Distribution of ranks across samples",
     xlab = "Ranks across samples (0-1)")

args<- list(formula = y ~  age  + sex  + age:pct_change + pct_change + mm_rank + age:mm_rank +(1|participant),
            family = glmmTMB::nbinom2())
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

#saveRDS(results, file = "./data/muscle_mass_model.RDS")

bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 7)) %>%
  filter(coef == "age") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues age in numbers")


bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 7)) %>%
  filter(coef == "pct_change") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues for pct_change")



bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 7)) %>%
  filter(coef == "age:pct_change") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues age on pct_change")


bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 7)) %>%
  filter(coef == "age:mm_rank") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues age on muscle mass rank")



bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 7)) %>%
  filter(coef == "mm_rank") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues muscle mass rank")
