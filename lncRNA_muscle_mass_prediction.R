source("./libraries.R")
source("./R/Trainome_functions.R")

#This dataset cntains only pre-exercise samples and the percentage change in muscle mass post exercise
meta_df <- readr::read_csv("./data/all_metadata_with_pct_change.csv") %>%
  dplyr::select("participant", "sample_id", "study",
                "sex", "age", "condition", "pct_change")%>%
  mutate(sex = factor(sex, levels = c("male", "female"))) %>%
  
  print()


#Load the lncRNA data
lncRNA  <- readr::read_csv("./data/lncRNAs.csv")

# Extract sample_ids that are common between lncRNAand metadata
lncRNA_intersect <- (intersect(colnames(lncRNA), meta_df$sample_id))


#subsett the lncRNA data to only include the intersects
lncRNA_data <- lncRNA %>%
  subset( select = c("transcript_name", lncRNA_intersect))%>%
  #Round counts to one significant figure
  mutate_at(2:131, ~ as.integer(round(., 0))) %>%
  print()
colnames(lncRNA_data)

meta_df <- meta_df %>%
  filter(sample_id %in% c(lncRNA_intersect)) %>%
  print()

#How many samples of each study are represented in the data
table(meta_df$study)
table(meta_df$sex)
#table(meta_df$time)

hist(meta_df$pct_change)

args<- list(formula = y ~  age  + sex  + age:pct_change + pct_change +(1|participant),
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

#saveRDS(results, file = "./data/pct_change_model.RDS")

bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 5)) %>%
  filter(coef == "age") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues age in numbers")


bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 5)) %>%
  filter(coef == "pct_change") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues for pct_change")



bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 5)) %>%
  filter(coef == "age:pct_change") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues age on pct_change")
