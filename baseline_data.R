source("./libraries.R")
source("./R/Trainome_functions.R")
#install.packages('dplyr', type = 'source')

#Load the lncRNA data and the metadata


lncRNA <- readr::read_csv("./data/lncRNAs.csv") #%>%



# mutate_at(2:619, ~ as.integer(round(., 0))) %>%
# print()

metadata <- readr::read_csv("./data/all_metadata.csv") %>%
   
  subset(time == "PreExc") %>%
  mutate(#time = factor(time, levels = c("PreExc", "MidExc",   "PostExc")), 
         age = factor(age, levels = c("young", "old")), 
         sex = factor(sex, levels = c("male", "female"))) %>%
  
  print()


# Extract sample_ids that are common between lncRNAand metadata
lncRNA_intersect <- (intersect(colnames(lncRNA), metadata$sample_id))


#subsett the lncRNA data to only include the intersects
lncRNA_data <- lncRNA %>%
  subset( select = c("transcript_name", lncRNA_intersect))%>%
  #Round counts to one significant figure
  mutate_at(2:158, ~ as.integer(round(., 0))) %>%
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
#table(meta_df$time)

colnames(lncRNA_data)
#Extract all except the first column
transcript <- lncRNA_data %>%
  pivot_longer(names_to = "sample_id",
               values_to = "y",
               cols = 2:158) %>%
  #mutate(counts = as.integer(round(counts,0))) %>%
  print()





args<- list(formula = y ~  age  + sex  +(1|participant),
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

#saveRDS(results, file = "./data/baseline_model.RDS")



bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 3)) %>%
  filter(coef == "ageold") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues age")+
  theme(axis.text = element_text(size = 15), text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))




bind_rows(results$model_summarises) %>%
  mutate(target = rep(names(results$model_summarises), each = 3)) %>%
  filter(coef == "sexfemale") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
  ggtitle("baseline lncRNA Pvalues sex")+
  theme(axis.text = element_text(size = 15), text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))





# bind_rows(results$model_summarises) %>%
#   mutate(target = rep(names(results$model_summarises), each = 4)) %>%
#   filter(coef == "ageold:sexfemale") %>%
#   ggplot(aes(Pr...z..)) + geom_histogram(bins = 80) +
#   ggtitle("baseline lncRNA Pvalues age on sex")









