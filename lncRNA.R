source("./libraries.R")
source("./R/Trainome_functions.R")


#Load the lncRNA data and the metadata

metadata <- readr::read_csv("./data/all_metadata.csv")
lncRNA <- readr::read_csv("./data/lncRNAs.csv")


colnames(lncRNA)

#Extract all except the first column
transcript <- lncRNA %>%
  pivot_longer(names_to = "sample_id",
               values_to = "counts",
               cols = X1023WVLL11:X99.subj40sample3) %>%
  mutate(counts = as.integer(round(counts,0))) %>%
  print()



form <- list(formula = list(counts ~  time + time:condition + age+ age:condition + sex:condition(1|participant)),
                       family = list(glmmTMB::nbinom2()))
                     


results <- seq_wrapper(fitting_fun = glmmTMB::glmmTMB,
                       arguments = form,
                       data = transcript,
                       metadata = metadata,
                       samplename = "sample_id",
                       cores = "max")
