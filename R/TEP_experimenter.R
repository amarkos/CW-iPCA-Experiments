# rm(list = ls())
TEP_experimenter <- function(n_sam = 5000, nblocks=10, myseed=1234, ncp = 5, rep=0){
  
  if (!require("dplyr")) install.packages("dplyr")
  if (!require("naniar")) install.packages("naniar")
  if (!require("ggcorrplot")) install.packages("ggcorrplot")
  if (!require("missMDA")) install.packages("missMDA")
  if (!require("factoextra")) install.packages("factoextra")
  if (!require("tidyverse")) install.packages("tidyverse")
  if (!require("FactoMineR")) install.packages("FactoMineR")
  if (!require("corrr")) install.packages("corrr")
  if (!require("skimr")) install.packages("skimr")
  
  source("R/sequential_corr_shrinker.R")
  source("R/weaken_corr_by_NA.R")
  source("R/Impute_Experiment_multi_block.R")
  
  load("TEP_results/TEP_preprocessed_data.RDATA")
  
  
  
  
  n_samd = n_sam/nblocks
  
  myseed = myseed + rep
  
  set.seed(myseed)
  sam_sel_tep_data = sel_tep_data %>% sample_n(n_sam)
  set.seed(myseed)
  sam_tep_data_m = tep_data_m %>% sample_n(n_sam)
  set.seed(myseed)
  sam_tep_data_full_NA = tep_data_full_NA %>% sample_n(n_sam)
  
  n=nrow(sam_tep_data_m)
  p=ncol(sam_tep_data_m)
  
  
  
  tep_data_MCAR_MNAR = sam_tep_data_m
  tep_data_MCAR_MNAR[1:n_samd, ] = sam_tep_data_full_NA[1:n_samd, ]
  
  
  
  outmis=list(X=as.matrix(sam_sel_tep_data),Xmis=as.matrix(tep_data_MCAR_MNAR))#,ncp=ncp)
  MNCAR_tep_out =  Impute_Experiment_multi_block(matList = outmis,n = n, p = p, ncp = ncp, 
                                                 nblocks = nblocks, prop_mis = prop.m, myseed=myseed)
  
  
  
  performance = as_tibble(rbind(MNCAR_tep_out$errors,
                                MNCAR_tep_out$indRVs,
                                MNCAR_tep_out$varRVs))
  names(performance)=c("iPCA", "CW-iPCA", "naive_CW-iPCA")
  
  performance = performance %>% gather(key="method",value="value" ,
                                       `iPCA`, `CW-iPCA`, `naive_CW-iPCA`) %>% 
    mutate(measure=rep(c("error","obs_RV","att_RV"),3),
           n_obs = n_sam,
           chunks = nblocks,
           ncp=ncp,
           replicate = rep) %>% 
    select(method, replicate, n_obs,chunks,ncp, measure, value)
  
  return(performance)
  
}



