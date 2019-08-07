rm(list = ls())


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

load("R/TEP_FaultFree_Training.RData")

tep_data_full=fault_free_training %>% 
  select(-faultNumber,-simulationRun, -sample) %>%
  scale()

tep_data_full=as.data.frame(tep_data_full)
n = nrow(tep_data_full)

#################################################################
#  PRE-PROCESSING ###############################################
#################################################################
tep_data_full_NA = tep_data_full


prop.m=.20
for (j in 1:ncol(tep_data_full_NA)) {
tep_data_full_NA[sample(1:n,size = round(n*prop.m),replace = FALSE),j] = NA
}

# miss_var_summary(data.frame(tep_data_full_NA))

tep_data_m = sequential_corr_shrinker(tep_data_full)

# MORE CLEANING
semi_final_vars=names(tep_data_m)
sel_tep_data = tep_data_full %>% select(semi_final_vars)

corr_vals=correlate(sel_tep_data,use = "pairwise.complete.obs",quiet=T)

sel_corr_vals_long=stretch(shave(corr_vals), na.rm = T) %>% 
  mutate(abs_r=abs(r)) %>% 
  filter(abs_r >= .35) %>% 
  arrange(desc(abs_r))

final_vars = unique(c(sel_corr_vals_long$x,sel_corr_vals_long$y))

sel_tep_data = sel_tep_data %>% select(final_vars)
tep_data_m = tep_data_m %>% select(final_vars)
tep_data_full_NA = tep_data_full_NA %>% select(final_vars)
# tep_data_full_final = as.data.frame(tep_data_full) %>% select(final_vars)

pdf(file="./CARME2019/CS_incremental_PCA_with_missings/figures/TEP_miss.pdf")
vis_miss(tep_data_m,sort_miss=FALSE,warn_large_data = F)+
  theme(axis.text.x=element_text(angle=90))
dev.off()

length(final_vars)
pdf(file="CARME2019/CS_incremental_PCA_with_missings/figures/TEP_corr.pdf")
corrplot(cor(sel_tep_data,use="pairwise.complete.obs"))
dev.off()

pdf(file="CARME2019/CS_incremental_PCA_with_missings/figures/TEP_weak_corr.pdf")
corrplot(cor(tep_data_m,use="pairwise.complete.obs"))
dev.off()

save(file="R/TEP_preprocessed_data.RDATA",sel_tep_data, tep_data_m, tep_data_full_NA)

  