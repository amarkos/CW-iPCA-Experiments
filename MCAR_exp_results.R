rm(list=ls())


if (!require("dplyr")) install.packages("dplyr")
if (!require("naniar")) install.packages("naniar")
if (!require("ggcorrplot")) install.packages("ggcorrplot")
if (!require("missMDA")) install.packages("missMDA")
if (!require("factoextra")) install.packages("factoextra")
if (!require("tidyverse")) install.packages("tidyverse")

# results_names=c(paste0("R/MCAR_experiments_results/MCAR_replicate_",c(1:10),"_exp_n_500_p_9_nblocks_",5,".RData"),
#                 paste0("R/MCAR_experiments_results/MCAR_replicate_",c(1:10),"_exp_n_500_p_9_nblocks_",10,".RData"),
#                 paste0("R/MCAR_experiments_results/MCAR_replicate_",c(1:10),"_exp_n_500_p_9_nblocks_",15,".RData"),
#                 paste0("R/MCAR_experiments_results/MCAR_replicate_",c(1:10),"_exp_n_500_p_9_nblocks_",20,".RData"),
#                 paste0("R/MCAR_experiments_results/MCAR_replicate_",c(1:10),"_exp_n_500_p_9_nblocks_",25,".RData"))

results_names=c(paste0("./Experiment_MCAR_results/MCAR_replicate_",1:10,"_exp_n_500_p_9_nblocks_",5,".RData"),
                paste0("./Experiment_MCAR_results/MCAR_replicate_",1:10,"_exp_n_500_p_9_nblocks_",10,".RData"),
                paste0("./Experiment_MCAR_results/MCAR_replicate_",1:10,"_exp_n_500_p_9_nblocks_",15,".RData"),
                paste0("./Experiment_MCAR_results/MCAR_replicate_",1:10,"_exp_n_500_p_9_nblocks_",20,".RData"),
                paste0("./Experiment_MCAR_results/MCAR_replicate_",1:10,"_exp_n_500_p_9_nblocks_",25,".RData"))


results_short_names=c(paste0("MCAR_replicate_",1:10,"_exp_n_500_p_9_nblocks_",5),
                      paste0("MCAR_replicate_",1:10,"_exp_n_500_p_9_nblocks_",10),
                      paste0("MCAR_replicate_",1:10,"_exp_n_500_p_9_nblocks_",15),
                      paste0("MCAR_replicate_",1:10,"_exp_n_500_p_9_nblocks_",20),
                      paste0("MCAR_replicate_",1:10,"_exp_n_500_p_9_nblocks_",25))

# results_short_names=c(paste0("MCAR_replicate_",c(1:10),"_exp_n_500_p_9_nblocks_",5),
#                       paste0("MCAR_replicate_",c(1:10),"_exp_n_500_p_9_nblocks_",10),
#                       paste0("MCAR_replicate_",c(1:10),"_exp_n_500_p_9_nblocks_",15),
#                       paste0("MCAR_replicate_",c(1:10),"_exp_n_500_p_9_nblocks_",20),
#                       paste0("MCAR_replicate_",c(1:10),"_exp_n_500_p_9_nblocks_",25))



for(i in 1:length(results_names)){
  load(results_names[i])
  assign(results_short_names[i],experiment_out)
}


error_results_names=paste0(results_short_names,"$errors")
errors_list=list()
for(j in 1:length(results_names)){
  errors_list[[j]]= as.tibble(eval(parse(text=error_results_names[j])))%>%
    mutate(scenario=error_results_names[j])
}

errors_tbl=bind_rows(errors_list)

tidy_errors=errors_tbl%>%
  separate(scenario,into=c("replicate","blocks"),sep="_exp_")%>%
  separate(replicate,into=c("replicate","replicate_num"),sep="e_",convert=TRUE)%>%
  separate(blocks,into=c("block_size","block_num"),sep="_nblocks_",convert=TRUE)%>%
  mutate(block_num = fct_inorder(as.factor(str_remove(block_num,pattern="\\$errors"))))%>%
  select(-replicate)%>%
  gather(key="method",value="error",imputePCA,inc_imp_naive,inc_imp_non_naive,factor_key=TRUE)


tidy_errors = tidy_errors%>%
  mutate(method=factor(fct_recode(method,
                           `CW-iPCA` = "inc_imp_non_naive",
                           `iPCA` = "imputePCA", 
                           `naive_CW-iPCA` = "inc_imp_naive"),
                           levels = c("CW-iPCA","iPCA","naive_CW-iPCA")),
         chunks=block_num,
         replicate=parse_factor(as.character(replicate_num))
  )

  error_boxes=  ggplot(data=tidy_errors, aes(x=chunks, y=error,fill=method))+
  geom_boxplot()+facet_grid(~method) + ylab("imputation error")+ theme_bw()+
  theme(legend.position="none") + theme(axis.text.y=element_text(size=10))

ggsave("MCAR_errors.pdf",error_boxes,width=5,height=5)
# ggsave("R/MCAR_experiment_results_errors_boxes.png",error_boxes)


######################################################
########### pull the indRVs #########################
######################################################
indRV_results_names=paste0(results_short_names,"$indRVs")
indRVs_list=list()
for(j in 1:length(results_names)){
  indRVs_list[[j]]= as.tibble(eval(parse(text=indRV_results_names[j])))%>%
    mutate(scenario=indRV_results_names[j])
}


indRVs_tbl=bind_rows(indRVs_list)

tidy_indRVs=indRVs_tbl%>%
  separate(scenario,into=c("replicate","blocks"),sep="_exp_")%>%
  separate(replicate,into=c("replicate","replicate_num"),sep="e_",convert=TRUE)%>%
  separate(blocks,into=c("block_size","block_num"),sep="_nblocks_",convert=TRUE)%>%
  mutate(block_num = fct_inorder(as.factor(str_remove(block_num,pattern="\\$indRVs"))))%>%
  select(-replicate)%>%
  gather(key="method",value="indRV",imputePCA,inc_imp_naive,inc_imp_non_naive,factor_key=TRUE)

resulting_plot=ggplot(data=tidy_indRVs,aes(x=block_num,y=indRV,group=method,color=method))+
  ylim(.45,1.2)+
  geom_smooth(method="loess")
ggsave("MCAR_experiment_results_indRVs.pdf",resulting_plot)


######################################################
########### pull the varRVs #########################
######################################################
varRV_results_names=paste0(results_short_names,"$varRVs")
varRVs_list=list()
for(j in 1:length(results_names)){
  varRVs_list[[j]]= as.tibble(eval(parse(text=varRV_results_names[j])))%>%
    mutate(scenario=varRV_results_names[j])
}


varRVs_tbl=bind_rows(varRVs_list)

tidy_varRVs=varRVs_tbl%>%
  separate(scenario,into=c("replicate","blocks"),sep="_exp_")%>%
  separate(replicate,into=c("replicate","replicate_num"),sep="e_",convert=TRUE)%>%
  separate(blocks,into=c("block_size","block_num"),sep="_nblocks_",convert=TRUE)%>%
  mutate(block_num = fct_inorder(as.factor(str_remove(block_num,pattern="\\$varRVs"))))%>%
  select(-replicate)%>%
  gather(key="method",value="varRV",imputePCA,inc_imp_naive,inc_imp_non_naive,factor_key=TRUE)

resulting_plot=ggplot(data=tidy_varRVs,aes(x=block_num,y=varRV,group=method,color=method))+
  ylim(.45,1.2)+
  geom_smooth(method="loess")
ggsave("MCAR_experiment_results_varRVs.pdf",resulting_plot)


