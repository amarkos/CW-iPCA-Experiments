load("./TEP_results/perf_n_2500_chunks_5.RData")
load("./TEP_results/perf_n_10000_chunks_20.RData")
load("./TEP_results/perf_n_7500_chunks_15.RData")
load("./TEP_results/perf_n_5000_chunks_10.RData")
load("./TEP_results/perf_n_12500_chunks_25.RData")

library("tidyverse")
overall_perf = full_join(perf_n_2500_chunks_5,perf_n_5000_chunks_10)
overall_perf = full_join(overall_perf, perf_n_7500_chunks_15)
overall_perf = full_join(overall_perf, perf_n_10000_chunks_20)
overall_perf = full_join(overall_perf, perf_n_12500_chunks_25)

error_plot = overall_perf %>% 
  filter(measure=="error") %>% 
  mutate(chunks=parse_factor(as.character(chunks)),
         replicate=parse_factor(as.character(replicate))) %>% 
  ggplot(aes(x=chunks, y=value,fill=method))+
  geom_boxplot()+facet_grid(~method) + ylab("imputation error") + theme_bw()+
  theme(legend.position="none")+ theme(axis.text.y=element_text(size=10))

ggsave(file="TEP_imputation_error.pdf",
       error_plot,width=5,height=5)


obs_RV = overall_perf %>% 
  filter(measure=="obs_RV") %>% 
  mutate(chunks=parse_factor(as.character(chunks)),
         replicate=parse_factor(as.character(replicate))) %>% 
  ggplot(aes(x=chunks, y=value,fill=method))+
  geom_boxplot()+facet_grid(~method) + ylab("objects RV")+ theme_bw()+
  theme(legend.position="none")+ theme(axis.text.y=element_text(size=10))


ggsave(file="TEP_object_scores_RV.pdf",obs_RV,width=5,height=5)

att_RV = overall_perf %>% 
  filter(measure=="att_RV") %>% 
  mutate(chunks=parse_factor(as.character(chunks)),
         replicate=parse_factor(as.character(replicate))) %>% 
  ggplot(aes(x=chunks, y=value,fill=method))+
  geom_boxplot()+facet_grid(~method) + ylab("attributes RV") + theme_bw()+
  theme(legend.position="none")+ theme(axis.text.y=element_text(size=10))


ggsave(file="TEP_attribute_scores_RV.pdf",att_RV,width=5,height=5)
