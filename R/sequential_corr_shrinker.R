sequential_corr_shrinker  <- function(data){
  library(tidyverse)
  library("corrplot")
  library("corrr")
  library("skimr")
  source("R/corr_shrinker.R")
  
  corr_vals=correlate(data,use = "pairwise.complete.obs",quiet=T)
  corr_vals_long=stretch(shave(corr_vals), na.rm = T)
  
  # corr_vals_long %>% select(r) %>% summarize(range) %>% print()
  
  disc_vars=corr_vals_long%>%
    filter(abs(r)>.95) %>% 
    pull(x) %>% unique()
  
  data = data %>% select(-disc_vars)
  
  # corr_vals=correlate(data,use = "pairwise.complete.obs",quiet=T)
  # corr_vals_long=stretch(shave(corr_vals), na.rm = T)
  # print(mean(corr_vals_long$r))
  
  cr_data=data
  
  criterion  = 100
  new_criterion  = 1
  # for(i in 1:100){
  it=0
    while(new_criterion <  criterion | it <=10 ){
    corr_vals=correlate(data,use = "pairwise.complete.obs",quiet=T)
it=it+1
# print(it)
    corr_vals_long=stretch(shave(corr_vals), na.rm = T)
    # print(mean(corr_vals_long$r))
    
    
    disc_vars=corr_vals_long%>%
      filter(abs(r)>.95) %>% 
      pull(x) %>% unique()
    
    data = data %>% select(-disc_vars)
    # cr_data = cr_data %>% select(-disc_vars)
    
    corr_vals_long = corr_vals_long%>%
      unite("var_pair",x,y,sep="_vs_") %>%
      mutate(var_pair=parse_factor(var_pair)) %>% 
      mutate(var_pair=fct_reorder(var_pair, desc(abs(r))))%>%
      separate(var_pair,into = c("x1","x2"),sep="_vs_") %>% 
      arrange(desc(abs(r)))
    
    criterion = new_criterion
    new_criterion  = mean(abs(corr_vals_long$r))
    
    pair_sel=corr_vals_long %>% slice(1) %>% select(x1,x2)
    
    out=corr_shrinker(data %>% pull(pair_sel$x1),data %>% pull(pair_sel$x2))
    
    # print(pair_sel)
    # print(out$old_corr)
    # print(out$new_corr)
    
    data = data %>% 
      mutate(!!pair_sel$x1 := out$new_vars$x1,
             !!pair_sel$x2 := out$new_vars$x2)  
    
    cr_data = cr_data %>%
      mutate(!!pair_sel$x1 := out$new_vars$x1,
             !!pair_sel$x2 := out$new_vars$x2)
    
    miss_out = map_df(data,function(x) prop_miss(x)) %>% 
      gather(key="vars",value="missings") %>%
      filter(missings>.45) %>% 
      pull(vars)
    
    # print(miss_out)
    
    data = data %>% select(-miss_out)
    # cr_data = cr_data %>% select(-miss_out)
    # print(length(names(data)))
  }
  return(cr_data)
  # return(data)
}