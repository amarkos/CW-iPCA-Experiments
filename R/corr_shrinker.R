corr_shrinker <- function(x1,x2,thresh=.75){
  library(tidyverse)
  
  out_vars=tibble(x1,x2) %>% 
    mutate(x12=x1*x2)
  
  r=as.numeric(correlate(tibble(x1,x2),use = "pairwise.complete.obs",quiet=T)[1,3]) #cause a tibble is created
                                                      # and the first col is for
                                                      # rownames, so it's equivalent
                                                      # to selecting the position 1,2
  # print(r)
  
  if(r>=.3){
    out_vars=out_vars %>% 
      mutate(wNA=case_when(x12<=quantile(x12,thresh,na.rm=T)~ x12,
                           TRUE ~  NA_real_))
  }else if(r<=-.3){
    out_vars=out_vars %>% 
      mutate(wNA=case_when(x12>quantile(x12,1-thresh,na.rm=T)~ x12,
                           TRUE ~  NA_real_))
  }else{
    
    out_vars=out_vars %>% 
    mutate(wNA=x12)
  }
  
  out_vars = out_vars %>% 
    mutate(uni=runif(n()),
           x1 = case_when(
             is.na(wNA) & uni >=.5 ~ NA_real_,
            TRUE ~ x1),
           x2 = case_when(
             is.na(wNA) & uni <.5 ~ NA_real_,
             TRUE ~ x2),
           )
  
  out=list()
  out$new_vars=out_vars
  out$old_corr = r
  new_corr = out_vars %>% filter(!is.na(wNA)) %>% select(x1,x2) %>% correlate(use="pairwise.complete.obs",quiet=T) 
  out$new_corr=as.numeric(new_corr[1,3])
  return(out)
}