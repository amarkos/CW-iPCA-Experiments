snap_lab_maker = function(ggsnap,bw=F){
  
  snap_build= ggplot_build(ggsnap)
  data_snap=snap_build$data[[1]]

  
  miss_vals = data_snap %>% 
    group_by(x,fill) %>%
    summarize(n())%>%
    filter(fill==ifelse(bw,"grey20","grey50")) %>%
    mutate(perc_miss = `n()`/12500)
  # %>%
  #   ungroup() %>% 
  #   select(x,perc_miss)
  
  print(miss_vals)
  
  new_xlabs = paste0("V",1:9,": ",round(miss_vals$perc_miss,4)*100,"%")
  
  return(new_xlabs)
}
