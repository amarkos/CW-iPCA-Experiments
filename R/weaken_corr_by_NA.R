########################################################################################
########################################################################################
#### HOW TO GENERATE MISSINGS THAT DESTROY THE CORRELATION STRUCTURE WITHIN A BLOCK ####

# To destroy the correlation between two variables x and y I compute the products of the 
# deviations from the mean, and take out (as missings) the top few. By removing, either in
# x or y, the values that produced the largest products of deviations, the correlation 
# coefficient drops down. 
# In a set of correlated variables, say V1, V2, ...,V4, all the unordered pairs are 
# considered. The merge of the 'top' position for a single variable, say V1 with all 
# the others, is made.


########################################################################################
########################################################################################


weaken_corr_by_NA<-function(mat,perc){
  
  library("missMDA")
  library("FactoMineR")
  library("purrr")
  library("dplyr")
  library("visdat")
  library("tidyr")
  library("tidyselect")
  p=ncol(mat)
  # print(mat)
  # mat=as_tibble(outmis$X[,1:4]) ## SELECT THE CORRELATED VARIABLE FROM THE SET_1
  
  
  #### SELECT THE POSITIONS OF THE NA'S SO THAT THE 
  mat=mat%>%
    map(function(x) x-mean(x))%>% # compute the deviations from the mean for each variable
    as_data_frame()
  nms=names(mat)
  
  mat_TF=as_tibble(model.matrix(data=mat,~(.)^2)) %>%
    mutate(perc=perc)%>%
    select(contains(":"))%>%
    map(function(x) x>quantile(x,perc))%>%
    as_tibble()
  
  
  for(i in 1:(p-1)){
    NApos=mat_TF%>%
      select(contains(paste0(names(mat)[i],":")))%>%
      apply(.,1,any)
    
    mat[NApos,i]=NA
  }
  
  return(mat)
}


