add_es <- function(eg,eg2,current_rank,ff=0,method=c("esm","isvd")){
source("R/add_svd.R")
source("R/add_eig.R")
    if (missing("current_rank")) {
    #full rank
    current_rank =  2
  }
  if(method=="esm"){
    out = add_eig(eg, eg2,current_rank=current_rank)}
  else{
    if (is.null(eg$m)) {
      #without orgn data is assumed as zero-mean
      m = dim(eg$u)[1]
    } 
    else
    {
      m = eg$m
    }
    B = eg2
    out = add_svd(eg,B,m,current_rank=current_rank,ff = ff)
  }
  return(out)
}
