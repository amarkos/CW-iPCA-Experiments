do_es <- function(data,init=1,mean_impute=NULL) {
 # require("corpcor")
  # data: data matrix
  out = list()
  
  data <- as.matrix(data)
  data_hat=data
  orgn = apply(data, 2, function(x)mean(x,na.rm=T))
  if (any(is.na(data))){
  w_mis <- which(is.na(data),arr.ind = T)
  mis_by_col=table(w_mis[,2])
  missings <- which(is.na(data),arr.ind=T)
  # data_hat=t(t(data_hat))-orgn
  data_hat[missings]=0
  }
  if (init < 1){ # initialization is greater than one, it initializes the missing entries with random normal numbers
    
    data_hat[missing] <- rnorm(length(missing))
  }
  
  # data_hat=t(t(data_hat))+orgn
  
  m = nrow(data_hat)  ## number of rows
  p = ncol(data_hat)  ## number of columns
  #if (is.cate == F) {
  
  
  orgn = apply(data, 2, function(x)mean(x,na.rm=T))
  
  orgn = t(t(orgn))
  # if (is.svd == T) {
  data_hat = t(data_hat)
  #  }
  oner = matrix(1, 1, m)
  
  
  #    
  
  cen.mat = data_hat - (orgn) %*% (oner)
#   print(t(cen.mat[1:5,1:6]))
 # svd.res = svd(cen.mat)
#  print(dim(svd.res$u))
  svd.res = fast.svd(cen.mat,0)
#  print(dim(svd.res$u))
  u = svd.res$u
  d = svd.res$d
  v = svd.res$v
  #  out = list()
  out$v = u
  #out$cen.mat = cen.mat
  #    }
  #out=list()
  out$m = m
  if (any(is.na(data))){
  out$true_m = m-mis_by_col
  }
  out$orgn = orgn
  out$u = v
  out$d = d
  out
} 
