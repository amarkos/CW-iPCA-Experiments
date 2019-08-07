exact_multi_block_impute_PCA <- function(Xlist,LabList=NULL, ncp = 4, scale = FALSE, pre_eg=NULL, method = "Regularized", threshold = 1e-06, seed = 1234, init = 1, maxiter = 1000, 
                                         row.w = NULL, coeff.ridge = 1, ...){
  source("R/add_es.R")
  source("R/do_es.r")
  source("R/block_impute_PCA.R")
  Xcom_list=list()
  imputed_vals=list()
  out=list()  
  # print(length(Xlist))
  Xcom=c()
  for(k in 1:length(Xlist)){
    
    # print(k)
    
    X=as.matrix(Xlist[[k]])
    
    
    # rowna=unique(which(is.na(X),arr.ind = T)[,1]) ## the NA's are identified first
    rowna=sort(unique(which(is.na(X),arr.ind = T)[,1]))
    
    X_non_mis = X[-rowna,] # then the complete part is isolated from the one containing missing
    block_eg=do_es(X_non_mis)
    
    # print(head(block_eg$u[,1:5]))
    
    X_mis = X[rowna,]
    misX=which(is.na(X_mis))
    
    if(k!=1){
      pre_eg=add_es(pre_eg,block_eg,current_rank = ncp,method = "esm")
      
      ####################################################################################################
      ####################################################################################################
      # HERE I HAVE TO APPEND THE PREVIOUSLY IMPUTED PART TO THE MISSING PART, SO FOR THE MISSING POSITION
      
      X_mis=rbind(Xcom,X_mis)
      # block_mis_pos=which(is.na(X_all_mis),arr.ind = T)
      pre_mis=imp_out$miss_pos
      # print(dim(pre_mis))
      ####################################################################################################
      ####################################################################################################
      
      
    }else{
      pre_eg=block_eg
      pre_mis=NULL
    }
    
    
    
    # print(dim(X_mis))
    # print(nrow(X))
    # 
    
    # print(dim(X_mis))
    
    imp_out=block_impute_PCA(X_mis, ncp = ncp, scale = scale, pre_eg=pre_eg, method = method, threshold = threshold, seed = seed, init = init, 
                             maxiter = maxiter, row.w = row.w, coeff.ridge = coeff.ridge,prev_missing = pre_mis , ...)
    
    # print(tail(imp_out$miss_pos))
    
    # imputed_vals[[k]]=imp_out$Xcom[imp_out$missing]
    ############################
    # pre_eg=imp_out$eg 
    pre_eg=imp_out$com_eg
    all_mis_pos=imp_out$mis_pos
    ############################
    
    Xcom_list[[k]]=imp_out$Xcom
    Xcom=imp_out$Xcom
  }
  out$Xcom_list=Xcom_list
  out$Xcom=Xcom
  out$imputed_vals=Xcom[imp_out$missing]
  return(out)
  
}
