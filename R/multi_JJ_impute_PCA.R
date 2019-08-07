multi_JJ_impute_PCA <- function(Xlist,LabList=NULL, ncp = 4, scale = FALSE, pre_eg=NULL, method = "Regularized", threshold = 1e-06, seed = 1234, init = 1, maxiter = 1000, 
                                row.w = NULL, coeff.ridge = 1, ...){
  source("R/add_es.R")
  source("R/do_es.r")
  source("R/block_impute_PCA.R")
  require("missMDA")
  Ximp_list=list()
  imputed_vals=list()
  imputed_labs=list()
  out=list()  
  Xcom=c()
  Ximputed=c()
  
  for(k in 1:length(Xlist)){
    
    # print(k)
    
    X=as.matrix(Xlist[[k]])
    
    # rowna=unique(which(is.na(X),arr.ind = T)[,1]) ## the NA's are identified first
    # 
    # X_com = X[-rowna,] # then the complete part is isolated from the one containing missing
    # X_mis = X[rowna,]
    #
    
    Ximp=rbind(Ximputed,X)
    misX=which(is.na(Ximp))
    
    if(is.null(LabList)==F){
      if(k<2){
        localLabs=LabList[[k]]
      }else{
        localLabs=c()
        for(j in 1:k){
          localLabs=rbind(localLabs,LabList[[j]])
        }
      }
      imputed_labs[[k]]=localLabs[misX]      
    }
    
    
    
    
    imp_out=imputePCA(Ximp, ncp = ncp, scale = scale, method = method, threshold = threshold, seed = seed, init = init, 
                      maxiter = maxiter, row.w = row.w, coeff.ridge = coeff.ridge, ...)
    
    
    
    imputed_vals[[k]]=imp_out$completeObs[misX]
    Ximp_list[[k]]=imp_out$completeObs
    Ximputed=imp_out$completeObs
    
    
  }
  
  out$Ximp_list=Ximp_list
  out$Ximputed=Ximputed
  out$imputed_vals_list=(imputed_vals)
  out$imputed_vals=unlist(imputed_vals)
  if(is.null(LabList)==F){
    out$imputed_lab_list=imputed_labs
    out$imputed_labs=unlist(imputed_labs)
    true_impute_ord=sort(unlist(imputed_labs),index.return=T)$ix
    out$LabList=LabList
    out$imputed_vals=unlist(imputed_vals)[true_impute_ord]
  }
  
  
  return(out)
  
  
  
}
