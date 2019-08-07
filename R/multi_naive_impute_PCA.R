multi_naive_impute_PCA <- function(Xlist,LabList=NULL, ncp = 4, scale = FALSE, pre_eg=NULL, method = "Regularized", threshold = 1e-06, seed = 1234, init = 1, maxiter = 1000, 
                                   row.w = NULL, coeff.ridge = 1, ...){
  source("R/add_es.R")
  source("R/do_es.r")
  source("R/block_impute_PCA.R")
  Xcom_list=list()
  out=list()  
  Xcom=c()
  imputed_vals=list()
  imputed_labs=list()
  for(k in 1:length(Xlist)){
    
    # print(k)
    
    X=as.matrix(Xlist[[k]])
    misX=which(is.na(X))
    
    if(is.null(LabList)==F){
      localLabs=LabList[[k]]
      imputed_labs[[k]]=localLabs[misX]
    }
    
    # rowna=unique(which(is.na(X),arr.ind = T)[,1]) ## the NA's are identified first
    # 
    # X_com = X[-rowna,] # then the complete part is isolated from the one containing missing
    # X_mis = X[rowna,]
    
    
    imp_out=naive_impute_PCA(X, ncp = ncp, scale = scale, method = method, threshold = threshold, seed = seed, init = init, 
                             maxiter = maxiter, row.w = row.w, coeff.ridge = coeff.ridge, ...)
    
    #lambda <- c(1,2,5)
    #U <- qr.Q(qr(matrix(rnorm(30),10,3)))
    #x <- U %*% diag(sqrt(lambda)) %*% rnorm(3) + rnorm(10, sd =.05)
    #x.na <- x
    #x.na[c(1,3,7)] <- NA
    #x.imputed <- impute(lambda,U,x.na)
    #cbind(x,x.imputed)
    
    Xcom_list[[k]]=imp_out$Xcom
    Xcom=rbind(Xcom,imp_out$Xcom)
    imputed_vals[[k]]=imp_out$Xcom[misX]
  }
  
  out$Xcom_list=Xcom_list
  out$Xcom=Xcom
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
