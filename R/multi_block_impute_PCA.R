multi_block_impute_PCA <- function(Xlist, LabList=NULL,ncp = 4, scale = FALSE, pre_eg=NULL, method = "Regularized", threshold = 1e-06, seed = 1234, init = 1, maxiter = 10000, 
                                   row.w = NULL, coeff.ridge = 1, ...){
  source("R/add_es.R")
  source("R/do_es.r")
  source("R/block_impute_PCA.R")
  Xcom_list=list()
  true_imputed_vals=list()
  imputed_vals=list()
  imputed_labs=list()
  out=list()  
  # print(length(Xlist))
  Xcom=c()
  Xfull=c()
  for(k in 1:length(Xlist)){
    
    # print(k)
    
    X=as.matrix(Xlist[[k]])
    ########################
    ########################
    Xfill=X
    ########################
    ########################
    # print(LabList)
    if(is.null(LabList)==F){
      localLabs=LabList[[k]]
    }
    
    # rowna=unique(which(is.na(X),arr.ind = T)[,1]) ## the NA's are identified first
    
    #### IT's all there: the check is on rowna
    
    rowna=sort(unique(which(is.na(X),arr.ind = T)[,1]) )
    # print(length(rowna))
    if ((nrow(X)-length(rowna)) < (nrow(X)*.05)){
      # print("entered")
      rowna=1:nrow(X)
    }else{
      X_com = X[-rowna,]
    } # then the complete part is isolated from the one containing missing
    
    X_mis = X[rowna,]
    
    
    if(is.null(LabList)==F){
      
      local_row_Labs=localLabs[rowna,]
    }
    
    misX=which(is.na(X_mis))
    
    if(is.null(LabList)==F){
      imputed_labs[[k]]=local_row_Labs[misX]
    }
    
    # print(dim(X_mis))
    # print(nrow(X))
    # 
    if ((nrow(X)-length(rowna)) >= (nrow(X)*.05)){ ## if the number of rows of the complete part is larger than the 
      ## 5% of the block size
      block_eg = do_es(X_com) # decompose the complete part
      
      if(k!=1){
        # the complete part eigenspace will be merged with the eigenspace of the previous blocks
        pre_eg = add_es(pre_eg,block_eg,current_rank = ncp,method = "esm") 
      }else{
        # if it's the first block, the complete part eigenspace will be act as previous eigenspace
        pre_eg=block_eg
      }
      
    }else{## if the number of rows of the complete part is SMALLER than the 
      ## 5% of the block size
      # print("re-entered")
      if(k==1){
        pre_eg=NULL
      }
    }
    # print(paste("k:",k)) 
    # print(str(pre_eg))
    imp_out=block_impute_PCA(X_mis, ncp = ncp, scale = scale, pre_eg=pre_eg, method = method, threshold = threshold, seed = seed, init = init, 
                             maxiter = maxiter, row.w = row.w, coeff.ridge = coeff.ridge, ...)
    ########################
    ########################
    Xfill[rowna,]=imp_out$Xcom
    ########################
    ########################
    
    # print(head(imp_out$Xcom))
    imputed_vals[[k]]=imp_out$Xcom[misX]
    true_imputed_vals[[k]]=imp_out$imputed
    
    pre_eg=imp_out$eg 
    Xcom_list[[k]]=imp_out$Xcom
    Xcom=rbind(Xcom,imp_out$Xcom)
    ########################
    ########################
    Xfull=rbind(Xfull,Xfill)
    ########################
    ########################
  }
  out$Xcom_list=Xcom_list
  out$Xcom=Xcom
  out$Xfull=Xfull
  out$Xlist=Xlist
  # out$imputed_vals=imputed_vals
  out$imputed_vals=unlist(imputed_vals)
  out$imputed_vals_list=imputed_vals
  out$true_imputed_vals_list=true_imputed_vals
  if(is.null(LabList)==F){
    out$imputed_lab_list=imputed_labs
    out$imputed_labs=unlist(imputed_labs)
    true_impute_ord=sort(unlist(imputed_labs),index.return=T)$ix
    out$LabList=LabList
    out$imputed_vals=unlist(imputed_vals)[true_impute_ord]
  }
  
  return(out)
  
}
