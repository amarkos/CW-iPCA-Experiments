block_impute_PCA <- function(X, ncp = 4, scale = FALSE, pre_eg=NULL, method = "Regularized", threshold = 1e-06, seed = 1234, init = 1, maxiter = 10000, 
                             row.w = NULL, coeff.ridge = 1,prev_missing=NULL, ...){
  source("R/add_es.R")
  source("R/do_es.r")
  
  
  nb.iter <- 1
  old <- Inf
  objective <- 0
  # print(seed)
  set.seed(seed)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  
  
  X <- as.matrix(X)
  rowna=sort(unique(which(is.na(X),arr.ind = T)[,1]))
  # print(rowna)
  # row.w=1/nrow(X)
  
  ncp <- min(ncp, ncol(X), nrow(X) - 1) # number of components
  
  
  miss_pos=which(is.na(X),arr.ind=T)
  #############################################
  #############################################
  #### if there are previously analysed blocks, 
  #### then the position of the NAs in the previous blocks are taken into account 
  if(is.null(prev_missing)==F){
    # print("in")
    miss_pos=rbind(prev_missing,miss_pos)
  } 
  # print(missing)
  nmis=nrow(X)
  pmis=ncol(X)
  
  missing=sort(miss_pos[,1]+(miss_pos[,2]-1)*nmis, decreasing=F)
  # missing <- which(is.na(X)) # identifies the missings  
  # print(cbind(missing,(missing2)))
  # print((missing2))
  # print(unique(missing2))
  ### position of the missing to be returned
  #############################################
  #############################################
  
  mean.p <- apply(X, 2,function (x) mean(x,na.rm=T))
  
  Xhat <- t(t(X) - mean.p) # centers the matrix; Xhat is the centered matrix
  # Xhat = X
  
  et <- apply(Xhat, 2, sd) 
  
  
  
  if (any(is.na(X))){ # if X contains any missing, it sets them to 0 (they will be imputed with some random values)
    
    Xhat[missing] <- 0 # e
  }
  
  #############################################
  #############################################
  #############################################
  #############################################    
  #############################################
  # TOTAL RANDOM IMPUTATION
  init=2
  #############################################    
  #############################################    
  #############################################
  #############################################    
  
  if (init > 1){ # initialization is greater than one, it initializes the missing entries with random normal numbers
    
    Xhat=X
    Xhat[missing] <- rnorm(length(missing))
    
  }
  # print(Xhat[200:264,])
  
  fittedX <- Xhat # then it initialize the fitted version of X with the centered version of X, with random completion
  
  if (ncp == 0) # if the number of components is set to 0, then so is the nb.iter (what's this, and why it goes to 0?)
    nb.iter = 0
  iteration_count=0
  
  
  
  while (nb.iter > 0) {# if nb.iter is greater than 0, then the fitted 
    iteration_count=iteration_count+1
    
    Xhat[missing] <- fittedX[missing]
    
    
    
    #### un-transform the data ########################################
    # if (scale) 
    #   Xhat = t(t(Xhat) * et) # you de-scale if you scaled before
    if(iteration_count>1){
      Xhat <- t(t(Xhat) + mean.p) # you re-add the mean to the data
    }
    
    
    ########################################
    
    ### update the means and the sds according to the new imputations and re-transform the data
    ########################################################################
    ########################################################################
    ########################################################################
    ########## REMOVE COMMENTS TO MAKE IT WORK WITH NO STARTING EG ################
    if(is.null(pre_eg) ){
      # print("in")
      big_eg=do_es(Xhat)
      # pre_eg=big_eg
    }else{
      # if((iteration_count==1) && (is.null(pre_eg)==FALSE) ){  
      eg=do_es(Xhat)
      com_rank=min(length(eg$d),length(pre_eg$d))
      big_eg = add_es(pre_eg,eg,current_rank=com_rank,method="esm")
    }
    ########################################################################
    ########################################################################
    ########################################################################
    
    # if((iteration_count>1) && (is.null(pre_eg)==FALSE)){  
    #   
    #   eg=do_es(Xhat)
    #   com_rank=min(length(eg$d),length(pre_eg$d))
    #   big_eg = add_es(pre_eg,eg,current_rank=com_rank,method="esm")
    # }
    
    #################################################################################
    
    svd.res=big_eg # this is the global eigenspace, meaning it is the eigenspace of all the data analysed insofar
    # if(iteration_count==1){
    #   print("nrow(big_eg$u)")
    #   print(nrow(big_eg$u))
    # }
    final_eg=big_eg
    com_eg=pre_eg
    
    svd.res$u=svd.res$u*sqrt(svd.res$m)
    
    svd.res$d=svd.res$d*sqrt(1/svd.res$m)
    
    mean.p=as.vector(svd.res$orgn)
    
    
    tot_rows=nrow(svd.res$u)
    
    act_rows=(tot_rows-nrow(Xhat)+1):tot_rows
    
    svd.res$u=svd.res$u[act_rows,]
    
    # if(iteration_count==27){
    # 
    # print("head(svd.res$u)")
    # print(head(svd.res$u[200:264,1:5]))
    # print("head(svd.res$v)")
    # print(head(svd.res$v[,1:5]))
    # print("head(svd.res$d)")
    # print(svd.res$d)
    # 
    # }
    # 
    
    if(length(svd.res$d)<=ncp){
      ncp=ncp-1
    }
    
    sigma2 <- svd.res$m * ncol(X)/min(ncol(X), svd.res$m - 1) * 
      sum((svd.res$d[-c(1:ncp)]^2)/((svd.res$m -  1) * 
                                      ncol(X) - (svd.res$m - 1) * ncp - ncol(X) * 
                                      ncp + ncp^2))
    
    # print(length(svd.res$d))
    # print(svd.res$d[-c(1:ncp)]^2)
    
    sigma2 <- min(sigma2 * coeff.ridge, svd.res$d[ncp + 1]^2)### coeff.ridge is the lambda value for the ridge regularization
    
    
    # if(iteration_count==116){
    #   
    #   print("sigma")
    #   print(sigma2)
    #   
    # }
    
    #### UPDATE 12/12/2017 : I just inserted brackets, should not change a thing
    if (method == "EM"){ ## when using em, there is no regularization
      sigma2 <- 0}
    
    # shrinked lambda
    lambda.shrinked = (svd.res$d[1:ncp]^2 - sigma2)/svd.res$d[1:ncp]
    
    
    ## This is the reconstruction formula to re-impute the missings in X
    row.w=1/svd.res$m
    
    fittedX = tcrossprod(t(t(svd.res$u[, 1:ncp, drop = FALSE] * row.w) * lambda.shrinked), svd.res$v[, 1:ncp, drop = FALSE])
    
    
    fittedX <- fittedX/row.w ## insert the weights, if needed.
    
    # print(fittedX[200:264,])
    # if(iteration_count==116){
    #   print("fittedX")
    #   print(fittedX[missing])
    #   #   
    #   #   # print("head(Xhat)")
    #   #   # print(head(Xhat))
    # }
    # 
    Xhat <- t(t(Xhat) - mean.p) #centered (by the old mean)
    
    # print(head(Xhat))
    # print(head(fittedX))
    
    diff <- Xhat - fittedX ## difference between the previous centered X and the newly imputed version
    
    mydiff=diff[missing]
    
    
    myobjective=sum(mydiff^2)/length(missing)
    
    mycriterion= myobjective
    
    
    diff[missing] <- 0 ## in the missings positions, a 0 is put.
    objective <- sum(diff^2 * row.w) ## squared differences between the previous and the following matrix
    criterion <- abs(1 - objective/old) # compute the improvement
    old <- objective # the old objective becomes the current one, and we're ready for a new iteration
    nb.iter <- nb.iter + 1
    
    # 
    ### condition to exit the while loop
    if (!is.nan(mycriterion)) {
      # if ((mycriterion < threshold) && (nb.iter > 5)) 
      #   nb.iter <- 0
      if ((mycriterion < threshold) && (nb.iter > 5)) 
        nb.iter <- 0
    }
    if (nb.iter > maxiter) {
      nb.iter <- 0
      warning(paste("Stopped after ", maxiter, " iterations"))
    }
  }# END WHILE LOOP
  
  
  ### de-transform the data
  # if (scale) 
  #   Xhat <- t(t(Xhat) * et) #de-scale
  # 
  
  Xhat <- t(t(Xhat) + mean.p)           ## de-center 
  
  
  completeObs <- X                      ## starting matrix
  completeObs[missing] <- Xhat[missing] ## insert in your matrix the imputed values
  
  ### de-transform the fittedX
  # if (scale) 
  #   fittedX <- t(t(fittedX) * et)
  # print("mean.p")
  # print(mean.p)
  
  fittedX <- t(t(fittedX) + mean.p)
  
  # print(fittedX[missing])
  ############### SLIGHTLY DIFFERENT FROM WHAT JJ DOES##############
  completeObs[missing] <- fittedX[missing]
  ##################################################################
  
  result <- list()
  result$Xcom <- completeObs
  result$miss_row <- rowna
  result$miss_pos <- miss_pos
  result$fittedX <- fittedX
  result$imputed=fittedX[missing]
  result$missing=missing
  result$iterations=iteration_count
  result$eg=final_eg
  result$com_eg=pre_eg
  
  return(result)
}

