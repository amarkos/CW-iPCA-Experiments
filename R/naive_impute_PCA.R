naive_impute_PCA<-function (X, ncp = 2, scale = TRUE, method = c("Regularized", 
                                                                  "EM"), row.w = NULL, coeff.ridge = 1, threshold = 1e-06, 
                             seed = NULL, nb.init = 1, maxiter = 1000, ...) 
{
  #the impute function does most of the job
  impute <- function(X, ncp = 4, scale = TRUE, method = NULL, 
                     threshold = 1e-06, seed = NULL, init = 1, maxiter = 1000, 
                     row.w = NULL, coeff.ridge = 1, ...) {
    
    #This function computes the mean of a variable, weighted by 'poids' and such that it ignores missings
    moy.p <- function(V, poids) {
      res <- sum(V * poids, na.rm = TRUE)/sum(poids[!is.na(V)])
    }
    
    #This function computes the sd of a centered variable, weighted by 'poids' and such that it ignores missings
    ec <- function(V, poids) {
      res <- sqrt(sum(V^2 * poids, na.rm = TRUE)/sum(poids[!is.na(V)]))
    }
    
    ## initialization stuff
    nb.iter <- 1
    old <- Inf
    myold <- Inf
    objective <- 0
    seed=1234
    set.seed(seed)
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    # print(seed)
    
    X <- as.matrix(X)
    
    ncp <- ncp #min(ncp, ncol(X), nrow(X) - 1) # number of components
    missing <- which(is.na(X)) # identifies the missings
    mean.p <- apply(X, 2, moy.p, row.w)# computes the weigthed mean, given the weights in row.w
    Xhat <- t(t(X) - mean.p) # centers the matrix; Xhat is the centered matrix
    et <- apply(Xhat, 2, ec, row.w) # computes the sd's, taking in input the centered Xhat
    # print("Xhat")
    # print(Xhat)
    if (scale) 
      Xhat <- t(t(Xhat)/et) # it scales Xhat, if required 
    
    if (any(is.na(X))) # if X contains any missing, it sets them to 0 (they will be imputed with some random values)
      Xhat[missing] <- 0 # e
    
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
    
    fittedX <- Xhat # then it initialize the fitted version of X with the centered version of X, with random completion
    
    if (ncp == 0) # if the number of components is set to 0, then so is the nb.iter (what's this, and why it goes to 0?)
      nb.iter = 0
    iteration_count=0
    criteria_s=data.frame(myobjective=0,objective=0,criterion=0)
    
    while (nb.iter > 0) {# if nb.iter is greater than 0, then the fitted 
      iteration_count=iteration_count+1
      
      Xhat[missing] <- fittedX[missing]
      # 
      #### un-transform the data ########################################
      if (scale) 
        Xhat = t(t(Xhat) * et) # you de-scale if you scaled before
      
      if(iteration_count>1){
        Xhat <- t(t(Xhat) + mean.p) # you re-add the mean to the data
      }
      
      # b_mini_mean.p <- apply(Xhat[1:18,], 2, mean)
      # a_mini_mean.p <- apply(Xhat[19:50,], 2, mean)
      # 
      mean.p <- apply(Xhat, 2, moy.p, row.w)
      
      
      
      Xhat <- t(t(Xhat) - mean.p)
      et <- apply(Xhat, 2, ec, row.w)
      if (scale) 
        Xhat <- t(t(Xhat)/et)
      ###################################################################################
      
      svd.res <- svd.triplet(Xhat, row.w = row.w, ncp = ncp) ## do the rank-ncp decomposition (take care of weights if needed)
      
      
      
      sigma2 <- nrow(X) * ncol(X)/min(ncol(X), nrow(X) - 1) * sum((svd.res$vs[-c(1:ncp)]^2)/((nrow(X) -  1) * ncol(X) - (nrow(X) - 1) * ncp - ncol(X) * 
                                                                                               ncp + ncp^2))
      
      sigma2 <- min(sigma2 * coeff.ridge, svd.res$vs[ncp + 
                                                       1]^2)### coeff.ridge is the lambda value for the ridge regularization
      
      
      if (method == "em") ## when using em, there is no regularization
        sigma2 <- 0
      
      # shrinked lambda
      lambda.shrinked = (svd.res$vs[1:ncp]^2 - sigma2)/svd.res$vs[1:ncp]
      
      ## This is the reconstruction formula to re-impute the missings in X
      fittedX = tcrossprod(t(t(svd.res$U[, 1:ncp, drop = FALSE] * 
                                 row.w) * lambda.shrinked), svd.res$V[, 1:ncp, 
                                                                      drop = FALSE])
      
      fittedX <- fittedX/row.w ## insert the weights, if needed.
      
      
      diff <- Xhat - fittedX ## difference between the previous centered X and the newly imputed version
      mydiff=diff[missing]
      
      myobjective=sum(mydiff^2)/length(missing)
      
      mycriterion= myobjective
      
      
      diff[missing] <- 0 ## in the missings positions, a 0 is put.
      objective <- sum(diff^2 * row.w) ## squared differences between the previous and the following matrix
      criterion <- abs(1 - objective/old) # compute the improvement
      criteria_s=rbind(criteria_s,c(myobjective,objective,criterion))
      
      myold<-myobjective
      old <- objective # the old objective becomes the current one, and we're ready for a new iteration
      nb.iter <- nb.iter + 1
      
      ### condition to exit the while loop
      if (!is.nan(mycriterion)) {
        if ((myobjective < threshold) && (nb.iter > 5)) 
          nb.iter <- 0
        if ((myobjective < threshold) && (nb.iter > 5)) 
          nb.iter <- 0
      }
      if (nb.iter > maxiter) {
        nb.iter <- 0
        warning(paste("Stopped after ", maxiter, " iterations"))
      }
      
    }# END OF WHILE LOOP
    
    
    ### de-transform the data
    if (scale) 
      Xhat <- t(t(Xhat) * et) #de-scale
    
    
    Xhat <- t(t(Xhat) + mean.p)           ## de-center 
    completeObs <- X                      ## starting matrix
    completeObs[missing] <- Xhat[missing] ## insert in your matrix the imputed values
    
    ### de-transform the fittedX
    if (scale) 
      fittedX <- t(t(fittedX) * et)
    
    fittedX <- t(t(fittedX) + mean.p)
    
    
    ############### SLIGHTLY DIFFERENT FROM WHAT JJ DOES##############
    completeObs[missing] <- fittedX[missing]
    ##################################################################
    
    result <- list()
    result$Xcom <- completeObs
    result$fittedX <- fittedX
    result$criteria=criteria_s[-1,]
    result$imputed=fittedX[missing]
    result$missing=missing
    result$iterations=iteration_count
    return(result)
  }## end Of the impute function
  
  
  
  method <- match.arg(method, c("Regularized", "regularized", 
                                "EM", "em"), several.ok = T)[1]
  obj = Inf
  method <- tolower(method)
  if (ncp > min(nrow(X) - 2, ncol(X) - 1)) 
    stop("ncp is too large")
  if (is.null(row.w)) 
    row.w = rep(1, nrow(X))/nrow(X)
  for (i in 1:nb.init) {
    if (!any(is.na(X))) 
      return(X)
    res.impute = impute(X, ncp = ncp, scale = scale, method = method, 
                        threshold = threshold, seed = if (!is.null(seed)) {
                          (seed * (i - 1))
                        }
                        else {
                          NULL
                        }, init = i, maxiter = maxiter, row.w = row.w, coeff.ridge = coeff.ridge)
    if (mean((res.impute$fittedX[!is.na(X)] - X[!is.na(X)])^2) < 
        obj) {
      res <- res.impute
      obj <- mean((res.impute$fittedX[!is.na(X)] - X[!is.na(X)])^2)
    }
  }
  return(res)
}
