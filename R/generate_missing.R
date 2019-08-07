#Generating normally distributed data with a specific correlation structure
# Based on Dray & Josse (2015) Principal component analysis with missing values:
# a comparative survey of methods.
generate_missing <- function(n, p, fix_means=0, prop_mis, mechanism, scenario,myseed) {
  set.seed(myseed)
  # p: number of variables {9, 18, 45}
  # n: number of individuals {20, 50, 100}
  # r: correlation within blocks
  # mechanism: missing data mechanism {MCAR, MNAR}
  # scenario: {M1, M2, M3}
  out = list()
  
  if (scenario == "M1") {
    #Scenario 1 (M1)
    r  = .8 #correlation within blocks
    c1 = 4*p/9 
    c2 = 3*p/9 
    c3 = 2*p/9 
    R1 = diag(c1)
    R1[R1 == 0] = r
    R2 = diag(c2)
    R2[R2 == 0] = r
    R3 = diag(c3)
    R3[R3 == 0] = r
    #R1 = matrix(cbind(1,r,r,r,  r,1,r,r,  r,r,1,r, r,r,r,1),nrow=4)
    #R2 = matrix(cbind(1,r,r,  r,1,r,  r,r,1),nrow=3)
    #R3 = matrix(cbind(1,r,  r,1),nrow=2)
    M = rbind(cbind(R1,matrix(0,c1,c2),matrix(0,c1,c3)),cbind(matrix(0,c2,c1),R2,matrix(0,c2,c3)), cbind(matrix(0,c3,c1),matrix(0,c3,c2),R3))
    ncp = 3
  } else if (scenario == "M2") {
    #Scenario 2 (M2)
    r  = .3
    c1 = 4*p/9 
    c2 = 3*p/9 
    c3 = 2*p/9 
    R1 = diag(c1)
    R1[R1 == 0] = r
    R2 = diag(c2)
    R2[R2 == 0] = r
    R3 = diag(c3)
    R3[R3 == 0] = r
    #R1 = matrix(cbind(1,r,r,r,  r,1,r,r,  r,r,1,r, r,r,r,1),nrow=4)
    #R2 = matrix(cbind(1,r,r,  r,1,r,  r,r,1),nrow=3)
    #R3 = matrix(cbind(1,r,  r,1),nrow=2)
    M = rbind(cbind(R1,matrix(0,c1,c2),matrix(0,c1,c3)),cbind(matrix(0,c2,c1),R2,matrix(0,c2,c3)), cbind(matrix(0,c3,c1),matrix(0,c3,c2),R3))
    ncp = 3
  } else if (scenario == "M3") {
    #Scenario 3 (M3)
    r  = .8
    M = diag(p)
    M[M == 0] = r
    #    M = matrix(cbind(1,r,r,r,r,r,r,r,r,  r,1,r,r,r,r,r,r,r,  r,r,1,r,r,r,r,r,r, r,r,r,1,r,r,r,r,r, r,r,r,r,1,r,r,r,r, r,r,r,r,r,1,r,r,r, r,r,r,r,r,r,1,r,r, r,r,r,r,r,r,r,1,r, r,r,r,r,r,r,r,r,1),nrow=9)
    ncp = 1
  }
  
  ## https://www.r-bloggers.com/simulating-random-multivariate-correlated-data-continuous-variables/
  #This approach takes an original X variable (or matrix) and uses the Cholesky transformation to create a new, correlated, Y variable. 
  # Exploit the Cholesky transformation
  U = t(chol(M))
  # random.normal = matrix(rnorm(p*n,0,1), nrow=p, ncol=n);
  random.normal = matrix(rnorm(p*n,0,1), nrow=p, ncol=n);
  X = U %*% random.normal
  X = t(X)
  ### Careful..
  
  X = scale(X,scale=TRUE,center=FALSE)
  
  ############
  if(is.null(fix_means)==F){
  fix_means=round(runif(p,min=2,max=20))  
  onen=t(t(rep(1,nrow(X))))
  meanX=fix_means
  MX=onen%*%meanX
  X=X+MX
  }
  ###########
  
  
  if (mechanism == "MCAR") {
    #Create the incomplete dataset
    # Case 1: MCAR
    #The idea: 
    #For 10% of your data missing completely at random (MCAR), 
    #draw a number from a uniform distribution each value of a specific variable 
    #and if it is <.10, replace the value with NA.
    Xmis = X
    prop.m = prop_mis  # % missingness (0.1, 0.2, 0.5)
    for (j in 1:ncol(X)) {
      # mcar   = runif(n, min=0, max=1)
      # Xmis[,j] = ifelse(mcar < prop.m, NA, X[,j])
      Xmis[sample(1:n,size = round(n*prop.m),replace = FALSE),j] = NA
    }
    #nb <- estim_ncpPCA(Xmis)$ncp  
  } else if (mechanism == "MNAR") {
    #Case 2: MNAR
    #MNAR 20% of the highest values of the 1st variable are missing
    #sort in ascending order
    Xmis = X
    sind = sort(X[,1],index.return=TRUE)$ix
    vals = round(prop_mis*n)
    sind = sind[-c(1:(length(sind)-vals-1))]
    Xmis[sind,1] = NA
    #nb <- estim_ncpPCA(Xmis)$ncp
  }
  out$Xmis = Xmis
  out$X = X
  out$ncp = ncp
  out
}
