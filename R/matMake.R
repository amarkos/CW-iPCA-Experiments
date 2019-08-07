matMake<-function(n,ps,cate=T,abs=F,prob=.2,Sigma=NULL,Mu=NULL,dims=3,alp=1){
  
  # n:  number of rows
  # ps: vecotor containing the number of attributes
  # probability of presences.
  # if one wants the standardized residual matrix (abs=T)
  #print(alp)
  require("MASS")
  source("R/RandomBinZmake.r")
  source("R/BurtSmart.r")
  Alist=list()
  for(i in 1:length(ps)){
    Alist[[i]]=list()
    if(cate==T){
      Zeta=RandomBinZmake(n,ps[i],probCol=prob)
      Z=Zeta$Z
      A=BurtSmart(Z)
      At=A
      Ar=A/sum(A)
      cm=apply(Ar,1,sum)
      eA=cm%*%t(cm)
      
      if(abs==F){
        A=(Ar-eA)/sqrt(eA)
        
        Alist[[i]]$A=A
      }else{
        
        Alist[[i]]$A=At
      } 
      Alist[[i]]$Z=Z
    }else{
      
      
      
      #Mu=runif(d)             # random mean
      #Sigma=diag(rep(sd2,d))   # variance matrix for the independent variables
      p=ps[i]
      d=dims
#       print(n)
#       print(Mu)
#       print(Sigma)
      X<-mvrnorm(n,Mu,Sigma)
      U<-matrix(rnorm(d*p),nrow=d,ncol=p)
      XU<-X %*% U;# XU has the desired size (n x p) and is d dimensional. 
      # We can add noise to this matrix to get "random" data with underlying d dimensional structure
      noise<-matrix(rnorm(n*p),nrow=n, ncol=p)
        # Controls amount of added noise. 
      # Should probably linked to number of columns (more columns, alp can be higher)
      XUn=XU+alp*noise;
      Xc=scale(XUn,center=F,scale=F)
      
      A=cor(Xc)
      Alist[[i]]$Z=Xc
      Alist[[i]]$A=A
    }
    
    
    
  }
  Alist
}