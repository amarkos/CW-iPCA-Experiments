Impute_Experiment_multi_block<-function(matList=NULL,n = 100,p = 9,ncp=5, prop_mis = .05, mechanism = "MCAR", scenario = "M1",nblocks,myseed=1234){
  
  source("R/naive_impute_PCA.R")#
  source("R/generate_missing.R")#
  source("R/my_generate_missing.R")#
  # source("R/my_imputePCA_myob.R")
  source("R/do_es.r")#
  source("R/add_es.r")#
  source("R/old_add_eig.R")#
  source("R/add_eig.R")#
  source("R/block_impute_PCA.R")#
  source("R/multi_naive_impute_PCA.R")#
  source("R/multi_block_impute_PCA.R")#
  source("R/exact_multi_block_impute_PCA.R")#
  source("R/multi_JJ_impute_PCA.R")#
  source("R/mat_split.r")#
  
  
  source("R/frob_rel_err.R")
  require(missMDA)
  require(FactoMineR)
  require(factoextra)
  require(corpcor)
  
  #set.seed(myseed)
  
  labels=matrix(1:(n*p),nrow=n,ncol=p)
  if(is.null(matList)==T){
    outmis = generate_missing(n = n, p = p, prop_mis = prop_mis, mechanism = mechanism, scenario = scenario,myseed=myseed)
    outmis$allXmis=outmis$Xmis
    outmis$allX=outmis$X
  }else{
    outmis=matList
    outmis$allXmis=outmis$Xmis
    outmis$allX=outmis$X
  } 
  # print(dim(outmis$Xmis))
  
  rowna=sort(unique(which(is.na(outmis$Xmis),arr.ind = T)[,1]))
  All_missings=which(is.na(outmis$Xmis))
  
  All_missings_full=which(is.na(outmis$allXmis))
  All_mis_lab=labels[All_missings]
  bigXmis=outmis$allXmis
  all_bigX=outmis$allX
  trueX=outmis$X
  # print(apply(outmis$X,2,mean))
  ############
  # print(nblocks)
  outX=mat_split(outmis$Xmis,nchunk = nblocks)
  Xlist=outX$splitMat
  outLabList=mat_split(labels,nchunk = nblocks)
  LabList=outLabList$splitMat
  
  # print(LabList)
  # print(length(outX$splitMat))
  #############
  # print(lapply(Xlist,dim))
  out_full_jj = imputePCA(bigXmis,scale=F,seed=1234,ncp=ncp,threshold=1e-12,maxiter = 10000) #  imputePCA
  # out_full = naive_impute_PCA(bigXmis,scale=F,seed=1234,ncp=ncp,threshold=1e-12,maxiter = 10000)#  imputePCAobj
  out_naive = multi_naive_impute_PCA(Xlist,LabList,scale=F,seed=1234,ncp=ncp,threshold=1e-12,maxiter = 10000) #naive
  # out_multi_jj = multi_JJ_impute_PCA(Xlist,LabList,scale=F,seed=1234,ncp=ncp,threshold=1e-12,maxiter = 10000) #inc_imp_JJ_PCA
  out_block = multi_block_impute_PCA(Xlist,LabList,scale=F,seed=1234,ncp=ncp+1,threshold=1e-12,maxiter = 10000) # non naive #1
  # out_ex_block = exact_multi_block_impute_PCA(Xlist,scale=F,seed=1234,ncp=ncp,threshold=1e-12,maxiter = 10000) #non naive #2
  
  # out_full$imputed_vals=out_full$Xcom[All_missings_full]
  out_full_jj$imputed_vals=out_full_jj$completeObs[All_missings_full]
  out=list()
  # out$multi_JJ=out_multi_jj
  out$no_naive=out_block
  out$naive=out_naive
  out$out_full_jj=out_full_jj
  # out$out_full_naive=out_full
  out$true_vals=trueX[All_missings]
  out$all_true_vals=all_bigX[All_missings_full]
  
  
  #############################################################################
  #############################################################################
  ############# computing imputating errors #################  
  ###########################################################
  imputed_mat=data.frame(imputePCA = out_full_jj$imputed_vals,
                         # imputePCA_myob = out_full$imputed_vals,
                         # inc_imp_ex_non_naive = out_ex_block$imputed_vals,
                         # inc_imp_JJ_PCA = out_multi_jj$imputed_vals,
                         inc_imp_non_naive = out_block$imputed_vals,
                         inc_imp_naive = out_naive$imputed_vals,
                         true_vals = out$true_vals
  )
  # JJvsNaive=cbind(out_multi_jj$imputed_vals,out_block$imputed_vals)
  errors=data.frame(
    imputePCA            = mean(abs(round(out_full_jj$imputed_vals,8)-round(out$all_true_vals,8))),
    # imputePCA_myob       = mean(abs(round(out_full$imputed_vals,8)-round(out$all_true_vals,8))),
    # inc_imp_JJ_PCA       = mean(abs(round(out_multi_jj$imputed_vals,8)-round(out$true_vals,8))),
    inc_imp_non_naive    = mean(abs(round(out_block$imputed_vals,8)-round(out$true_vals,8))),
    inc_imp_naive        = mean(abs(round(out_naive$imputed_vals,8)-round(out$true_vals,8)))
    # inc_imp_ex_non_naive = mean(abs(round(out_ex_block$imputed_vals,8)-round(out$true_vals,8)))
    
    
    #    imputePCA=mean(abs(out_full_jj$imputed_vals)-abs(out$true_vals)),
    #    imputePCA_myob=mean(abs(out_full$imputed_vals)-abs(out$true_vals)),
    #    inc_imp_JJ_PCA=mean(abs(out_multi_jj$imputed_vals)-abs(out$true_vals)),
    #    inc_imp_naive=mean(abs(out_naive$imputed_vals)-abs(out$true_vals)),
    #    inc_imp_non_naive=mean(abs(out_block$imputed_vals)-abs(out$true_vals)),
    #    inc_imp_ex_non_naive=mean(abs(out_ex_block$imputed_vals)-abs(out$true_vals))
    
  )
  
  # print("ENTERED RV COMPUTATION")
  #### RV calculation - INDIVIDUALS ####
  truePCA = PCA(trueX,scale.unit = FALSE,graph=FALSE,ncp = ncp)$ind$coord
  out_fullPCA = PCA(out_full_jj$completeObs,scale.unit = FALSE,graph=FALSE,ncp = ncp)$ind$coord
  # out_fullPCAobj = PCA(out_full$Xcom,scale.unit = FALSE,graph=FALSE,ncp = ncp)$ind$coord
  # out_JJPCA = PCA(out_multi_jj$Ximputed,scale.unit = FALSE,graph=FALSE,ncp = ncp)$ind$coord
  out_nonnaivePCA = PCA(out_block$Xfull,scale.unit = FALSE,graph=FALSE,ncp = ncp)$ind$coord
  out_naivePCA = PCA(out_naive$Xcom,scale.unit = FALSE,graph=FALSE,ncp = ncp)$ind$coord
  # out_exnonnaivePCA = PCA(out_ex_block$Xcom,scale.unit = FALSE,graph=FALSE,ncp = ncp)$ind$coord
  
  indRVs=data.frame(
    imputePCA            = coeffRV(out_fullPCA,truePCA)$rv, 
    # imputePCA_myob       = coeffRV(out_fullPCAobj,truePCA)$rv, 
    # inc_imp_JJ_PCA       = coeffRV(out_JJPCA,truePCA)$rv, 
    inc_imp_non_naive    = coeffRV(out_nonnaivePCA,truePCA)$rv,
    inc_imp_naive        = coeffRV(out_naivePCA,truePCA)$rv
    # inc_imp_ex_non_naive = coeffRV(out_fullPCA,truePCA)$rv #same with imputePCA  
    
  )
  
  # #### RV calculation - ATTRIBUTES ####
  truePCA = PCA(trueX,scale.unit = FALSE,graph=FALSE,ncp = ncp)$var$coord
  out_fullPCA = PCA(out_full_jj$completeObs,scale.unit = FALSE,graph=FALSE,ncp = ncp)$var$coord
  # out_fullPCAobj = PCA(out_full$Xcom,scale.unit = FALSE,graph=FALSE,ncp = ncp)$var$coord
  # out_JJPCA = PCA(out_multi_jj$Ximputed,scale.unit = FALSE,graph=FALSE,ncp = ncp)$var$coord
  out_nonnaivePCA = PCA(out_block$Xfull,scale.unit = FALSE,graph=FALSE,ncp = ncp)$var$coord
  out_naivePCA = PCA(out_naive$Xcom,scale.unit = FALSE,graph=FALSE,ncp = ncp)$var$coord
  # out_exnonnaivePCA = PCA(out_ex_block$Xcom,scale.unit = FALSE,graph=FALSE,ncp = ncp)$ind$coord
  # 
  varRVs=data.frame(
    imputePCA            = coeffRV(out_fullPCA,truePCA)$rv, 
    # imputePCA_myob       = coeffRV(out_fullPCAobj,truePCA)$rv, 
    # inc_imp_JJ_PCA       = coeffRV(out_JJPCA,truePCA)$rv, 
    inc_imp_non_naive    = coeffRV(out_nonnaivePCA,truePCA)$rv,
    inc_imp_naive        = coeffRV(out_naivePCA,truePCA)$rv
    # inc_imp_ex_non_naive = coeffRV(out_fullPCA,truePCA)$rv #same with imputePCA  
    
  )
  
  ###############################################################################
  ###############################################################################
  #  print(out_multi_jj$Ximputed)#[1:10,1:2])
  #  print(out_block$Xcom)#[1:10,1:2])
  # print(round(abs(out_multi_jj$imputed_vals-out_block$imputed_vals),2))
  
  out$errors=errors
  out$indRVs=indRVs
  out$varRVs=varRVs
  out$imputed_mat=imputed_mat
  out$All_mis_lab=All_mis_lab
  # out$JJvsNaive=JJvsNaive
  out$Xlist=Xlist
  out$LabList=LabList
  ###########################################################
  ###########################################################
  
  return(out)
  
}
