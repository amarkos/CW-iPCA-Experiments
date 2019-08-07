MNCAR_missings_experiment=function(nn=500, nblocks=5, ncp=3,p=9,myseed=1234,MCAR=FALSE,
                                   full_rand_vals=NULL,perc=.75,Vsets=NULL,pointer_blocks=NULL,
                                   hi_cor=c(.75,.6,.8),low_cor=.15,
                                   prop.m = .2 ){
  
  library("naniar")
  library("visdat")
  library("purrr")
  library("ggcorrplot")
  source("R/weaken_corr_by_NA.R")
  source("R/my_generate_missing.R")
  source("R/corr_mat_make.R")
  source("R/Impute_Experiment_multi_block.R")  
  # fullCor=NULL
  # prop_mis=.25
  # % missingness (0.1, 0.2, 0.5)
  set.seed(myseed)
  n=nn*nblocks
  # Vsets within correlation
  
  if(is.null(Vsets)){
    Vsets=list(paste0("V",1:4),
               paste0("V",5:7),
               paste0("V",8:9)
    )
  }
  
  if(is.null(pointer_blocks)){
    pointer_blocks=tibble(blocks_to= c(3,3,4,5,5,6,6,7,8,9,9,10),
                          vsets_to = c(1,3,2,1,2,2,2,3,2,2,3,3))
    pointer_blocks=pointer_blocks[1:nblocks,]
  }
  
  
  
  
  
  
  custom_corr_mat = corr_mat_make(p = p, vsets_list = Vsets, hi_rho = hi_cor, low_rho = low_cor)
  U = t(chol(custom_corr_mat))
  random.normal = matrix(rnorm(p*n,0,1), nrow=p, ncol=n);
  X = U %*% random.normal
  X = t(X)
  sdX=apply(X,2,sd)
  X=as_tibble(X/sdX)
  
  blocks = factor(rep(paste0("block",1:nblocks),each=nn),levels=paste0("block",1:nblocks))
  X_by_block=split(x=X,f = blocks)
  miss_X_by_block=split(x=X,f = blocks)
  
  
  if(MCAR){
    mcar_blocks = 1:nblocks
  }else{
    mcar_blocks= setdiff(1:nblocks,as.matrix(pointer_blocks[,1]))
  }
  for(i in 1:length(mcar_blocks)){
    for (j in 1:ncol(miss_X_by_block[[mcar_blocks[i]]])) {
      miss_X_by_block[[mcar_blocks[i]]][sample(1:nn,size = round(nn*prop.m),replace = FALSE),j] = NA
    }
  }
  
  if(!MCAR){
    for(i in 1:nrow(pointer_blocks)){
      bb=as.numeric(pointer_blocks[i,1])
      vv=as.numeric(pointer_blocks[i,2])
      miss_X_by_block[[bb]][,Vsets[[vv]]]=weaken_corr_by_NA(X_by_block[[bb]][,Vsets[[vv]]],perc)
    }
  }
  Xmis=bind_rows(miss_X_by_block)
  xmis_snapshot=vis_dat(Xmis)
  xmis_snapshot_bw=vis_miss(Xmis,sort_miss=TRUE)
  xmis_missingness=gg_miss_var(Xmis)
  outmis=list(X=as.matrix(X),Xmis=as.matrix(Xmis),ncp=ncp)
  custom_corplot=ggcorrplot(custom_corr_mat)
  # ggcorrplot(cor(outmis$Xmis,use="complete.obs"))
  # ggcorrplot(cor(outmis$X))
  # ggcorrplot(cor(miss_X_by_block[[5]],use="complete.obs"))
  
  
  
  
  exp_out = Impute_Experiment_multi_block(matList = outmis,n = nn*nblocks, p = p, ncp = ncp, 
                                          nblocks = nblocks, prop_mis = prop_mis, myseed=myseed) 
  # mechanism = mechanism, scenario = scenario)  
  # save(update_out,file="R/pdf_update_result.RData")
  exp_out$error
  exp_out$indRVs
  exp_out$varRVs
  exp_out$cor_structure=custom_corplot
  exp_out$data_snapshot=xmis_snapshot
  exp_out$data_snapshot_bw=xmis_snapshot_bw
  exp_out$data_missingness=xmis_missingness
  
  
  return(exp_out)
}
