corr_mat_make<-function(p,vsets_list,hi_rho=.7,low_rho=0){
  
  n_vsets=length(vsets_list)
  if(length(hi_rho)==1){hi_rho=rep(hi_rho,n_vsets)}
  
  corr_mat=matrix(rep(0,p^2),nrow=p)
  diag(corr_mat)=1
  rownames(corr_mat)=paste0("V",1:p)
  colnames(corr_mat)=paste0("V",1:p)
  
  for(j in 1:n_vsets){
  vcomb=t(combn(vsets_list[[j]],2))
  tvcomb=cbind(vcomb[,2],vcomb[,1])
  corr_mat[vcomb]=hi_rho[j]
  corr_mat[tvcomb]=hi_rho[j]
  }
  corr_mat[corr_mat==0]=low_rho
  return(corr_mat)
}