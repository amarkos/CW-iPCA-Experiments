frob_rel_err=function(trueV,V){
  
  err=norm(trueV%*%t(trueV) - V%*%t(V),type="F")
  err=err/norm(trueV%*%t(trueV),type="F")
  return(err)
}