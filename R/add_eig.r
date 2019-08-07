add_eig <- function(m1, m2, current_rank) {
  
  #output
  out = list()
  
  m1$u=m1$u[,1:current_rank]
  m1$v=m1$v[,1:current_rank]
  m1$d=m1$d[1:current_rank]
  m2$u=m2$u[,1:current_rank]
  
  m2$v=m2$v[,1:current_rank]
  m2$d=m2$d[1:current_rank]
  
  bign = m1$m + m2$m

  if(is.null(m1$true_m) & is.null(m2$true_m)){
    bigOrgn = ((as.vector(m1$orgn) * m1$m) + (as.vector(m2$orgn) * m2$m))/bign  
  }
  if(!is.null(m1$true_m)){
    big_true_n = m1$true_m + m2$m
    bigOrgn = ((as.vector(m1$orgn) * as.vector(m1$true_m)) + as.vector(m2$orgn * m2$m))/big_true_n
  }
  if(!is.null(m2$true_m)){
    big_true_n = m1$m + m2$true_m
    bigOrgn = (as.vector(m1$orgn * m1$m) + (as.vector(m2$orgn) * as.vector(m2$true_m)))/big_true_n
  }
  dorgn = as.vector(m1$orgn) - as.vector(m2$orgn)
  G = t(m1$v) %*% m2$v
  H = m2$v - m1$v %*% G
  g = t(m1$v) %*% dorgn
  h = dorgn - m1$v %*% g
  eps = 2.2204e-16
  HH = apply(H * H, 2, sum)
  hh = apply(h * h, 2, sum)
  #if (any(HH > eps) | any(hh > eps)) {
  #   source("orth.r")
  #  Hnz = which(HH > eps)
  #  hnz = which(hh > eps)
  #  nu = orth(cbind(HH[, Hnz], hh[, hnz]))
  #  H = t(m1$v) %*% nu
  #  nu = nu[, -Hnz]
  #} else {
  nu = matrix(0, length(m1$orgn), 0)
  # }
  resn1 = ncol(m1$v) - nrow(m1$v)
  #  if (resn1 > 0) {
  #    rpern1 = m1.resudue/resn1
  #  } else {
  rpern1 = 0
  #  }
  resn2 = ncol(m2$v) - nrow(m2$v)
  #  if (resn2 > 0) {
  #    rpern2 = m2.resudue/resn2
  #  } else {
  rpern2 = 0
  #  }
  nnu = nrow(nu)
  mnu = ncol(nu)
  Gamma = t(nu) %*% m2$v
  
  D = t(t(G) * as.vector(m2$d))
  E = t(t(Gamma) * as.vector(m2$d))
  gamma = t(nu) %*% dorgn
  
  A11 = cbind(m1$d * t(m1$u), G %*% (m2$d * t(m2$u)) )
  A12b = (t(nu) %*% m2$v) %*% (m2$d * t(m2$u))
  t = nrow(A12b)
  A12a = matrix(0, t, (m1$m))
  A12 = cbind(A12a, A12b)
  A1 = rbind(A11, A12)
  onerX = matrix(1, 1, m1$m)
  onerY = matrix(1, 1, m2$m)
  #dorgnX = m1$orgn - bigOrgn
  dorgnX = as.double(m1$orgn) - as.double(bigOrgn)
  dorgnY = as.double(m2$orgn) - as.double(bigOrgn)
  A21a = (t(m1$v) %*% as.double(dorgnX) %*% onerX)
  A21b = (t(m1$v) %*% as.double(dorgnY) %*% onerY)
  
  A21 = cbind(A21a, A21b)
  A22a = t(nu) %*% as.double(dorgnX) %*% onerX
  A22b = t(nu) %*% as.double(dorgnY) %*% onerY
  A22 = cbind(A22a, A22b)
  A2 = rbind(A21, A22)
  A = A1 + A2
  svd.res = fast.svd(A,0)
  
  v = svd.res$u
  
  v = cbind(m1$v, nu) %*% v
  d = svd.res$d
  u = svd.res$v
  out$u = u
  
  out$m = bign
  out$orgn = bigOrgn
  out$v = v
  out$d = d
  # out$A = A
  out
}  






# add_eig <- function(m1, m2,current_rank=5) {
#   # print('in')
#   #output
#   out = list()
#   
#   m1$u=m1$u[,1:current_rank]
#   m1$v=m1$v[,1:current_rank]
#   m1$d=m1$d[1:current_rank]
#   
#   m2$u=m2$u[,1:current_rank]
#   m2$v=m2$v[,1:current_rank]
#   m2$d=m2$d[1:current_rank]
#   
#   
#   bign = m1$m + m2$m
#   if(is.null(m1$true_m)&is.null(m2$true_m)){
#     bigOrgn = ((m1$orgn * m1$m) + (m2$orgn * m2$m))/bign  
#   }
#   if(!is.null(m1$true_m)){
#     big_true_n = m1$true_m + m2$m
#   bigOrgn = ((as.vector(m1$orgn) * as.vector(m1$true_m)) + as.vector(m2$orgn * m2$m))/big_true_n
#   }
#   if(!is.null(m2$true_m)){
#     big_true_n = m1$m + m2$true_m
#     bigOrgn = (as.vector(m1$orgn * m1$m) + (as.vector(m2$orgn) * as.vector(m2$true_m)))/big_true_n
#   }
#   
#   
#   
#   dorgn = m1$orgn - m2$orgn
#   G = t(m1$v) %*% m2$v
#   H = m2$v - m1$v %*% G
#   g = t(m1$v) %*% dorgn
#   h = dorgn - m1$v %*% g
#   eps = 2.2204e-16
#   HH = apply(H * H, 2, sum)
#   hh = apply(h * h, 2, sum)
#   #if (any(HH > eps) | any(hh > eps)) {
#     #   source("orth.r")
#   #  Hnz = which(HH > eps)
#   #  hnz = which(hh > eps)
#   #  nu = orth(cbind(HH[, Hnz], hh[, hnz]))
#   #  H = t(m1$v) %*% nu
#   #  nu = nu[, -Hnz]
#   #} else {
#     nu = matrix(0, length(m1$orgn), 0)
#   # }
#   resn1 = ncol(m1$v) - nrow(m1$v)
#   #  if (resn1 > 0) {
#   #    rpern1 = m1.resudue/resn1
#   #  } else {
#   rpern1 = 0
#   #  }
#   resn2 = ncol(m2$v) - nrow(m2$v)
#   #  if (resn2 > 0) {
#   #    rpern2 = m2.resudue/resn2
#   #  } else {
#   rpern2 = 0
#   #  }
#   nnu = nrow(nu)
#   mnu = ncol(nu)
#   Gamma = t(nu) %*% m2$v
#   
#   D = t(t(G) * as.vector(m2$d))
#   E = t(t(Gamma) * as.vector(m2$d))
#   gamma = t(nu) %*% dorgn
#  
#   A11 = cbind(m1$d * t(m1$u), G %*% (m2$d * t(m2$u)) )
#   A12b = (t(nu) %*% m2$v) %*% (m2$d * t(m2$u))
#   t = nrow(A12b)
#   A12a = matrix(0, t, (m1$m))
#   A12 = cbind(A12a, A12b)
#   A1 = rbind(A11, A12)
#   onerX = matrix(1, 1, m1$m)
#   onerY = matrix(1, 1, m2$m)
#   #dorgnX = m1$orgn - bigOrgn
#   dorgnX = as.double(m1$orgn) - as.double(bigOrgn)
#   dorgnY = as.double(m2$orgn) - as.double(bigOrgn)
#   A21a = (t(m1$v) %*% as.double(dorgnX) %*% onerX)
#   A21b = (t(m1$v) %*% as.double(dorgnY) %*% onerY)
#  
#   A21 = cbind(A21a, A21b)
#   A22a = t(nu) %*% as.double(dorgnX) %*% onerX
#   A22b = t(nu) %*% as.double(dorgnY) %*% onerY
#   A22 = cbind(A22a, A22b)
#   A2 = rbind(A21, A22)
#   A = A1 + A2
#   #library(svd)
#   # print(dim(A))
#   svd.res = fast.svd(A,0)
#   
#   #   svd.res = irlba(A,nrow(A),nrow(A))
#   #svd.res = svd(A)
#   v = svd.res$u
# 
#   v = cbind(m1$v, nu) %*% v
#   #v=sign.matcher(m1$v,v)
#   d = svd.res$d
#   u = svd.res$v
#   # print(dim(u))
#   out$u = u
#   #  }
#   
#   out$m = bign
#   out$orgn = bigOrgn
#   out$v = v
#   out$d = d
#   # out$A = A
#   out
# }  
# 
