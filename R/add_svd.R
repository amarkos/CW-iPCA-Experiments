add_svd <-  function(eg, B, m, current_rank, ff = 0) {
	
	B = t(B)
	N = nrow(B)
	n = ncol(B)
	ff = 1 -ff
	
	if (missing("current_rank")) {
		#full rank
		current_rank = N
	}
	
	orgn = eg$orgn
	U0 = eg$v[,1:current_rank]
	D0 = eg$d[1:current_rank]  
	V0 = eg$u[,1:current_rank]
	orgnb = rowMeans(B) 
	#center data
	B = B - as.matrix(orgnb) %*% as.matrix(t(rep(1,n)))
	
	B <- cbind(B,sqrt(n*m/(n+m))*as.matrix((orgnb-orgn)))
	#mean update
	orgnc <- (ff*m*orgn + n*orgnb)/(n+ff*m)
	
	B_proj = t(U0)%*%B 
	B_res = B - U0%*%B_proj
	qrstr = qr(B_res)
	q = qr.Q(qrstr,complete=T)
	Q = cbind(U0,q)
	R = rbind(cbind(ff*diag(D0),B_proj),cbind(matrix(0,nrow(B),length(D0)),t(q)%*%B_res))
	
	eg12 = fast.svd(R, 0) 
	
	D = eg12$d[1:current_rank]
	U = Q %*% eg12$u[, 1:current_rank]
	eg12$v = eg12$v[,1:current_rank]
	#these left eigenvectors are not exact
	V = rbind(V0 %*% eg12$v[1:current_rank,],eg12$v[(current_rank+1):(dim(eg12$v)[1]-1),])   #Exploit structure to compute this fast  Vp = [ Vp ; tVp( current_rank+1:size(tVp,1), : ) ];
	
	m = n + ff*m
	
	#force_orthogonality
	# qrUp = qr(U)
	# UQ = qr.Q(qrUp)
	# UR = qr.R(qrUp)
	# qrVp = qr(V)
	# VQ = qr.Q(qrVp)
	# VR = qr.R(qrVp)
	# # 
	# eg = fast.svd(UR%*%diag(D)%*%t(VR),0)
	# tUp = eg$u
	# tSp = eg$d
	# tVp = eg$v
	# U = UQ %*% tUp
	# V = VQ %*% tVp
	# D = tSp
	
	out = list()
	out$u <- V
	out$d <- D
	out$v <- U
	out$m <- m
	out$orgn <- orgnc
	out
}
