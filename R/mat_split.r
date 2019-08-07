mat_split<- function(mat,nchunk){
  mat.story=list()
  n=nrow(mat)
  
  #case when input is number of chunks 
  if (length(nchunk) == 1) {
    # seq.chu[nchunk]=n 
    # print(nc)
    
    
    seq.chu=floor(seq(from=1,to=n,length.out =nchunk+1))
     
      
    
    # seq.chu[nchunk]=n 
    #seq.chu=seq.chu[1:nchunk]
    # print(length(seq.chu))
    
    
  #  print(seq.chu)
    # print(length(seq.chu))
  }
  #case when input is specific block sizes
  if (length(nchunk) > 1){
  #  nc = length(nchunk)
    seq.chu=0
    seq.chu=c(seq.chu,cumsum(nchunk))+1
    #fix last value
    seq.chu[length(seq.chu)] = n
  }
  
  kk=length(seq.chu)
  
  if (kk > 2)
  {
    mat.story[[1]]=mat[1:seq.chu[2],]
  
    for(k in 2:(kk-1)){
      
      a=(seq.chu[k]+1)
      b=(seq.chu[k+1])
      # print("a")
      # print(a)
      # print("b")
      # print(b)
      mat.story[[k]]=mat[a : b, ]
      }
    
  } else { #case of a single chunk
    mat.story[[1]] = mat  
  }
  out=list()
  out$splitMat=mat.story
  
  out
}


#     mat.story[[2]]=mat2[92 : 181, ]
#print(mat.story)
#
