
source("R/MNCAR_missings_experiment.R")

if (!require("naniar")) install.packages("naniar")
if (!require("ggcorrplot")) install.packages("ggcorrplot")
if (!require("missMDA")) install.packages("missMDA")
if (!require("factoextra")) install.packages("factoextra")


Vsets=list(paste0("V",1:4),
           paste0("V",5:7),
           paste0("V",8:9)
)

# nb=10
nn=500
ncp=3
p=9
hi_cor=c(.75,.6,.8)
low_cor=.15
number_of_replicates=10
nb_vec=c(5,10)#,15,20,25)

for(reps in 1:number_of_replicates){
  print(paste("within the replicate",reps))
  #### generate pointer_blocks
  myseed=123+reps
  set.seed(myseed)
  sam_size=36
  mysam=sample(1:3,size=sam_size,replace=T,prob = c(.34,.33,.33))
  # mysam[3]=3
  diff_mysam=c(0,diff(mysam))>0
  
  vsets_to=mysam
  blocks_to=rep(0,sam_size)
  blocks_to[!diff_mysam]=3:(sum(!diff_mysam)+2)
  wzero=which(blocks_to==0)
  w_not_zero=which(blocks_to!=0)
  for(i in 1:length(wzero)){
    # print(tail(which(wzero[i]>w_not_zero),1))
    # print(blocks_to[tail(which(wzero[i]>w_not_zero),1)])
    blocks_to[wzero[i]]=blocks_to[w_not_zero[tail(which(wzero[i]>w_not_zero),1)]]  
  }
  
  
  pointer_blocks=tibble(blocks_to,vsets_to)
  
  
  waiting_time<-system.time({
    for(jj in 1:length(nb_vec)){
      print("entered nblocks")
      print(nb_vec[jj])
      nb=nb_vec[jj]
      sub_pointer_blocks=pointer_blocks[pointer_blocks$blocks_to<=nb,]
      
      
      
      experiment_out=MNCAR_missings_experiment(nn=nn, nblocks=nb, ncp=3,p=9,myseed=myseed,MCAR=TRUE,
                                               perc=.75,Vsets=Vsets,
                                               hi_cor=hi_cor,low_cor=low_cor,
                                               pointer_blocks=sub_pointer_blocks,
                                               prop.m = .2 )
      #full_rand_vals=NULL,
      exp_name=paste("MCAR_replicate",reps,"exp","n",nn,"p",p,"nblocks",nb,sep="_")
      assign(exp_name[1],experiment_out)
      save(file=paste0("R/",exp_name,".RData"),experiment_out)
    }
  })
  
}

