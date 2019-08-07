rm(list = ls())
source("R/TEP_experimenter.r")

nb=20
tot_size = 500 * nb

perf = TEP_experimenter(n_sam = tot_size, nblocks=nb, myseed=1234, ncp = 5, rep=1)
for(i in 2:10){
rep_perf = TEP_experimenter(n_sam = tot_size, nblocks=nb, myseed=1234, ncp = 5, rep=i)
perf = full_join(perf, rep_perf) 
}

perf_name = paste0("perf_n_",tot_size,"_chunks_",nb)
perf_fname = paste0("R/",perf_name,".RData")
assign(perf_name,perf)
save(list=perf_name, file=perf_fname)

