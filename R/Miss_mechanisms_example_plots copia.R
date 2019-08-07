#load the first replicate of each of the scenarios, 25 blocks
rm(list=ls())
library("dplyr")
library("tidyverse")
library("gridExtra")


load("R/MNCAR_bad_start_experiments_results/MNCAR_bad_start_replicate_1_exp_n_500_p_9_nblocks_25.RData")
assign("MNCAR_bs_example",experiment_out)                

load("R/MNCAR_experiments_results/MNCARreplicate_1_exp_n_500_p_9_nblocks_25.RData")
assign("MNCAR_example",experiment_out)                

load("R/MCAR_experiments_results/MCAR_replicate_1_exp_n_500_p_9_nblocks_25.RData")
assign("MCAR_example",experiment_out)                

#correlation structures (they're the same)
MNCAR_bs_example$cor_structure
MNCAR_example$cor_structure
MCAR_example$cor_structure
ggsave("./CARME2019/Corr_structure.pdf",MCAR_example$cor_structure,height=8,width=8)

#missingness description
require(tidyverse)
source('~/Dropbox/idm_Stanford/R/snap_lab_maker.R')
source('~/Dropbox/idm_Stanford/R/block_highlight.R')

MCAR_ggsnap = MCAR_example$data_snapshot_bw
MCAR_labs = snap_lab_maker(MCAR_ggsnap,bw=T)



MCAR_ggsnap = MCAR_ggsnap + scale_x_discrete(labels = MCAR_labs)+
  theme(axis.text.x = element_text(angle=-70))
tgt_chunk = 1:25
tgt_block = rep(1,25)
outgg = block_highlight(MCAR_ggsnap,tgt_chunk , tgt_block, block_fill="red",border_box="grey40")
tgt_block = rep(2,25)
outgg = block_highlight(outgg,tgt_chunk , tgt_block, block_fill="green",border_box="grey40")
tgt_block = rep(3,25)
outgg = block_highlight(outgg,tgt_chunk , tgt_block, block_fill="blue",border_box="grey40")

ggsave("./CARME2019/CS_incremental_PCA_with_missings/figures/MCAR_exampleA.pdf",outgg,height=8,width=6)
####################################################################################################

MCAR_ggsnap = MCAR_example$data_snapshot
snap_build= ggplot_build(MCAR_ggsnap)
data_snap = snap_build_bw$data[[1]]

# MNCAR
####################################################################################################
####################################################################################################
MNCAR_ggsnap = MNCAR_example$data_snapshot_bw
MNCAR_labs = snap_lab_maker(MNCAR_ggsnap)

MNCAR_ggsnap = MNCAR_ggsnap + scale_x_discrete(labels = MNCAR_labs)+
  theme(axis.text.x = element_text(angle=-70))+theme(legend.position="none")
tgt_chunk = c(1,1,1,2,2,2)
tgt_block = c(1,2,3,1,2,3)
outgg = block_highlight(MNCAR_ggsnap,tgt_chunk , tgt_block, block_fill="orange",border_box="grey40")

ggsave("./CARME2019/CS_incremental_PCA_with_missings/figures/MNCAR_exampleA.pdf",outgg,height=8,width=6)


####################################################################################################
MNCAR_bs_ggsnap = MNCAR_bs_example$data_snapshot_bw
MNCAR_bs_labs = snap_lab_maker(MNCAR_bs_ggsnap,bw=T)

MNCAR_bs_ggsnap = MNCAR_bs_ggsnap + scale_x_discrete(labels = MNCAR_bs_labs)+
  theme(axis.text.x = element_text(angle=-70)) + theme(legend.position="none")

tgt_chunk = c(3,3,3,9,9,9,12,12,12,18,18,18,23,23,23)
tgt_block = c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)
outgg = block_highlight(MNCAR_bs_ggsnap,tgt_chunk , tgt_block, block_fill="orange",border_box="grey40")

ggsave("./CARME2019/CS_incremental_PCA_with_missings/figures/MNCAR_exampleB.pdf",outgg,height=8,width=6)

