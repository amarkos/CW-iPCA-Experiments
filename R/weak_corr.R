rm(list=ls(all=T))
load("R/pdf_update_result.RData")
library("corrplot")
library("corrr")
str(update_out)

pdf(file="CARME2019/CS_incremental_PCA_with_missings/figures/weak_corr.pdf",lowcorr)
corrplot(cor(as.data.frame(update_out$Xlist[[5]]),use = "pairwise.complete.obs"))
dev.off()

pdf(file="CARME2019/CS_incremental_PCA_with_missings/figures/strong_corr.pdf",lowcorr)
corrplot(cor(as.data.frame(update_out$no_naive$Xcom_list[[5]]),use = "pairwise.complete.obs"))
dev.off()
