rm(list=ls(all=T))
load("R/pdf_update_result.RData")
str(update_out)

lowcorr=ggcorrplot(cor(update_out$Xlist[[6]],use = "pairwise.complete.obs"))
ggsave("~/Dropbox/idm_Stanford/CARME2019/weak_corr.pdf",lowcorr,height=8,width=8)


library('tidyverse')
library('MASS')
library('xtable')

set.seed(1234)
r=.75
samples=10

data = as_tibble(round(mvrnorm(n=samples, mu=c(0, 0), Sigma=matrix(c(1, r, r, 1), nrow=2), empirical=TRUE),2))
names(data) = c("Xa","Xb")
cor(data)

my_data=data %>% 
  mutate(Xab=Xa * Xb,
         Xa_miss= replace(Xa,Xab > quantile(Xab,.75),NA),
         Q3=quantile(Xab,.75)
  )

my_data %>%
  select(Xa_miss, Xb)%>%
  cor(use="complete")

  
  xtable(
    my_data %>% select(-Q3)
    
  )
  


