## --  Create the plot for the paper-- ##

#-- packages --#
library(ggplot2)
library(dplyr)
library(stats)
library(ggpubr)

nsim <- 5
psim <- 50*3
nmeth <- psim*nsim

res_df <-
  data.frame(method=rep(0,nmeth), sensitivity = rep(0,nmeth), specificity = rep(0,nmeth),
             MSE = rep(0,nmeth), MSEP = rep(0,nmeth), MCC = rep(0,nmeth), p=rep(0,nmeth))

pvals <- c(20,40,80,100,200)


for(i in 1:nsim){
  load(paste0("results/data_n_40p_",pvals[i],"gam_0.4q_20m_6R2_0.85.rda"))

  res_df[1:psim + psim*(i-1),1:6] <- df[,c(6, 1, 3, 4, 5, 2)]
  res_df[1:psim + psim*(i-1),7] <- pvals[i]
}
res_df[,1] <- rep(df[,6], nsim)
res_df$p <- as.factor(res_df$p)

p1<- ggplot(data = res_df, aes(x = p, y = sensitivity, fill = method)) +
  geom_boxplot()

p2<-ggplot(data = res_df, aes(x = p, y = specificity, fill = method)) +
  geom_boxplot() #+ coord_cartesian(ylim = c(0.85,1))

p3<-ggplot(data = res_df, aes(x = p, y = MSE, fill = method)) +
  geom_boxplot()

p4<-ggplot(data = res_df, aes(x = p, y = MSEP, fill = method)) +
  geom_boxplot() #+ coord_cartesian(ylim = c(0.3,0.46))

fig <- ggarrange(p1,p2,p3,p4, common.legend = T)
ggsave(filename = "results/figure.pdf", plot = fig)
