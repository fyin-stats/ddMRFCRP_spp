########################
########### summary of simulation results
#########################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("dplyr","fossil","xtable","ggplot2","fields","lattice","gridExtra")
ipak(packages)
##########################
grid <- c(20,20)
##########################
### get the rand index
### 
rlt_list <- rlt_setting1
dic_matrix <- matrix(NA, nrow=nsim, ncol = length(eta_vec))
lpml_matrix <- matrix(NA, nrow=nsim, ncol=length(eta_vec))
bic_matrix <- matrix(NA, nrow=nsim, ncol=length(eta_vec))
for(ii in 1:nsim){
  dic_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$dic))
  lpml_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$lpml))
  bic_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$bic))
}
######################
minbic_model_id <- apply(bic_matrix,1,function(x) which(x==min(x)))
mindic_model_id <- apply(dic_matrix,1,function(x) which(x==min(x)))
minlpml_model_id <- apply(lpml_matrix,1,function(x) which(x==min(x)))
######################
###########
bic_ri <- rep(NA,nsim)
dic_ri <- rep(NA,nsim)
lpml_ri <- rep(NA,nsim)
for(ii in 1:nsim){
  bic_ri[ii] <- rlt_list[[ii]][[minbic_model_id[ii]]]$ri
  dic_ri[ii] <- rlt_list[[ii]][[mindic_model_id[ii]]]$ri
  lpml_ri[ii] <- rlt_list[[ii]][[minlpml_model_id[ii]]]$ri
  # dic_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$dic))
  # lpml_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$lpml))
  # bic_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$bic))
}
#######################
bic_Khat <- rep(NA,nsim)
dic_Khat <- rep(NA,nsim)
lpml_Khat <- rep(NA,nsim)
for(ii in 1:nsim){
  bic_Khat[ii] <- rlt_list[[ii]][[minbic_model_id[ii]]]$Khat
  dic_Khat[ii] <- rlt_list[[ii]][[mindic_model_id[ii]]]$Khat
  lpml_Khat[ii] <- rlt_list[[ii]][[minlpml_model_id[ii]]]$Khat
  # dic_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$dic))
  # lpml_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$lpml))
  # bic_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$bic))
}
#######################
#######################
#######################
### get the estimated lambda
setwd("/Users/fan/Documents/dd_spp_grf/")
simulation_results_setting1 <- readRDS("./simulation_results_setting1.rds")
simulation_results_setting2 <- readRDS("./simulation_results_setting2.rds")
simulation_results_setting3 <- readRDS("./simulation_results_setting3.rds")
simulation_results_setting4 <- readRDS("./simulation_results_setting4.rds")
# simulation_results_setting4_queen <- readRDS("./simulation_results_setting4_queen.rds")
#
MFM_simulation_results_setting1 <- readRDS("./MFM_simulation_results_setting1.rds")
MFM_simulation_results_setting2 <- readRDS("./MFM_simulation_results_setting2.rds")
MFM_simulation_results_setting3 <- readRDS("./MFM_simulation_results_setting3.rds")
MFM_simulation_results_setting4 <- readRDS("./MFM_simulation_results_setting4.rds")
# 
nsim = length(simulation_results_setting1$bic_ri)
#######################
#######################
## 
# RI_summary_table <- rbind(c(0.679,0.725,0.872), c(0.804,0.805,0.846))
RI_summary_table <- rbind(c(mean(simulation_results_setting1$eta0_ri),
                            mean(simulation_results_setting1$bic_ri),
                            mean(simulation_results_setting1$dic_ri),
                            mean(simulation_results_setting1$lpml_ri),
                            mean(MFM_simulation_results_setting1$ri)), 
                          c(mean(simulation_results_setting2$eta0_ri),
                            mean(simulation_results_setting2$bic_ri),
                            mean(simulation_results_setting2$dic_ri),
                            mean(simulation_results_setting2$lpml_ri),
                            mean(MFM_simulation_results_setting2$ri)),
                          c(mean(simulation_results_setting3$eta0_ri),
                            mean(simulation_results_setting3$bic_ri),
                            mean(simulation_results_setting3$dic_ri),
                            mean(simulation_results_setting3$lpml_ri),
                            mean(MFM_simulation_results_setting3$ri)),
                          c(mean(simulation_results_setting4$eta0_ri),
                            mean(simulation_results_setting4$bic_ri),
                            mean(simulation_results_setting4$dic_ri),
                            mean(simulation_results_setting4$lpml_ri),
                            mean(MFM_simulation_results_setting4$ri)))
colnames(RI_summary_table) <- c("eta=0","BIC","DIC","LPML","MFM")
rownames(RI_summary_table) <- c("Setting 1", "Setting 2","Setting 3","Setting 4")
xtable(RI_summary_table, digits = 3)
######################
# MFM_simulation_results_setting1$ri %>% mean()
# MFM_simulation_results_setting2$ri %>% mean()
# MFM_simulation_results_setting3$ri %>% mean()
# MFM_simulation_results_setting4$ri %>% mean()
# RI_summary_table <- rbind(c(mean(simulation_results_setting1$bic_adj_ri),
#                             mean(simulation_results_setting1$dic_adj_ri),
#                             mean(simulation_results_setting1$lpml_adj_ri)), 
#                           c(mean(simulation_results_setting2$bic_adj_ri),
#                             mean(simulation_results_setting2$dic_adj_ri),
#                             mean(simulation_results_setting2$lpml_adj_ri)),
#                           c(mean(simulation_results_setting3$bic_adj_ri),
#                             mean(simulation_results_setting3$dic_adj_ri),
#                             mean(simulation_results_setting3$lpml_adj_ri)),
#                           c(mean(simulation_results_setting4$bic_adj_ri),
#                             mean(simulation_results_setting4$dic_adj_ri),
#                             mean(simulation_results_setting4$lpml_adj_ri)))
# colnames(RI_summary_table) <- c("BIC","DIC","LPML")
# rownames(RI_summary_table) <- c("setting 1", "setting 2","setting 3","setting 4")
# xtable(RI_summary_table, digits = 3)
##
# setting 1: true K = 3
# setting 2: true K = 5
# setting 3: true K = 4
# setting 4: true K = 6
Khat_summary_table <- rbind(c(mean(simulation_results_setting1$eta0_Khat==3),
                              mean(simulation_results_setting1$bic_Khat==3),
                            mean(simulation_results_setting1$dic_Khat==3),
                            mean(simulation_results_setting1$lpml_Khat==3),
                            mean(MFM_simulation_results_setting1$Khat==3)), 
                          c(mean(simulation_results_setting2$eta0_Khat==5),
                            mean(simulation_results_setting2$bic_Khat==5),
                            mean(simulation_results_setting2$dic_Khat==5),
                            mean(simulation_results_setting2$lpml_Khat==5),
                            mean(MFM_simulation_results_setting2$Khat==5)),
                          c(mean(simulation_results_setting3$eta0_Khat==4),
                            mean(simulation_results_setting3$bic_Khat==4),
                            mean(simulation_results_setting3$dic_Khat==4),
                            mean(simulation_results_setting3$lpml_Khat==4),
                            mean(MFM_simulation_results_setting3$Khat==4)),
                          c(mean(simulation_results_setting4$eta0_Khat==6),
                            mean(simulation_results_setting4$bic_Khat==6),
                            mean(simulation_results_setting4$dic_Khat==6),
                            mean(simulation_results_setting4$lpml_Khat==6),
                            mean(MFM_simulation_results_setting4$Khat==6)))
colnames(Khat_summary_table) <- c("eta=0","BIC","DIC","LPML","MFM")
rownames(Khat_summary_table) <- c("Setting 1", "Setting 2","Setting 3", "Setting 4")
xtable(Khat_summary_table, digits = 2)

##############################
### 
xtable(cbind(Khat_summary_table, RI_summary_table), digits = 3)

##### barplot
# Khat_bic <- data.frame(Khat = c(simulation_results_setting1$bic_Khat,simulation_results_setting2$bic_Khat,simulation_results_setting3$bic_Khat,
#                                 simulation_results_setting4$bic_Khat),
# setting = rep(c("setting 1", "setting 2", "setting 3", "setting 4"), rep(nsim,4) ), Ktrue = rep(c(3,5,4,6), rep(nsim,4) )  )
Khat_bic_optimal <- data.frame(Khat = c(simulation_results_setting1$bic_Khat,simulation_results_setting2$bic_Khat,
                                        simulation_results_setting3$bic_Khat),
                       setting = rep(c("Setting 1", "Setting 2", "Setting 3"), rep(nsim,3) ), Ktrue = rep(c(3,5,4), rep(nsim,3) ),
                       method = "BIC")
Khat_dic_optimal <- data.frame(Khat = c(simulation_results_setting1$dic_Khat,simulation_results_setting2$dic_Khat,simulation_results_setting3$dic_Khat),
                               setting = rep(c("Setting 1", "Setting 2", "Setting 3"), rep(nsim,3) ), Ktrue = rep(c(3,5,4), rep(nsim,3) ),
                               method = "DIC")
Khat_lpml_optimal <- data.frame(Khat = c(simulation_results_setting1$lpml_Khat,
                                         simulation_results_setting2$lpml_Khat,
                                         simulation_results_setting3$lpml_Khat),
                               setting = rep(c("Setting 1", "Setting 2", "Setting 3"), rep(nsim,3) ), Ktrue = rep(c(3,5,4), rep(nsim,3) ),
                               method = "LPML")
Khat_eta0 <- data.frame(Khat = c(simulation_results_setting1$eta0_Khat,simulation_results_setting2$eta0_Khat,
                                     simulation_results_setting3$eta0_Khat),
                       setting = rep(c("Setting 1", "Setting 2", "Setting 3"), rep(nsim,3) ), Ktrue = rep(c(3,5,4), rep(nsim,3) ),
                       method = "eta=0")
Khat_MFM <- data.frame(Khat = c(MFM_simulation_results_setting1$Khat, MFM_simulation_results_setting2$Khat,
                                MFM_simulation_results_setting3$Khat),
                        setting = rep(c("Setting 1", "Setting 2", "Setting 3"), rep(nsim,3) ), Ktrue = rep(c(3,5,4), rep(nsim,3) ),
                        method = "MFM-NHPP")
#####
Khat_df <- rbind(Khat_bic_optimal,Khat_dic_optimal, Khat_lpml_optimal, Khat_eta0, Khat_MFM)
Khat_df$method <- factor(Khat_df$method, levels = c("eta=0","BIC","DIC","LPML","MFM-NHPP"))
#####
# p<-ggplot(data=Khat_bic, aes(x=Khat)) +
#   geom_bar(stat="identity", fill="steelblue") + facet_wrap(~setting)
# p<-ggplot(data= Khat_bic %>% filter(type %in% c("BIC","eta=0")), aes(Khat, fill=type)) +
#   geom_bar() + facet_wrap(~setting+type, nrow=3) + xlab(expression(hat(K))) +
#   geom_vline(data = Khat_bic %>% filter(type %in% c("eta=0","BIC")), aes(xintercept = Ktrue), col="blue", size=0.5)
p<-ggplot(data= Khat_df, aes(Khat, fill=method)) +
  geom_bar(width=2) + facet_grid(setting~method) + xlab(expression(hat(K))) +
  geom_vline(data = Khat_df, aes(xintercept = Ktrue), col="blue", size=0.25) + theme(strip.text = element_text(size=12),
                                                                                      legend.text = element_text(size=12),
                                                                                      legend.title = element_text(size=12),
                                                                                      legend.position="bottom")
pdf("K_histogram.pdf")
p
dev.off()
########################
########################
########################
########################
########################
# est_lambda_summary_table <- rbind(c(mean(simulation_results_setting1$bic_Khat),
#                               mean(simulation_results_setting1$dic_Khat),
#                               mean(simulation_results_setting1$lpml_Khat)), 
#                             c(mean(simulation_results_setting2$bic_Khat),
#                               mean(simulation_results_setting2$dic_Khat),
#                               mean(simulation_results_setting2$lpml_Khat)),
#                             c(mean(simulation_results_setting3$bic_Khat),
#                               mean(simulation_results_setting3$dic_Khat),
#                               mean(simulation_results_setting3$lpml_Khat)))
# colnames(Khat_summary_table) <- c("BIC","DIC","LPML")
# rownames(Khat_summary_table) <- c("setting 1", "setting 2","setting 3")
# xtable(Khat_summary_table, digits = 2)
################################
true_lambda_setting1 <- readRDS("./true_lambda_setting1.rds")
true_lambda_setting2 <- readRDS("./true_lambda_setting2.rds")
true_lambda_setting3 <- readRDS("./true_lambda_setting3.rds")
true_lambda_setting4 <- readRDS("./true_lambda_setting4.rds")
### MRF mean
est_lambda_mean_bic_setting1 <- apply(simulation_results_setting1$bic_est_lambda,
                                      c(1,2), mean)
est_lambda_mean_bic_setting2 <- apply(simulation_results_setting2$bic_est_lambda,
                                      c(1,2), mean)
est_lambda_mean_bic_setting3 <- apply(simulation_results_setting3$bic_est_lambda,
                                      c(1,2), mean)
est_lambda_mean_bic_setting4 <- apply(simulation_results_setting4$bic_est_lambda,
                                      c(1,2), mean)
# median
est_lambda_median_bic_setting1 <- apply(simulation_results_setting1$bic_est_lambda,
                                      c(1,2), median)
est_lambda_median_bic_setting2 <- apply(simulation_results_setting2$bic_est_lambda,
                                      c(1,2), median)
est_lambda_median_bic_setting3 <- apply(simulation_results_setting3$bic_est_lambda,
                                      c(1,2), median)
est_lambda_median_bic_setting4 <- apply(simulation_results_setting4$bic_est_lambda,
                                      c(1,2), median)
####
## MFM mean
est_lambda_mean_MFM_setting1 <- apply(MFM_simulation_results_setting1$est_lambda,
                                      c(1,2), mean)
est_lambda_mean_MFM_setting2 <- apply(MFM_simulation_results_setting2$est_lambda,
                                      c(1,2), mean)
est_lambda_mean_MFM_setting3 <- apply(MFM_simulation_results_setting3$est_lambda,
                                      c(1,2), mean)
est_lambda_mean_MFM_setting4 <- apply(MFM_simulation_results_setting4$est_lambda,
                                      c(1,2), mean)

#
est_lambda_median_MFM_setting1 <- apply(MFM_simulation_results_setting1$est_lambda,
                                      c(1,2), median)
est_lambda_median_MFM_setting2 <- apply(MFM_simulation_results_setting2$est_lambda,
                                      c(1,2), median)
est_lambda_median_MFM_setting3 <- apply(MFM_simulation_results_setting3$est_lambda,
                                      c(1,2), median)
est_lambda_median_MFM_setting4 <- apply(MFM_simulation_results_setting4$est_lambda,
                                      c(1,2), median)

####
# relative bias
#### MRF
relative_bias_lambda_setting1 <- (est_lambda_mean_bic_setting1-true_lambda_setting1)/true_lambda_setting1
relative_bias_lambda_setting2 <- (est_lambda_mean_bic_setting2-true_lambda_setting2)/true_lambda_setting2
relative_bias_lambda_setting3 <- (est_lambda_mean_bic_setting3-true_lambda_setting3)/true_lambda_setting3
relative_bias_lambda_setting4 <- (est_lambda_mean_bic_setting4-true_lambda_setting4)/true_lambda_setting4
#### MFM
relative_bias_lambda_MFM_setting1 <- (est_lambda_mean_MFM_setting1-true_lambda_setting1)/true_lambda_setting1
relative_bias_lambda_MFM_setting2 <- (est_lambda_mean_MFM_setting2-true_lambda_setting2)/true_lambda_setting2
relative_bias_lambda_MFM_setting3 <- (est_lambda_mean_MFM_setting3-true_lambda_setting3)/true_lambda_setting3
relative_bias_lambda_MFM_setting4 <- (est_lambda_mean_MFM_setting4-true_lambda_setting4)/true_lambda_setting4

#
relative_bias_median_lambda_setting1 <- (est_lambda_mean_bic_setting1-true_lambda_setting1)/true_lambda_setting1
relative_bias_median_lambda_setting2 <- (est_lambda_mean_bic_setting2-true_lambda_setting2)/true_lambda_setting2
relative_bias_median_lambda_setting3 <- (est_lambda_mean_bic_setting3-true_lambda_setting3)/true_lambda_setting3
relative_bias_median_lambda_setting4 <- (est_lambda_mean_bic_setting4-true_lambda_setting4)/true_lambda_setting4
#### MFM
relative_bias_median_lambda_MFM_setting1 <- (est_lambda_median_MFM_setting1-true_lambda_setting1)/true_lambda_setting1
relative_bias_median_lambda_MFM_setting2 <- (est_lambda_median_MFM_setting2-true_lambda_setting2)/true_lambda_setting2
relative_bias_median_lambda_MFM_setting3 <- (est_lambda_median_MFM_setting3-true_lambda_setting3)/true_lambda_setting3
relative_bias_median_lambda_MFM_setting4 <- (est_lambda_median_MFM_setting4-true_lambda_setting4)/true_lambda_setting4

####
tool_matrix <- matrix(1, nrow = grid[1], grid[2])
#
X <- data.frame(cbind(which(tool_matrix==1, arr.ind = T),1:prod(grid)))
colnames(X) <- c("x", "y", "linear_index")
#
relative_bias_setting1 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_setting1), setting = "Setting 1")
relative_bias_setting2 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_setting2), setting = "Setting 2")
relative_bias_setting3 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_setting3), setting = "Setting 3")
relative_bias_setting4 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_setting4), setting = "Setting 4")
#
relative_bias_MFM_setting1 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_MFM_setting1), setting = "Setting 1")
relative_bias_MFM_setting2 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_MFM_setting2), setting = "Setting 2")
relative_bias_MFM_setting3 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_MFM_setting3), setting = "Setting 3")
relative_bias_MFM_setting4 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_MFM_setting4), setting = "Setting 4")

#
# relative_bias_setting1 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_median_lambda_setting1), setting = "Setting 1")
# relative_bias_setting2 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_median_lambda_setting2), setting = "Setting 2")
# relative_bias_setting3 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_median_lambda_setting3), setting = "Setting 3")
# relative_bias_setting4 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_median_lambda_setting4), setting = "Setting 4")
# #
# relative_bias_MFM_setting1 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_median_lambda_MFM_setting1), setting = "Setting 1")
# relative_bias_MFM_setting2 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_median_lambda_MFM_setting2), setting = "Setting 2")
# relative_bias_MFM_setting3 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_median_lambda_MFM_setting3), setting = "Setting 3")
# relative_bias_MFM_setting4 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_median_lambda_MFM_setting4), setting = "Setting 4")

#
relative_bias_MFM <- data.frame(rbind(relative_bias_MFM_setting1,
                                      relative_bias_MFM_setting2,
                                      relative_bias_MFM_setting3,
                                      relative_bias_MFM_setting4), method = "MFM-NHPP")
#
relative_bias_MRF <- data.frame(rbind(relative_bias_setting1,
                                      relative_bias_setting2,
                                      relative_bias_setting3,
                                      relative_bias_setting4), method = "MRF-DPM-NHPP")

#
relative_bias <- rbind(relative_bias_MRF, relative_bias_MFM)
#
relative_bias %>% filter(setting != "Setting 4") %>% ggplot(aes(x,y, fill=abs(value))) + geom_tile() + 
  scale_fill_gradient(low="white", high="red") + 
  facet_grid(method~setting)  +
  theme_bw() + labs(fill="|relative bias|") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        legend.position="bottom") + coord_fixed() # default aspect ratio = 1

##### relative bias of posterior median
# ggtitle("Relative bias")
# low="darkgoldenrod", high="darkolivegreen1"
# scale_fill_gradientn(colours = rev(terrain.colors(100)))
# tool_matrix <- matrix(1, nrow = grid[1], grid[2])
# #
# X <- data.frame(cbind(which(tool_matrix==1, arr.ind = T),1:prod(grid)))
# colnames(X) <- c("x", "y", "linear_index")
# #
# relative_bias_setting1 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_setting1), setting = "Setting 1")
# relative_bias_setting2 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_setting2), setting = "Setting 2")
# relative_bias_setting3 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_setting3), setting = "Setting 3")
# relative_bias_setting4 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_setting4), setting = "Setting 4")
# #
# relative_bias_MFM_setting1 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_MFM_setting1), setting = "Setting 1")
# relative_bias_MFM_setting2 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_MFM_setting2), setting = "Setting 2")
# relative_bias_MFM_setting3 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_MFM_setting3), setting = "Setting 3")
# relative_bias_MFM_setting4 <- data.frame(y = X$y, x = X$x, value = c(relative_bias_lambda_MFM_setting4), setting = "Setting 4")
# #
# relative_bias_MFM <- data.frame(rbind(relative_bias_MFM_setting1,
#                                       relative_bias_MFM_setting2,
#                                       relative_bias_MFM_setting3,
#                                       relative_bias_MFM_setting4), method = "MFM-NHPP")
# #
# relative_bias_MRF <- data.frame(rbind(relative_bias_setting1,
#                                       relative_bias_setting2,
#                                       relative_bias_setting3,
#                                       relative_bias_setting4), method = "MRF-DPM-NHPP")
# #
# relative_bias <- rbind(relative_bias_MRF, relative_bias_MFM)
# #
# relative_bias %>% filter(setting != "Setting 4") %>% ggplot(aes(x,y, fill=abs(value))) + geom_tile() + 
#   scale_fill_gradient(low="white", high="red") + 
#   facet_grid(method~setting)  +
#   theme_bw() + labs(fill="|relative bias|") +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(), 
#         strip.text = element_text(size=12),
#         legend.text = element_text(size=12),
#         legend.title = element_text(size=12),
#         plot.title = element_text(hjust = 0.5))
# ####
# image(x=1:grid[1],y=1:grid[2], z=relative_bias_lambda_setting1,
#       useRaster=TRUE,
#       xlab="",ylab="")
# par(mfrow=c(1,2))
# # image.plot(x=1:grid[1],y=1:grid[2], z=relative_bias_lambda_setting1,
# #       useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
# #       xlab="",ylab="")
# # image.plot(x=1:grid[1],y=1:grid[2], z=true_lambda_setting1,
# #       useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
# #       xlab="",ylab="")
# image.plot(x=1:grid[1],y=1:grid[2], z=true_lambda_setting1,
#            useRaster=TRUE,
#            xlab="",ylab="",main="True value")
# image.plot(x=1:grid[1],y=1:grid[2], z=relative_bias_lambda_setting1,
#            useRaster=TRUE,
#            xlab="",ylab="",main="Relative bias")
# ######################
# par(mfrow=c(1,2))
# image.plot(x=1:grid[1],y=1:grid[2], z=true_lambda_setting2,
#            useRaster=TRUE,
#            xlab="",ylab="",main="True value")
# image.plot(x=1:grid[1],y=1:grid[2], z=relative_bias_lambda_setting2,
#            useRaster=TRUE,
#            xlab="",ylab="",main="Relative bias")
# #######################
# par(mfrow=c(1,2))
# image.plot(x=1:grid[1],y=1:grid[2], z=true_lambda_setting3,
#            useRaster=TRUE,
#            xlab="",ylab="",main="True value")
# image.plot(x=1:grid[1],y=1:grid[2], z=relative_bias_lambda_setting3,
#            useRaster=TRUE,
#            xlab="",ylab="",main="Relative bias")
### 2.5% quantile
est_lambda_0025_bic_setting1 <- apply(simulation_results_setting1$bic_est_lambda,
                                      c(1,2), function(x) quantile(x,0.025))
est_lambda_0025_bic_setting2 <- apply(simulation_results_setting2$bic_est_lambda,
                                     c(1,2), function(x) quantile(x,0.025))
est_lambda_0025_bic_setting3 <- apply(simulation_results_setting3$bic_est_lambda,
                                     c(1,2), function(x) quantile(x,0.025))
est_lambda_0025_bic_setting4 <- apply(simulation_results_setting4$bic_est_lambda,
                                     c(1,2), function(x) quantile(x,0.025))
### 97.5% quantile
est_lambda_0975_bic_setting1 <- apply(simulation_results_setting1$bic_est_lambda,
                                     c(1,2), function(x) quantile(x,0.975))
est_lambda_0975_bic_setting2 <- apply(simulation_results_setting2$bic_est_lambda,
                                     c(1,2), function(x) quantile(x,0.975))
est_lambda_0975_bic_setting3 <- apply(simulation_results_setting3$bic_est_lambda,
                                     c(1,2), function(x) quantile(x,0.975))
est_lambda_0975_bic_setting4 <- apply(simulation_results_setting4$bic_est_lambda,
                                     c(1,2), function(x) quantile(x,0.975))
### 50% quantile (median)
est_lambda_050_bic_setting1 <- apply(simulation_results_setting1$bic_est_lambda,
                                     c(1,2), function(x) quantile(x,0.5))
est_lambda_050_bic_setting2 <- apply(simulation_results_setting2$bic_est_lambda,
                                     c(1,2), function(x) quantile(x,0.5))
est_lambda_050_bic_setting3 <- apply(simulation_results_setting3$bic_est_lambda,
                                     c(1,2), function(x) quantile(x,0.5))
est_lambda_050_bic_setting4 <- apply(simulation_results_setting4$bic_est_lambda,
                                     c(1,2), function(x) quantile(x,0.5))

### 
tool_matrix <- matrix(1, nrow = grid[1], grid[2])
#
X <- data.frame(cbind(which(tool_matrix==1, arr.ind = T),1:prod(grid)))
colnames(X) <- c("x", "y", "linear_index")
# setting 1
# image(x=1:20,y=1:20, z=true_lambda_setting1, 
#       useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
#       xlab="",ylab="")
# image(z=true_lambda_setting1), axes=FALSE)
# image(z=true_lambda_setting1)
# lambda0_lower <- colQuantiles(lambda0_sample, probs = 0.025)
# lambda0_upper <- colQuantiles(lambda0_sample, probs = 0.975)
# lambda0_median <- colQuantiles(lambda0_sample, probs = 0.5)
###

############## results for setting 1
inten0fit_dat_setting1 <- data.frame(y = X$y, x = X$x, value = c(est_lambda_0025_bic_setting1), type = "0.025 quantile")
inten0fit_dat_setting1 <- rbind(inten0fit_dat_setting1, data.frame(y = X$y, x = X$x, value = c(est_lambda_050_bic_setting1), type = "median"))
inten0fit_dat_setting1 <- rbind(inten0fit_dat_setting1, data.frame(y = X$y, x = X$x, value = c(est_lambda_0975_bic_setting1), type = "0.975 quantile"))
inten0fit_dat_setting1 <- rbind(inten0fit_dat_setting1, data.frame(y = X$y, x = X$x, value = c(true_lambda_setting1), type = "Truth"))
#############
# ggsave(p, file = "../manuscript/plots/IntenFit.pdf", height = 5, width = 4)
levelplot(value~x*y,filter(inten0fit_dat_setting1, type=="Truth"), 
          col.regions = terrain.colors(100)[length(terrain.colors(100)):1], 
          xlab = "X", ylab = "Y",
          scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
#############
breaks <- seq(0, max(c(inten0fit_dat_setting1$value))+1, length.out = 100)
cols <- terrain.colors(100)[length(terrain.colors(100)):1]
# colorRampPalette(c("grey", "yellow", "green"))(length(breaks) - 1)
######
p1_1 <- levelplot(value~x*y,filter(inten0fit_dat_setting1, type=="Truth"), 
                  col.regions = cols, at = breaks, 
                  xlab = "X", ylab = "Y", main = "Truth",
                  scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
# 
p2_1 <- levelplot(value~x*y,filter(inten0fit_dat_setting1, type=="0.025 quantile"), 
                  col.regions = cols, at = breaks, colorkey=FALSE,
                  xlab = "X", ylab = "Y", main = "2.5% quantile",
                  scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
#
p3_1 <- levelplot(value~x*y,filter(inten0fit_dat_setting1, type=="median"), 
                  col.regions = cols, at = breaks, 
                  xlab = "X", ylab = "Y", main = "Median",
                  scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
#
p4_1 <- levelplot(value~x*y,filter(inten0fit_dat_setting1, type=="0.975 quantile"), 
                  col.regions = cols, at = breaks, colorkey=FALSE,
                  xlab = "X", ylab = "Y", main = "97.5% quantile",
                  scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
#
grid.arrange(p1_1, p2_1, p3_1, p4_1, ncol=2)

############ results from setting 2
############
inten0fit_dat_setting2 <- data.frame(y = X$y, x = X$x, value = c(est_lambda_0025_bic_setting2), type = "0.025 quantile")
inten0fit_dat_setting2 <- rbind(inten0fit_dat_setting2, data.frame(y = X$y, x = X$x, value = c(est_lambda_050_bic_setting2), type = "median"))
inten0fit_dat_setting2 <- rbind(inten0fit_dat_setting2, data.frame(y = X$y, x = X$x, value = c(est_lambda_0975_bic_setting2), type = "0.975 quantile"))
inten0fit_dat_setting2 <- rbind(inten0fit_dat_setting2, data.frame(y = X$y, x = X$x, value = c(true_lambda_setting2), type = "Truth"))
#############
# ggsave(p, file = "../manuscript/plots/IntenFit.pdf", height = 5, width = 4)
levelplot(value~x*y,filter(inten0fit_dat_setting2, type=="Truth"), 
          col.regions = terrain.colors(100)[length(terrain.colors(100)):1], 
          xlab = "X", ylab = "Y",
          scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
##############
#############
breaks <- seq(0, max(c(inten0fit_dat_setting2$value))+1, length.out = 100)
cols <- terrain.colors(100)[length(terrain.colors(100)):1]
# colorRampPalette(c("grey", "yellow", "green"))(length(breaks) - 1)
######
p1_2 <- levelplot(value~x*y,filter(inten0fit_dat_setting2, type=="Truth"), 
                  col.regions = cols, at = breaks, 
                  xlab = "X", ylab = "Y", main = "Truth",
                  scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
# 
p2_2 <- levelplot(value~x*y,filter(inten0fit_dat_setting2, type=="0.025 quantile"), 
                  col.regions = cols, at = breaks, colorkey=FALSE,
                  xlab = "X", ylab = "Y", main = "2.5% quantile",
                  scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
#
p3_2 <- levelplot(value~x*y,filter(inten0fit_dat_setting2, type=="median"), 
                  col.regions = cols, at = breaks, 
                  xlab = "X", ylab = "Y", main = "Median",
                  scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
#
p4_2 <- levelplot(value~x*y,filter(inten0fit_dat_setting2, type=="0.975 quantile"), 
                  col.regions = cols, at = breaks, colorkey=FALSE,
                  xlab = "X", ylab = "Y", main = "97.5% quantile",
                  scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
#
grid.arrange(p1_2, p2_2, p3_2, p4_2, ncol=2)

##########################
### results from setting 3
##########################
inten0fit_dat_setting3 <- data.frame(y = X$y, x = X$x, value = c(est_lambda_0025_bic_setting3), type = "0.025 quantile")
inten0fit_dat_setting3 <- rbind(inten0fit_dat_setting3, data.frame(y = X$y, x = X$x, value = c(est_lambda_050_bic_setting3), type = "median"))
inten0fit_dat_setting3 <- rbind(inten0fit_dat_setting3, data.frame(y = X$y, x = X$x, value = c(est_lambda_0975_bic_setting3), type = "0.975 quantile"))
inten0fit_dat_setting3 <- rbind(inten0fit_dat_setting3, data.frame(y = X$y, x = X$x, value = c(true_lambda_setting3), type = "Truth"))
#############
# ggsave(p, file = "../manuscript/plots/IntenFit.pdf", height = 5, width = 4)
# par(mfrow=c(1,4))
breaks <- seq(0, max(c(inten0fit_dat_setting3$value))+1, length.out = 100)
cols <- terrain.colors(100)[length(terrain.colors(100)):1]
# colorRampPalette(c("grey", "yellow", "green"))(length(breaks) - 1)
######
p1_3 <- levelplot(value~x*y,filter(inten0fit_dat_setting3, type=="Truth"), 
          col.regions = cols, at = breaks, 
          xlab = "X", ylab = "Y", main = "Truth",
          scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
# 
p2_3 <- levelplot(value~x*y,filter(inten0fit_dat_setting3, type=="0.025 quantile"), 
          col.regions = cols, at = breaks, colorkey=FALSE,
          xlab = "X", ylab = "Y", main = "2.5% quantile",
          scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
#
p3_3 <- levelplot(value~x*y,filter(inten0fit_dat_setting3, type=="median"), 
          col.regions = cols, at = breaks, 
          xlab = "X", ylab = "Y", main = "Median",
          scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
#
p4_3 <- levelplot(value~x*y,filter(inten0fit_dat_setting3, type=="0.975 quantile"), 
          col.regions = cols, at = breaks, colorkey=FALSE,
          xlab = "X", ylab = "Y", main = "97.5% quantile",
          scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
#
grid.arrange(p1_3, p2_3, p3_3, p4_3, ncol=2)
#############################
#############################
### results from setting 4
##########################
inten0fit_dat_setting4 <- data.frame(y = X$y, x = X$x, value = c(est_lambda_0025_bic_setting4), type = "0.025 quantile")
inten0fit_dat_setting4 <- rbind(inten0fit_dat_setting4, data.frame(y = X$y, x = X$x, value = c(est_lambda_050_bic_setting4), type = "median"))
inten0fit_dat_setting4 <- rbind(inten0fit_dat_setting4, data.frame(y = X$y, x = X$x, value = c(est_lambda_0975_bic_setting4), type = "0.975 quantile"))
inten0fit_dat_setting4 <- rbind(inten0fit_dat_setting4, data.frame(y = X$y, x = X$x, value = c(true_lambda_setting4), type = "Truth"))
#############
# ggsave(p, file = "../manuscript/plots/IntenFit.pdf", height = 5, width = 4)
breaks <- seq(0, max(c(inten0fit_dat_setting4$value))+1, length.out = 100)
cols <- terrain.colors(100)[length(terrain.colors(100)):1]
# colorRampPalette(c("grey", "yellow", "green"))(length(breaks) - 1)
######
p1_4 <- levelplot(value~x*y,filter(inten0fit_dat_setting4, type=="Truth"), 
                  col.regions = cols, at = breaks, 
                  xlab = "X", ylab = "Y", main = "Truth",
                  scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
# 
p2_4 <- levelplot(value~x*y,filter(inten0fit_dat_setting4, type=="0.025 quantile"), 
                  col.regions = cols, at = breaks, colorkey=FALSE,
                  xlab = "X", ylab = "Y", main = "2.5% quantile",
                  scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
#
p3_4 <- levelplot(value~x*y,filter(inten0fit_dat_setting4, type=="median"), 
                  col.regions = cols, at = breaks, 
                  xlab = "X", ylab = "Y", main = "Median",
                  scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
#
p4_4 <- levelplot(value~x*y,filter(inten0fit_dat_setting4, type=="0.975 quantile"), 
                  col.regions = cols, at = breaks, colorkey=FALSE,
                  xlab = "X", ylab = "Y", main = "97.5% quantile",
                  scales = list(x = list(at = seq(0, 20, length.out = 5)), y = list(at = seq(0, 20, length.out = 5))))
#
grid.arrange(p1_4, p2_4, p3_4, p4_4, ncol=2)
##############################
##############################

#### put them together using ggplot2
inten0fit_dat_MRF <- data.frame(rbind(data.frame(inten0fit_dat_setting1, setting = "Setting 1"),
                                      data.frame(inten0fit_dat_setting2, setting = "Setting 2"),
                                      data.frame(inten0fit_dat_setting3, setting = "Setting 3"),
                                      data.frame(inten0fit_dat_setting4, setting = "Setting 4")), 
                                method = "MRF-DPM-NHPP, BIC")
#
inten0fit_dat_MRF$type <- factor(inten0fit_dat_MRF$type, levels = c("Truth",
                                                                    "0.025 quantile",
                                                                    "median",
                                                                    "0.975 quantile"),
                                 labels = c("Truth", "2.5% quantile", "Median",
                                            "97.5% quantile"))

# plot using geom_tile
inten0fit_dat_MRF %>% filter(setting != "Setting 4") %>% ggplot(aes(x,y, fill=value)) + geom_tile() + 
  scale_fill_gradientn(colours = rev(terrain.colors(100))) +
  facet_grid(setting~type) +
  theme_bw() + labs(fill=expression(lambda)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        strip.text = element_text(size=12),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position="bottom", plot.margin=grid::unit(c(0,0,0,0), "mm")) + coord_fixed()

# options(repr.plot.width = 1.25, repr.plot.height = 1.25)
#   scale_fill_gradient(low="darkgoldenrod", high="blueviolet") + 


# ggplot(data, aes(X, Y, fill= Z)) +
#   geom_tile() + theme_bw() +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())





####################################################
### MFM-NHPP intensity surface
####################################################
true_lambda_setting1 <- readRDS("./true_lambda_setting1.rds")
true_lambda_setting2 <- readRDS("./true_lambda_setting2.rds")
true_lambda_setting3 <- readRDS("./true_lambda_setting3.rds")
true_lambda_setting4 <- readRDS("./true_lambda_setting4.rds")
### mean
est_lambda_mean_MFM_setting1 <- apply(MFM_simulation_results_setting1$est_lambda,
                                      c(1,2), mean)
est_lambda_mean_MFM_setting2 <- apply(MFM_simulation_results_setting2$est_lambda,
                                      c(1,2), mean)
est_lambda_mean_MFM_setting3 <- apply(MFM_simulation_results_setting3$est_lambda,
                                      c(1,2), mean)
est_lambda_mean_MFM_setting4 <- apply(MFM_simulation_results_setting4$est_lambda,
                                      c(1,2), mean)

### 
### 2.5% quantile
est_lambda_0025_MFM_setting1 <- apply(MFM_simulation_results_setting1$est_lambda,
                                      c(1,2), function(x) quantile(x,0.025))
est_lambda_0025_MFM_setting2 <- apply(MFM_simulation_results_setting2$est_lambda,
                                      c(1,2), function(x) quantile(x,0.025))
est_lambda_0025_MFM_setting3 <- apply(MFM_simulation_results_setting3$est_lambda,
                                      c(1,2), function(x) quantile(x,0.025))
est_lambda_0025_MFM_setting4 <- apply(MFM_simulation_results_setting4$est_lambda,
                                      c(1,2), function(x) quantile(x,0.025))
### 97.5% quantile
est_lambda_0975_MFM_setting1 <- apply(MFM_simulation_results_setting1$est_lambda,
                                      c(1,2), function(x) quantile(x,0.975))
est_lambda_0975_MFM_setting2 <- apply(MFM_simulation_results_setting2$est_lambda,
                                      c(1,2), function(x) quantile(x,0.975))
est_lambda_0975_MFM_setting3 <- apply(MFM_simulation_results_setting3$est_lambda,
                                      c(1,2), function(x) quantile(x,0.975))
est_lambda_0975_MFM_setting4 <- apply(MFM_simulation_results_setting4$est_lambda,
                                      c(1,2), function(x) quantile(x,0.975))
### 50% quantile (median)
est_lambda_050_MFM_setting1 <- apply(MFM_simulation_results_setting1$est_lambda,
                                     c(1,2), function(x) quantile(x,0.5))
est_lambda_050_MFM_setting2 <- apply(MFM_simulation_results_setting2$est_lambda,
                                     c(1,2), function(x) quantile(x,0.5))
est_lambda_050_MFM_setting3 <- apply(MFM_simulation_results_setting3$est_lambda,
                                     c(1,2), function(x) quantile(x,0.5))
est_lambda_050_MFM_setting4 <- apply(MFM_simulation_results_setting4$est_lambda,
                                     c(1,2), function(x) quantile(x,0.5))

#### setting 1
inten0fit_dat_MFM_setting1 <- data.frame(y = X$y, x = X$x, value = c(est_lambda_0025_MFM_setting1), type = "0.025 quantile")
inten0fit_dat_MFM_setting1 <- rbind(inten0fit_dat_setting1, data.frame(y = X$y, x = X$x, value = c(est_lambda_050_MFM_setting1), type = "median"))
inten0fit_dat_MFM_setting1 <- rbind(inten0fit_dat_setting1, data.frame(y = X$y, x = X$x, value = c(est_lambda_0975_MFM_setting1), type = "0.975 quantile"))
inten0fit_dat_MFM_setting1 <- rbind(inten0fit_dat_setting1, data.frame(y = X$y, x = X$x, value = c(true_lambda_setting1), type = "Truth"))

#### setting 2
inten0fit_dat_MFM_setting2 <- data.frame(y = X$y, x = X$x, value = c(est_lambda_0025_MFM_setting2), type = "0.025 quantile")
inten0fit_dat_MFM_setting2 <- rbind(inten0fit_dat_setting2, data.frame(y = X$y, x = X$x, value = c(est_lambda_050_MFM_setting2), type = "median"))
inten0fit_dat_MFM_setting2 <- rbind(inten0fit_dat_setting2, data.frame(y = X$y, x = X$x, value = c(est_lambda_0975_MFM_setting2), type = "0.975 quantile"))
inten0fit_dat_MFM_setting2 <- rbind(inten0fit_dat_setting2, data.frame(y = X$y, x = X$x, value = c(true_lambda_setting2), type = "Truth"))

#### setting 3
inten0fit_dat_MFM_setting3 <- data.frame(y = X$y, x = X$x, value = c(est_lambda_0025_MFM_setting3), type = "0.025 quantile")
inten0fit_dat_MFM_setting3 <- rbind(inten0fit_dat_setting3, data.frame(y = X$y, x = X$x, value = c(est_lambda_050_MFM_setting3), type = "median"))
inten0fit_dat_MFM_setting3 <- rbind(inten0fit_dat_setting3, data.frame(y = X$y, x = X$x, value = c(est_lambda_0975_MFM_setting3), type = "0.975 quantile"))
inten0fit_dat_MFM_setting3 <- rbind(inten0fit_dat_setting3, data.frame(y = X$y, x = X$x, value = c(true_lambda_setting3), type = "Truth"))

#### setting 4
inten0fit_dat_MFM_setting4 <- data.frame(y = X$y, x = X$x, value = c(est_lambda_0025_MFM_setting4), type = "0.025 quantile")
inten0fit_dat_MFM_setting4 <- rbind(inten0fit_dat_setting4, data.frame(y = X$y, x = X$x, value = c(est_lambda_050_MFM_setting4), type = "median"))
inten0fit_dat_MFM_setting4 <- rbind(inten0fit_dat_setting4, data.frame(y = X$y, x = X$x, value = c(est_lambda_0975_MFM_setting4), type = "0.975 quantile"))
inten0fit_dat_MFM_setting4 <- rbind(inten0fit_dat_setting4, data.frame(y = X$y, x = X$x, value = c(true_lambda_setting4), type = "Truth"))

####
inten0fit_dat_MFM <- data.frame(rbind(data.frame(inten0fit_dat_MFM_setting1, setting = "Setting 1"),
                                      data.frame(inten0fit_dat_MFM_setting2, setting = "Setting 2"),
                                      data.frame(inten0fit_dat_MFM_setting3, setting = "Setting 3"),
                                      data.frame(inten0fit_dat_MFM_setting4, setting = "Setting 4")), 
                                method = "MFM-NHPP")
#
inten0fit_dat_MFM$type <- factor(inten0fit_dat_MFM$type, levels = c("Truth",
                                                                    "0.025 quantile",
                                                                    "median",
                                                                    "0.975 quantile"),
                                 labels = c("Truth", "2.5% quantile", "Median",
                                            "97.5% quantile"))

#
inten0fit_dat_MFM %>% filter(setting != "Setting 4") %>% ggplot(aes(x,y, fill=value)) + geom_tile() + 
  scale_fill_gradientn(colours = rev(terrain.colors(100))) + 
  facet_grid(setting~type) +
  theme_bw()  +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position="bottom") + coord_fixed(ratio = 1)



########### visualize their difference
















##############################
# relative bias
# relative_bias_lambda_setting1 <- (est_lambda_mean_bic_setting1-true_lambda_setting1)/true_lambda_setting1
# relative_bias_lambda_setting2 <- (est_lambda_mean_bic_setting2-true_lambda_setting2)/true_lambda_setting2
# relative_bias_lambda_setting3 <- (est_lambda_mean_bic_setting3-true_lambda_setting3)/true_lambda_setting3
# relative_bias_lambda_setting4 <- (est_lambda_mean_bic_setting4-true_lambda_setting4)/true_lambda_setting4
# # relative bias lambda 
# relative_bias_lambda_df <- data.frame(y = X$y, x = X$x, value = c(est_lambda_002_bic_setting4), 
#                                       type = "0.025 quantile")

# true values
################################
########### traceplot for RI
################################
ri_trace_df_setting1 <- data.frame(rbind(simulation_results_setting1$ri_trace_df,
                                         MFM_simulation_results_setting1$ri_trace_df), setting = "Setting 1")
ri_trace_df_setting2 <- data.frame(rbind(simulation_results_setting2$ri_trace_df,
                                         MFM_simulation_results_setting2$ri_trace_df), setting = "Setting 2")
ri_trace_df_setting3 <- data.frame(rbind(simulation_results_setting3$ri_trace_df,
                                         MFM_simulation_results_setting3$ri_trace_df), setting = "Setting 3")
ri_trace_df_setting4 <- data.frame(rbind(simulation_results_setting4$ri_trace_df,
                                         MFM_simulation_results_setting4$ri_trace_df), setting = "Setting 4")
###############################
# ri_trace_df <- rbind(ri_trace_df_setting1, ri_trace_df_setting2, ri_trace_df_setting3, ri_trace_df_setting4)
ri_trace_df <- rbind(ri_trace_df_setting1, ri_trace_df_setting2, ri_trace_df_setting3)
ri_trace_df$method[which(ri_trace_df$method == "MFM")] <- "MFM-NHPP"
ri_trace_df$method <- factor(ri_trace_df$method, levels = c("eta=0","BIC","DIC","LPML","MFM-NHPP"))
ri_trace_df_summary <- data.frame(ri_trace_df %>% group_by(method, setting, iter) %>% 
                                    summarise(mean_RI = mean(RI)))
#

###

p_ri_trace <- ggplot(data=ri_trace_df, aes(x=iter, y=RI)) + geom_line(colour="darkgrey", size=0.25) + 
  facet_grid(setting~method) + geom_line(data=ri_trace_df_summary, aes(x=iter,y=mean_RI), colour="red", size=0.5) + theme(strip.text = element_text(size=12),
                                                                                                                                 legend.text = element_text(size=12),
                                                                                                                                 legend.title = element_text(size=12))

p_ri_trace_scalefree <- ggplot(data=ri_trace_df, aes(x=iter, y=RI)) + geom_line(colour="darkgrey", size=0.25) + 
  facet_grid(setting~method, scales="free_y") + 
  geom_line(data=ri_trace_df_summary, aes(x=iter,y=mean_RI), colour="red", size=0.5) + 
  theme(strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))


# RI.rep=matrix(0,1000,totalrep)
# for(j in 1:totalrep){
#   for(i in 1:1000){
#     RI.rep[i,j]=fossil::rand.index(asm,fit.rep[[j]]$Iterates[[i]]$zout)
#   }
# }
# 
# df <- data.frame(method=c(rep("MFM",1000)),
#                  iteration=c(1:1000),
#                  RI=apply(RI.rep,1,mean))
# 
# RIplot=ggplot(data=df, aes(x=iteration, y=RI, group=method)) +
#   geom_line(aes(color=method))+
#   labs(x="Iteration", y = "Rand Index")+ theme(legend.position = c(0.8, 0.2))
# 
# 
# df2 <- data.frame(method=c(rep("REP",1000*50)),
#                   iteration=c(rep(1:1000,50)),
#                   RI=c(RI.rep))
# RIplot <- ggplot(data=df2, aes(x=iteration, y=RI, group=method)) +
#   geom_line(aes(color=method),colour="darkgrey")+
#   labs(x="Iteration", y = "Rand Index")+ 
#   geom_line(data = df,aes(x = iteration,y = RI,group = method),colour = "red")+
#   theme(legend.position = c(0.8, 0.2),
#         text = element_text(size = 16)) + theme_bw()
# 
# ggsave(plot = RIplot, width = 6, height = 3,
#        filename = "~/Desktop/simtrace.pdf")

##############################
#### summary of tuning parameters
##############################