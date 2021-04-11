################### 
################### simulation studies
###################
###################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("spatstat","doParallel","MHadaptive","purrr","dplyr")
ipak(packages)
# library(spatstat)
# library(doParallel)
####################
# # uniform Poisson process with intensity 1 in a 10 x 10 square
# pp <- rpoispp(1, win=owin(c(0,10),c(0,10)))
# pp <- rpoispp(1, win=owin(c(0,10),c(0,10)))
source("./code/ddMRFCRP_gibbs.R")
###################
### setting 1: three clusters
### (0.2,5,20)
### 10*9: 1 (lambda = 12)
### 5*18: 1 (lambda = 6)
### irregular: 1 (lambda = 0.2)
grid <- c(20,20)
eta_vec <- seq(0,8,by=0.5)
# gamma prior
a <- 1
b <- 1
# dirichlet prior
alpha <- 1
# initial number of clusters
initNCluster <- 5
niterations <- 4000
burnin <- 2000
thin <- 10
##########
# generate inhomogenous point patterns according to Z
pp_region_list <- readRDS("./pp_region_list_setting1.rds")
true_cluster_membership_setting <- readRDS("./true_cluster_membership_setting1.rds")
##########
nsim <- 100
###########
time1 <- Sys.time()
registerDoParallel(cores=50)
rlt_setting1 <- foreach(ii = 1:nsim, .errorhandling="pass") %dopar%{
  # generate point patterns
  # set.seed(ii)
  # pp <- rpoispp(Z)
  # generate region number
  pp_region <- pp_region_list[[ii]] 
  # pp2region(pp,nrow=grid[1],ncol=grid[2])
  results <- vector("list", length = length(eta_vec))
  # run the algorithm under different values of inverse temperature (eta)
  for(k in 1:length(eta_vec)){
      sample <- pp_Gibbs_ddMRFCRP(n = prod(grid), index_grid = pp_region$pp_region_number,
                               a = a, b = b, alpha = alpha, A = 1, eta = eta_vec[k],
                               neighbor_list = pp_region$neighbor_list,
                               initNCluster = initNCluster,
                               niterations = niterations)
      # n, index_grid, a, b, eta, alpha, A, neighbor_list,
      # initNCluster, niterations
      post_sample <- sample$Iterates[seq(from = (burnin+1), to = niterations, by = thin)]
      # ri, adj_ri trace (after burn-in)
      ri_trace <- rep(NA, length(sample$Iterates))
      adj_ri_trace <- rep(NA, length(sample$Iterates))
      for(jj in 1:length(ri_trace)){
        ri_trace[jj] <- fossil::rand.index(true_cluster_membership_setting, sample$Iterates[[jj]]$index_cluster)
        adj_ri_trace[jj] <- fossil::adj.rand.index(true_cluster_membership_setting, sample$Iterates[[jj]]$index_cluster)
      }
      ######
      Dahl_index <- dahl(post_sample)$cLS
      z_fitted <- post_sample[[Dahl_index]]$index_cluster
      # ri <- RI(z_fitted, sample$Iterates)
      ri <- fossil::rand.index(true_cluster_membership_setting, z_fitted)
      adj_ri <- fossil::adj.rand.index(true_cluster_membership_setting, z_fitted)
      #
      bic <- BIC(samples = post_sample,
                 Dahl_index = Dahl_index, index_grid = pp_region$pp_region_number, A = 1)

      lpml <- LPML(samples = post_sample,
                   index_grid = pp_region$pp_region_number, A = 1)
      dic <- DIC(samples = post_sample,
                 index_grid = pp_region$pp_region_number, grid = grid)$DIC
      # results[[k]] <- list(sample = sample, simudata = simudata, ri = ri, Dahl_index = Dahl_index,
      #                      bic = bic, lpml = lpml, dic = dic)
      results[[k]] <- list(ri = ri, adj_ri = adj_ri, Dahl_index = Dahl_index, Khat = max(z_fitted),
                           z_fitted=z_fitted, post_sample = post_sample, sample = sample,
                           bic = bic, lpml = -lpml, dic = dic, ri_trace = ri_trace, adj_ri_trace = adj_ri_trace)
  }
  results
}
#####################
time2 <- Sys.time()
######################
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
eta0_ri <- rep(NA, nsim)
###########
bic_adj_ri <- rep(NA,nsim)
dic_adj_ri <- rep(NA,nsim)
lpml_adj_ri <- rep(NA,nsim)
eta0_adj_ri <- rep(NA, nsim)
for(ii in 1:nsim){
  bic_ri[ii] <- rlt_list[[ii]][[minbic_model_id[ii]]]$ri
  dic_ri[ii] <- rlt_list[[ii]][[mindic_model_id[ii]]]$ri
  lpml_ri[ii] <- rlt_list[[ii]][[minlpml_model_id[ii]]]$ri
  eta0_ri[ii] <- rlt_list[[ii]][[1]]$ri
  #
  bic_adj_ri[ii] <- rlt_list[[ii]][[minbic_model_id[ii]]]$adj_ri
  dic_adj_ri[ii] <- rlt_list[[ii]][[mindic_model_id[ii]]]$adj_ri
  lpml_adj_ri[ii] <- rlt_list[[ii]][[minlpml_model_id[ii]]]$adj_ri
  eta0_adj_ri[ii] <- rlt_list[[ii]][[1]]$adj_ri
  # dic_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$dic))
  # lpml_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$lpml))
  # bic_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$bic))
}
#######################
##########
bic_Khat <- rep(NA,nsim)
dic_Khat <- rep(NA,nsim)
lpml_Khat <- rep(NA,nsim)
eta0_Khat <- rep(NA,nsim)
for(ii in 1:nsim){
  bic_Khat[ii] <- rlt_list[[ii]][[minbic_model_id[ii]]]$Khat
  dic_Khat[ii] <- rlt_list[[ii]][[mindic_model_id[ii]]]$Khat
  lpml_Khat[ii] <- rlt_list[[ii]][[minlpml_model_id[ii]]]$Khat
  eta0_Khat[ii] <- rlt_list[[ii]][[1]]$Khat
  # dic_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$dic))
  # lpml_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$lpml))
  # bic_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$bic))
}

#####################
##### lambda estimates
#####################
#### stored in an array (20*20*50) # 50 reps
######################
#####################
# rlt_list[[1]][[1]]$post_sample[[rlt_list[[1]][[1]]$Dahl_index]]$lambda
# bic_est_lambda
bic_est_lambda <- array(NA, dim=c(grid,nsim))
dic_est_lambda <- array(NA, dim=c(grid,nsim))
lpml_est_lambda <- array(NA, dim=c(grid,nsim))
# dic_ri <- rep(NA,nsim)
# lpml_ri <- rep(NA,nsim)
for(ii in 1:nsim){
  # 
  # figure out z_fitted
  # rlt_list[[ii]][[ minbic_model_id[ii] ]]$z_fitted
  bic_est_lambda[,,ii] <- matrix(rlt_list[[ii]][[ minbic_model_id[ii] ]]$post_sample[[ rlt_list[[ii]][[minbic_model_id[ii]]]$Dahl_index ]]$lambda[rlt_list[[ii]][[ minbic_model_id[ii] ]]$z_fitted],
                                  nrow=20,ncol=20)
  dic_est_lambda[,,ii] <- matrix(rlt_list[[ii]][[ mindic_model_id[ii] ]]$post_sample[[ rlt_list[[ii]][[mindic_model_id[ii]]]$Dahl_index ]]$lambda[rlt_list[[ii]][[ mindic_model_id[ii] ]]$z_fitted],
                                 nrow=20,ncol=20)
  lpml_est_lambda[,,ii] <- matrix(rlt_list[[ii]][[ minlpml_model_id[ii] ]]$post_sample[[ rlt_list[[ii]][[minlpml_model_id[ii]]]$Dahl_index ]]$lambda[rlt_list[[ii]][[ minlpml_model_id[ii] ]]$z_fitted],
                                 nrow=20,ncol=20)

  # bic_ri[ii] <- rlt_list[[ii]][[minbic_model_id[ii]]]$ri
}
######################
######################
######################
# ri trace, adj ri trace
### data frame, columns: iterations, RI (adj RI), method
ri_trace_df_eta0 <- data.frame(iter=NA, RI = NA, method = NA)
ri_trace_df_bic <- data.frame(iter=NA, RI = NA, method = NA)
ri_trace_df_dic <- data.frame(iter=NA, RI = NA, method = NA)
ri_trace_df_lpml <- data.frame(iter=NA, RI = NA, method = NA)
#
adj_ri_trace_df_eta0 <- data.frame(iter=NA, RI = NA, method = NA)
adj_ri_trace_df_bic <- data.frame(iter=NA, RI = NA, method = NA)
adj_ri_trace_df_dic <- data.frame(iter=NA, RI = NA, method = NA)
adj_ri_trace_df_lpml <- data.frame(iter=NA, RI = NA, method = NA)
#
for(ii in 1:nsim){
  ri_trace_df_eta0 <- rbind(ri_trace_df_eta0,
                            data.frame(iter=1:length(rlt_list[[ii]][[1]]$ri_trace),
                                       RI=rlt_list[[ii]][[1]]$ri_trace,
                                       method=rep("eta=0",length(rlt_list[[ii]][[1]]$ri_trace))))
  
  ri_trace_df_bic <- rbind(ri_trace_df_bic,
                           data.frame(iter=1:length(rlt_list[[ii]][[ minbic_model_id[ii] ]]$ri_trace),
                                      RI=rlt_list[[ii]][[ minbic_model_id[ii] ]]$ri_trace,
                                      method=rep("BIC",length(rlt_list[[ii]][[ minbic_model_id[ii] ]]$ri_trace))))
  
  ri_trace_df_dic <- rbind(ri_trace_df_dic,
                           data.frame(iter=1:length(rlt_list[[ii]][[ mindic_model_id[ii] ]]$ri_trace),
                                      RI=rlt_list[[ii]][[ mindic_model_id[ii] ]]$ri_trace,
                                      method=rep("DIC",length(rlt_list[[ii]][[ mindic_model_id[ii] ]]$ri_trace))))
  
  ri_trace_df_lpml <- rbind(ri_trace_df_lpml,
                            data.frame(iter=1:length(rlt_list[[ii]][[ minlpml_model_id[ii] ]]$ri_trace),
                                       RI=rlt_list[[ii]][[ minlpml_model_id[ii] ]]$ri_trace,
                                       method=rep("LPML",length(rlt_list[[ii]][[ minlpml_model_id[ii] ]]$ri_trace))))
  
  # adj_ri
  adj_ri_trace_df_eta0 <- rbind(adj_ri_trace_df_eta0,
                                data.frame(iter=1:length(rlt_list[[ii]][[1]]$adj_ri_trace),
                                           RI=rlt_list[[ii]][[1]]$adj_ri_trace,
                                           method=rep("eta=0",length(rlt_list[[ii]][[1]]$adj_ri_trace))))
  
  adj_ri_trace_df_bic <- rbind(adj_ri_trace_df_bic,
                               data.frame(iter=1:length(rlt_list[[ii]][[ minbic_model_id[ii] ]]$adj_ri_trace),
                                          RI=rlt_list[[ii]][[ minbic_model_id[ii] ]]$adj_ri_trace,
                                          method=rep("BIC",length(rlt_list[[ii]][[ minbic_model_id[ii] ]]$adj_ri_trace))))
  
  adj_ri_trace_df_dic <- rbind(adj_ri_trace_df_dic,
                               data.frame(iter=1:length(rlt_list[[ii]][[ mindic_model_id[ii] ]]$adj_ri_trace),
                                          RI=rlt_list[[ii]][[ mindic_model_id[ii] ]]$adj_ri_trace,
                                          method=rep("DIC",length(rlt_list[[ii]][[ mindic_model_id[ii] ]]$adj_ri_trace))))
  
  adj_ri_trace_df_lpml <- rbind(adj_ri_trace_df_lpml,
                                data.frame(iter=1:length(rlt_list[[ii]][[ minlpml_model_id[ii] ]]$adj_ri_trace),
                                           RI=rlt_list[[ii]][[ minlpml_model_id[ii] ]]$adj_ri_trace,
                                           method=rep("LPML",length(rlt_list[[ii]][[ minlpml_model_id[ii] ]]$adj_ri_trace))))
  
}
##### combine them together
ri_trace_df <- rbind(ri_trace_df_eta0, ri_trace_df_bic, ri_trace_df_dic, ri_trace_df_lpml)
adj_ri_trace_df <- rbind(adj_ri_trace_df_eta0, adj_ri_trace_df_bic, adj_ri_trace_df_dic, adj_ri_trace_df_lpml)
###### remove the NA's
ri_trace_df <- ri_trace_df[complete.cases(ri_trace_df),]
adj_ri_trace_df <- ri_trace_df[complete.cases(adj_ri_trace_df),]

# save the results: Khat, RI, est_lambda
######################
######################
simulation_results <- list(bic_ri=bic_ri,
                           dic_ri=dic_ri,
                           lpml_ri=lpml_ri,
                           eta0_ri = eta0_ri,
                           bic_est_lambda=bic_est_lambda,
                           dic_est_lambda=dic_est_lambda,
                           lpml_est_lambda=lpml_est_lambda,
                           bic_Khat=bic_Khat,
                           dic_Khat=dic_Khat,
                           lpml_Khat=lpml_Khat,
                           eta0_Khat = eta0_Khat,
                           minbic_model_id = minbic_model_id,
                           mindic_model_id = mindic_model_id,
                           minlpml_model_id = minlpml_model_id,
                           ri_trace_df = ri_trace_df,
                           adj_ri_trace_df = adj_ri_trace_df)
######################
saveRDS(simulation_results, "simulation_results_setting1.rds")
#####################
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