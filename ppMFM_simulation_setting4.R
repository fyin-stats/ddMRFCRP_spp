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
source("./code/MFM_example.R")
###################
### setting 1: three clusters
### (0.2,5,20)
### 10*9: 1 (lambda = 12)
### 5*18: 1 (lambda = 6)
### irregular: 1 (lambda = 0.2)
grid <- c(20,20)
n <- prod(grid)
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
pp_region_list <- readRDS("./pp_region_list_setting4.rds")
true_cluster_membership_setting <- readRDS("./true_cluster_membership_setting4.rds")
##########
nsim <- 100
###########
time1 <- Sys.time()
registerDoParallel(cores=50)
rlt_setting4 <- foreach(ii = 1:nsim, .errorhandling="pass") %dopar%{
  # generate point patterns
  # set.seed(ii)
  # pp <- rpoispp(Z)
  # generate region number
  pp_region <- pp_region_list[[ii]] 
  # pp2region(pp,nrow=grid[1],ncol=grid[2])
  # results <- vector("list", length = length(eta_vec))
  
  # data preprocess
  index_grid <- pp_region$pp_region_number
  count_grid0 <- as.matrix(table(index_grid)) # may have missing grid
  count_grid <- rep(0, n)
  count_grid[as.numeric(rownames(count_grid0))] <- count_grid0
  
  # run MCMC algorithm
  # CDMFM_new1 <- function(data, niterations, alpha, beta, GAMMA, LAMBDA, initNClusters,VN)
  temp_result <- CDMFM_new1(data = count_grid, niterations = niterations, alpha = a, beta = b, 
                            GAMMA = alpha, LAMBDA = 1, 
                            initNClusters = initNCluster, VN = VN)
  #
  ri_trace <- rep(NA, length(temp_result$Iterates))
  adj_ri_trace <- rep(NA, length(temp_result$Iterates))
  for(jj in 1:length(ri_trace)){
    ri_trace[jj] <- fossil::rand.index(temp_result$Iterates[[jj]]$zout, true_cluster_membership_setting)
    adj_ri_trace[jj] <- fossil::adj.rand.index(temp_result$Iterates[[jj]]$zout, true_cluster_membership_setting)
  }
  
  # dahl method
  temp_Dahl_result <- getDahl(temp_result, burn = burnin)
  
  # ri
  z_fitted <- temp_Dahl_result$zout
  ri <- fossil::rand.index(true_cluster_membership_setting, z_fitted)
  adj_ri <- fossil::adj.rand.index(true_cluster_membership_setting, z_fitted)
  # 
  # result <- temp_MFM_Dahl_rlt$zout
  result <- list(K_hat = max(z_fitted), ri = ri, adj_ri = adj_ri, samples = temp_result,
                 z_fitted = z_fitted, Dahl_result = temp_Dahl_result,
                 ri_trace = ri_trace, adj_ri_trace = adj_ri_trace, 
                 est_lambda = matrix(temp_Dahl_result$phiout[temp_Dahl_result$zout], 
                                     nrow = grid[1], ncol = grid[2]))
  result
}
#####################
time2 <- Sys.time()
# ######################
rlt_list <- rlt_setting4
# ######################
ri_trace_df_MFM <- data.frame(iter=NA, RI = NA, method = NA)
adj_ri_trace_df_MFM <- data.frame(iter=NA, RI = NA, method = NA)
#######################
for(ii in 1:nsim){
  ri_trace_df_MFM <- rbind(ri_trace_df_MFM,
                           data.frame(iter=1:length(rlt_list[[ii]]$ri_trace), 
                                      RI=rlt_list[[ii]]$ri_trace, 
                                      method=rep("MFM",length(rlt_list[[ii]]$ri_trace)))  ) 
  adj_ri_trace_df_MFM <- rbind(adj_ri_trace_df_MFM,
                               data.frame(iter=1:length(rlt_list[[ii]]$adj_ri_trace), 
                                          RI=rlt_list[[ii]]$adj_ri_trace, 
                                          method=rep("MFM",length(rlt_list[[ii]]$adj_ri_trace))))
  
}
########################
# remove NA's
ri_trace_df_MFM <- ri_trace_df_MFM[complete.cases(ri_trace_df_MFM),]
adj_ri_trace_df_MFM <- adj_ri_trace_df_MFM[complete.cases(adj_ri_trace_df_MFM),]
# ######################
# # save the results: Khat, RI, est_lambda
# ######################
est_lambda <- array(NA, dim=c(grid,nsim))
for(ii in 1:nsim){
  # 
  # figure out z_fitted
  # rlt_list[[ii]][[ minbic_model_id[ii] ]]$z_fitted
  est_lambda[,,ii] <- rlt_setting4[[ii]]$est_lambda
  
  # matrix(rlt_list[[ii]][[ minbic_model_id[ii] ]]$post_sample[[ rlt_list[[ii]][[minbic_model_id[ii]]]$Dahl_index ]]$lambda[rlt_list[[ii]][[ minbic_model_id[ii] ]]$z_fitted],
  #                              nrow=20,ncol=20)
  # dic_est_lambda[,,ii] <- matrix(rlt_list[[ii]][[ mindic_model_id[ii] ]]$post_sample[[ rlt_list[[ii]][[mindic_model_id[ii]]]$Dahl_index ]]$lambda[rlt_list[[ii]][[ mindic_model_id[ii] ]]$z_fitted],
  #                                nrow=20,ncol=20)
  # lpml_est_lambda[,,ii] <- matrix(rlt_list[[ii]][[ minlpml_model_id[ii] ]]$post_sample[[ rlt_list[[ii]][[minlpml_model_id[ii]]]$Dahl_index ]]$lambda[rlt_list[[ii]][[ minlpml_model_id[ii] ]]$z_fitted],
  #                                 nrow=20,ncol=20)
  
  # bic_ri[ii] <- rlt_list[[ii]][[minbic_model_id[ii]]]$ri
}
# ######################
simulation_results <- list(ri = do.call("c",lapply(rlt_list, function(x) x$ri)),
                           adj_ri = do.call("c",lapply(rlt_list, function(x) x$adj_ri)),
                           Khat = do.call("c",lapply(rlt_list, function(x) x$K_hat)),
                           ri_trace_df = ri_trace_df_MFM,
                           adj_ri_trace_df = adj_ri_trace_df_MFM,
                           est_lambda = est_lambda)
# ####################
######################
saveRDS(simulation_results, "MFM_simulation_results_setting4.rds")