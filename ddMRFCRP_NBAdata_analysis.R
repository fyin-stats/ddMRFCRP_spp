########################
############ real data analysis
########################
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
########################################
# focus on 50*35 grid
# setwd("/Users/fan/Documents/dd_spp_grf")
# curry_data <- read.csv("./data/curry.csv")
# curry_data <- subset(curry_data, y >= -47.5 & y <= 302.5)
# curry_data$x <- (curry_data$x+250)/10 # rescale the x coordinate value
# curry_data$y <- (curry_data$y+47.5)/10 # rescale the y coordinate value
# grid <- c(35, 50)
n <- prod(grid)
# #colindex <- ceiling(curry_data$x)
# #rowindex <- ceiling(curry_data$y)
# #index_grid <- (rowindex - 1) * grid[2] + colindex # data
# pp_region_curry <- pp2region(curry_data, nrow=35, ncol=50)
# saveRDS(pp_region_curry, "pp_region_curry.rds")
################################
# read the data
# pp_region <- readRDS("./data/pp_region_curry.rds")
##
a <- 1
b <- 1
alpha <- 1
initNCluster <- 5
niterations <- 4000
burnin <- 2000
thin <- 10
eta_vec <- seq(0,7,by=0.5)
# h <- seq(0, 0.1, by = 0.01)
# neighbor <- 5
# distance <- dist_mat(grid) # 1750 by 1750 matrix, because we have 1750 grids
# results <- vector("list", length = length(eta_vec))
# run the algorithm under different values of inverse temperature (eta)
registerDoParallel(cores = length(eta_vec))
#
time3 <- Sys.time()
results <- foreach(k = 1:length(eta_vec)) %dopar% {
  sample <- pp_Gibbs_ddMRFCRP(n = prod(grid), index_grid = pp_region$pp_region_number,
                              a = a, b = b, alpha = alpha, A = 1, eta = eta_vec[k],
                              neighbor_list = pp_region$neighbor_list,
                              initNCluster = initNCluster,
                              niterations = niterations)
  # n, index_grid, a, b, eta, alpha, A, neighbor_list,
  # initNCluster, niterations
  post_sample <- sample$Iterates[seq(from = (burnin+1), to = niterations, by = thin)]
  # ###
  # ri_trace <- rep(NA, length(sample$Iterates)-burnin)
  # adj_ri_trace <- rep(NA, length(sample$Iterates)-burnin)
  # for(jj in 1:length(ri_trace)){
  #   ri_trace[jj] <- fossil::rand.index(true_cluster_membership_setting, sample$Iterates[[burnin + jj]]$index_cluster)
  #   adj_ri_trace[jj] <- fossil::adj.rand.index(true_cluster_membership_setting, sample$Iterates[[burnin + jj]]$index_cluster)
  # }
  ###
  Dahl_index <- dahl(post_sample)$cLS
  z_fitted <- post_sample[[Dahl_index]]$index_cluster
  ri <- RI(z_fitted, sample$Iterates)
  bic <- BIC(samples = post_sample,
             Dahl_index = Dahl_index, index_grid = pp_region$pp_region_number, A = 1)
  
  lpml <- LPML(samples = post_sample,
               index_grid = pp_region$pp_region_number, A = 1)
  dic <- DIC(samples = post_sample,
             index_grid = pp_region$pp_region_number, grid = grid)$DIC
  
  temp_result <- list(ri = ri, Dahl_index = Dahl_index, Khat = max(z_fitted),
                      z_fitted=z_fitted, post_sample = post_sample,
                      bic = bic, lpml = -lpml, dic = dic)
  temp_result
}
###########
# 
# rlt_list <- results
dic_vec <- do.call("c", lapply(results, function(x) x$dic))
bic_vec <- do.call("c", lapply(results, function(x) x$bic))
lpml_vec <- do.call("c", lapply(results, function(x) x$lpml))
#
minbic_model_id <- which(bic_vec == min(bic_vec))
mindic_model_id <- which(dic_vec == min(dic_vec))
minlpml_model_id <- which(lpml_vec == min(lpml_vec))
# figure out the Khat estimate and the lambda estimate and the ri
results_summary <- list(bic_Khat = results[[minbic_model_id]]$Khat,
                        dic_Khat = results[[mindic_model_id]]$Khat,
                        lpml_Khat = results[[minlpml_model_id]]$Khat,
                        bic_est_lambda = matrix(results[[ minbic_model_id ]]$post_sample[[ results[[minbic_model_id]]$Dahl_index ]]$lambda[results[[ minbic_model_id ]]$z_fitted],
                                                nrow=grid[1],ncol=grid[2]),
                        dic_est_lambda = matrix(results[[ mindic_model_id ]]$post_sample[[ results[[mindic_model_id]]$Dahl_index ]]$lambda[results[[ mindic_model_id ]]$z_fitted],
                                                nrow=grid[1],ncol=grid[2]),
                        lpml_est_lambda = matrix(results[[ minlpml_model_id ]]$post_sample[[ results[[minlpml_model_id]]$Dahl_index ]]$lambda[results[[ minlpml_model_id ]]$z_fitted],
                                                nrow=grid[1],ncol=grid[2]),
                        bic_ri =  results[[minbic_model_id]]$ri,
                        dic_ri =  results[[mindic_model_id]]$ri,
                        lpml_ri =  results[[minlpml_model_id]]$ri, 
                        minbic_model_id = minbic_model_id,
                        mindic_model_id = mindic_model_id,
                        minlpml_model_id = minlpml_model_id)
time4 <- Sys.time()
###################
sink("./results_summary.txt", append = TRUE)
print(data_filenames[i])
print(results_summary)
print(time4-time3)
cat("\n")
sink()
###################
# dic_matrix <- matrix(NA, nrow=nsim, ncol = length(eta_vec))
# lpml_matrix <- matrix(NA, nrow=nsim, ncol=length(eta_vec))
# bic_matrix <- matrix(NA, nrow=nsim, ncol=length(eta_vec))
# for(ii in 1:nsim){
#   dic_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$dic))
#   lpml_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$lpml))
#   bic_matrix[ii,] <- do.call("c",lapply(rlt_list[[ii]], function(x) x$bic))
# }
# ######################
# minbic_model_id <- apply(bic_matrix,1,function(x) which(x==min(x)))
# mindic_model_id <- apply(dic_matrix,1,function(x) which(x==min(x)))
# minlpml_model_id <- apply(lpml_matrix,1,function(x) which(x==min(x)))
# saveRDS(results,"./results_neighbor5_neweta.rds")