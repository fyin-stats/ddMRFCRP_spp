########################
########################
### main function for NBA data analysis
########################
########################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("spatstat","doParallel","MHadaptive","purrr","dplyr","stringi")
ipak(packages)
# library(spatstat)
# library(doParallel)
####################
# # uniform Poisson process with intensity 1 in a 10 x 10 square
# pp <- rpoispp(1, win=owin(c(0,10),c(0,10)))
# pp <- rpoispp(1, win=owin(c(0,10),c(0,10)))
source("./code/ddMRFCRP_gibbs.R")
# source("./code/ddMRFCRP_gibbs.R")
source("./code/MFM_example.R")
#####################
# data_filenames <- list.files("./data/")
# data_filenames <- c("pp_region_curry.rds", "pp_region_durant.rds", "pp_region_james.rds", 
#                     "pp_region_jordan.rds", "pp_region_harden.rds", "pp_region_Embiid.rds")
data_filenames <- c("pp_region_curry.rds",
                    "pp_region_jordan.rds",
                    "pp_region_durant.rds",
                    "pp_region_james.rds",
                    "pp_region_harden.rds",
                    "pp_region_Embiid.rds",
                    "pp_region_Antetokounmpo.rds",
                    "pp_region_DeRozan.rds",
                    "pp_region_george.rds",
                    "pp_region_Lillard.rds", 
                    "pp_region_griffin.rds",
                    "pp_region_towns.rds", 
                    "pp_region_westbrook.rds", 
                    "pp_region_howard.rds", 
                    "pp_region_paul.rds", 
                    "pp_region_Aldridge.rds",
                    "pp_region_Porzingis.rds", 
                    "pp_region_irving.rds",
                    "pp_region_Butler.rds", 
                    "pp_region_thompson.rds")
# data_filenames <- c("pp_region_curry.rds", "pp_region_jordan.rds")
# results_summary_list <- vector(mode="list", length = length(data_filenames))
# results_list <- vector(mode="list", length = length(data_filenames))
# 
grid <- c(50,35)
n <- prod(grid)
a <- 1
b <- 1
# dirichlet prior
alpha <- 1
# initial number of clusters
initNCluster <- 5
niterations <- 4000
burnin <- 2000
thin <- 10
#########################
registerDoParallel(cores = length(data_filenames))
#########################
time1 <- Sys.time()
results_list <- foreach(i = 1:length(data_filenames)) %dopar% {
  pp_region <- readRDS(paste0("./rdsdata/",data_filenames[i]) )
  # source("./ddMRFCRP_NBAdata_analysis.R")
  # source("./ppMFM_NBAdata_analysis.R")
  # results_summary_list[[i]] <- result
  # results_list[[i]] <- result
  index_grid <- pp_region$pp_region_number
  count_grid0 <- as.matrix(table(index_grid)) # may have missing grid
  count_grid <- rep(0, n)
  count_grid[as.numeric(rownames(count_grid0))] <- count_grid0
  # run MCMC algorithm
  # CDMFM_new1 <- function(data, niterations, alpha, beta, GAMMA, LAMBDA, initNClusters,VN)
  temp_result <- CDMFM_new1(data = count_grid, 
                            niterations = niterations, 
                            alpha = a, beta = b, 
                            GAMMA = alpha, LAMBDA = 1, 
                            initNClusters = initNCluster, 
                            VN = VN)
  #
  # ri_trace <- rep(NA, length(temp_result$Iterates))
  # adj_ri_trace <- rep(NA, length(temp_result$Iterates))
  # for(jj in 1:length(ri_trace)){
  #   ri_trace[jj] <- fossil::rand.index(temp_result$Iterates[[jj]]$zout, true_cluster_membership_setting)
  #   adj_ri_trace[jj] <- fossil::adj.rand.index(temp_result$Iterates[[jj]]$zout, true_cluster_membership_setting)
  # }
  # dahl method
  temp_Dahl_result <- getDahl(temp_result, burn = burnin)
  # ri
  z_fitted <- temp_Dahl_result$zout
  # ri <- fossil::rand.index(true_cluster_membership_setting, z_fitted)
  # adj_ri <- fossil::adj.rand.index(true_cluster_membership_setting, z_fitted)
  # result <- temp_MFM_Dahl_rlt$zout
  # ri = ri, 
  # adj_ri = adj_ri, 
  # ri_trace = ri_trace, adj_ri_trace = adj_ri_trace
  result <- list(K_hat = max(z_fitted), 
                 samples = temp_result,
                 z_fitted = z_fitted, Dahl_result = temp_Dahl_result,
                 est_lambda = matrix(temp_Dahl_result$phiout[temp_Dahl_result$zout], 
                                     nrow = grid[1], ncol = grid[2]))
  result
}
time2 <- Sys.time()
##########################
# saveRDS(results_summary_list, "ppMFM_NBAdata_results_summary_list_new.rds")
saveRDS(results_list, "ppMFM_NBAdata_results_list_new.rds")
# ppMFM_NBAdata_results_list_new.rds
##########################