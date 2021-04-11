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
#####################
# data_filenames <- list.files("./data/")
# data_filenames <- c("pp_region_curry.rds", "pp_region_durant.rds", "pp_region_james.rds", 
#                     "pp_region_jordan.rds", "pp_region_harden.rds", "pp_region_Embiid.rds")
data_filenames <- c("pp_region_Porzingis.rds", "pp_region_irving.rds")
# data_filenames <- c("pp_region_curry.rds", "pp_region_jordan.rds")
results_summary_list <- vector(mode="list", length = length(data_filenames))
results_list <- vector(mode="list", length = length(data_filenames))
# 
grid <- c(50,35)
#########################
time1 <- Sys.time()
for(i in 1:length(data_filenames)){
  pp_region <- readRDS(paste0("./rdsdata/",data_filenames[i]) )
  source("./ddMRFCRP_NBAdata_analysis.R")
  results_summary_list[[i]] <- results_summary 
  results_list[[i]] <- results
}
time2 <- Sys.time()
##########################
saveRDS(results_summary_list, "NBAdata_results_summary_list_new4.rds")
saveRDS(results_list, "NBAdata_results_list_new4.rds")