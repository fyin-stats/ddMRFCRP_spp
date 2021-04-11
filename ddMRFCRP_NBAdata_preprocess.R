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
data_filenames <- list.files("./data/")
# data_filenames <- c("curry.csv", "Antetokounmpo.csv", "james.csv", "jordan.csv", "harden.csv")
# results_summary_list <- vector(mode="list", length = length(data_filenames))
#########################
for(i in 1:length(data_filenames)){
  temp_data <- read.csv(paste0("./data/",data_filenames[i]))
  # read.csv("./data/curry.csv")
  temp_data <- subset(temp_data, y >= -47.5 & y <= 302.5)
  temp_data$x <- (temp_data$x+250)/10 # rescale the x coordinate value
  temp_data$y <- (temp_data$y+47.5)/10 # rescale the y coordinate value
  grid <- c(50, 35)
  n <- prod(grid)
  # #colindex <- ceiling(curry_data$x)
  # #rowindex <- ceiling(curry_data$y)
  # #index_grid <- (rowindex - 1) * grid[2] + colindex # data
  pp_region <- pp2region(temp_data, nrow=grid[1], ncol=grid[2])
  # 
  # plot the location of these points
  # temp_x <- rep(NA, length(pp_region$pp_region_id))
  # temp_y <- rep(NA, length(pp_region$pp_region_id))
  # for(k in 1:length(pp_region$pp_region_id)){
  #   temp_x[k] <- as.numeric(stri_sub(pp_region$pp_region_id[k], from = 1L, to = which(strsplit(pp_region$pp_region_id[k], "")[[1]]==":")-1))
  #   temp_y[k] <- as.numeric(stri_sub(pp_region$pp_region_id[k], from = which(strsplit(pp_region$pp_region_id[k], "")[[1]]==":")+1, to = -1L))
  # }
  
  #
  saveRDS(pp_region, paste0("./rdsdata/pp_region_", stri_sub(data_filenames[i], from = 1L, to = which(strsplit(data_filenames[i], "")[[1]]==".")) ,"rds") )
  # source("./ddMRFCRP_NBAdata_analysis.R")
  # results_summary_list[[i]] <- results_summary 
}
# ########################
# for(i in 1:length(data_filenames)){
#   temp_data <- read.csv(paste0("./data/",data_filenames[i]))
#   # read.csv("./data/curry.csv")
#   temp_data <- subset(temp_data, y >= -47.5 & y <= 352.5)
#   temp_data$x <- (temp_data$x+250)/10 # rescale the x coordinate value
#   temp_data$y <- (temp_data$y+47.5)/10 # rescale the y coordinate value
#   grid <- c(50, 40)
#   n <- prod(grid)
#   # #colindex <- ceiling(curry_data$x)
#   # #rowindex <- ceiling(curry_data$y)
#   # #index_grid <- (rowindex - 1) * grid[2] + colindex # data
#   pp_region <- pp2region(temp_data, nrow=grid[1], ncol=grid[2])
#   #
#   saveRDS(pp_region, paste0("./rdsdata/pp_region_50_47_", stri_sub(data_filenames[i], from = 1L, to = which(strsplit(data_filenames[i], "")[[1]]==".")) ,"rds") )
#   # source("./ddMRFCRP_NBAdata_analysis.R")
#   # results_summary_list[[i]] <- results_summary
# }