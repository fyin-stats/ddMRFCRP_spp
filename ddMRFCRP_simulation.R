##################
########## generating data locally
###################
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
packages <- c("spatstat","doParallel","MHadaptive","purrr","spdep","expp","dplyr")
ipak(packages)
# library(spatstat)
# library(doParallel)
####################
# # uniform Poisson process with intensity 1 in a 10 x 10 square
# pp <- rpoispp(1, win=owin(c(0,10),c(0,10)))
# pp <- rpoispp(1, win=owin(c(0,10),c(0,10)))
source("./code/ddMRFCRP_gibbs.R")
####################
# setting 1
grid <- c(20, 20)
z <- matrix(0.2,nrow=grid[1],ncol=grid[2])
#
z[10:20,7:15] <- 12
z[1:5,1:18] <- 4
saveRDS(z, "true_lambda_setting1.rds")
# z[10:20,1:3] <- 8
# z[10:20,18:20] <- 8
stuff <- list(x=0.5:(grid[1]-0.5), y=0.5:(grid[2]-0.5), z=z)
Z <- as.im(stuff)
# true cluster membership
##########
# generate inhomogenous point patterns according to Z
nsim <- 100
pp_region_list_setting1 <- vector(mode="list",length=nsim)
for(ii in 1:nsim){
  set.seed(ii)
  pp <- rpoispp(Z)
  # generate region number
  pp_region_list_setting1[[ii]] <- pp2region(pp,nrow=grid[1],ncol=grid[2])
}
#####################
unique_cluster_intensity <- unique(c(z))
true_cluster_membership_setting1 <- rep(NA, prod(grid))
# 
for(i in 1:length(unique_cluster_intensity)){
    temp_arr_ind <- which(z == unique_cluster_intensity[i], arr.ind=T)
    temp_arr_ind <- data.frame(id=paste0(temp_arr_ind[,1],":",temp_arr_ind[,2]))
    # table join to figure out its number
    temp_arr_number <- (temp_arr_ind %>% inner_join(pp_region_list_setting1[[1]]$nbrt_regionid_number, by=c("id"="id")))$number
    true_cluster_membership_setting1[temp_arr_number] <- i
}
#####################
saveRDS(pp_region_list_setting1, "pp_region_list_setting1.rds")
saveRDS(true_cluster_membership_setting1, "true_cluster_membership_setting1.rds")
#####################
#####################
# setting 2
######################
grid <- c(20, 20)
z <- matrix(0.2,nrow=grid[1],ncol=grid[2])
#
z[1:8,11:20] <- 8
z[1:8,1:10] <- 1
z[9:16, 11:20] <- 4 
z[9:16, 1:10] <- 16
saveRDS(z, "true_lambda_setting2.rds")
###
stuff <- list(x=0.5:(grid[1]-0.5), y=0.5:(grid[2]-0.5), z=z)
Z <- as.im(stuff)
########################
nsim <- 100
pp_region_list_setting2 <- vector(mode="list",length=nsim)
for(ii in 1:nsim){
  set.seed(ii)
  pp <- rpoispp(Z)
  # generate region number
  pp_region_list_setting2[[ii]] <- pp2region(pp,nrow=grid[1],ncol=grid[2])
}
#####################
unique_cluster_intensity <- unique(c(z))
true_cluster_membership_setting2 <- rep(NA, prod(grid))
# 
for(i in 1:length(unique_cluster_intensity)){
  temp_arr_ind <- which(z == unique_cluster_intensity[i], arr.ind=T)
  temp_arr_ind <- data.frame(id=paste0(temp_arr_ind[,1],":",temp_arr_ind[,2]))
  # table join to figure out its number
  temp_arr_number <- (temp_arr_ind %>% inner_join(pp_region_list_setting2[[1]]$nbrt_regionid_number, by=c("id"="id")))$number
  true_cluster_membership_setting2[temp_arr_number] <- i
}
#####################
saveRDS(pp_region_list_setting2, "pp_region_list_setting2.rds")
saveRDS(true_cluster_membership_setting2, "true_cluster_membership_setting2.rds")
#####################
### setting 3
#####################
grid <- c(20, 20)
z <- matrix(0.2,nrow=grid[1],ncol=grid[2])
z[10:20,7:15] <- 20
z[1:5,1:18] <- 4
z[10:20,1:3] <- 10
z[10:20,18:20] <- 10
saveRDS(z, "true_lambda_setting3.rds")
#######################
stuff <- list(x=0.5:(grid[1]-0.5), y=0.5:(grid[2]-0.5), z=z)
Z <- as.im(stuff)
plot(Z)
##################
nsim <- 100
pp_region_list_setting3 <- vector(mode="list",length=nsim)
for(ii in 1:nsim){
  set.seed(ii)
  pp <- rpoispp(Z)
  # generate region number
  pp_region_list_setting3[[ii]] <- pp2region(pp,nrow=grid[1],ncol=grid[2])
}
#####################
unique_cluster_intensity <- unique(c(z))
true_cluster_membership_setting3 <- rep(NA, prod(grid))
# 
for(i in 1:length(unique_cluster_intensity)){
  temp_arr_ind <- which(z == unique_cluster_intensity[i], arr.ind=T)
  temp_arr_ind <- data.frame(id=paste0(temp_arr_ind[,1],":",temp_arr_ind[,2]))
  # table join to figure out its number
  temp_arr_number <- (temp_arr_ind %>% inner_join(pp_region_list_setting3[[1]]$nbrt_regionid_number, by=c("id"="id")))$number
  true_cluster_membership_setting3[temp_arr_number] <- i
}
######################
saveRDS(pp_region_list_setting3, "pp_region_list_setting3.rds")
saveRDS(true_cluster_membership_setting3, "true_cluster_membership_setting3.rds")
#######################

#######################
### setting 4: 6 clusters
#######################
### cluster 1: 3*3 (intensity: 16)
### cluster 2: 6*6 - 3*3 (intensity: 10)
### cluster 3: 9*9 - 6*6 - 3*3 (intensity: 6)
### cluster 4: 12*12 - 9*9 - 6*6 - 3*3 (intensity: 3)
### cluster 5: 16 * 20 - 12*12 - 9*9 - 6*6 - 3*3 (intensity: 1)
### cluster 6: 4 * 20 (intensity: 0.2)
######################
grid <- c(20, 20)
z <- matrix(NA,nrow=grid[1],ncol=grid[2])
z[1:4,1:20] <- 0.2
# 16*20 - ...
z[5:20,1:20] <- 4
# 12*15 - ...
z[8:20,3:18] <- 1
# 9*9 - ...
z[11:20,6:15] <- 20
# 5*5 - ...
z[15:20,8:13] <- 10
# 2*2 - ...
z[18:20,9:11] <- 40
# change the shape
z[11:15,6:15] <- 10
saveRDS(z, "true_lambda_setting4.rds")
#######################
stuff <- list(x=0.5:(grid[1]-0.5), y=0.5:(grid[2]-0.5), z=z)
Z <- as.im(stuff)
plot(Z)
#######################
nsim <- 100
pp_region_list_setting4 <- vector(mode="list",length=nsim)
for(ii in 1:nsim){
  set.seed(ii)
  pp <- rpoispp(Z)
  # generate region number
  pp_region_list_setting4[[ii]] <- pp2region(pp,nrow=grid[1],ncol=grid[2])
}
#####################
unique_cluster_intensity <- unique(c(z))
true_cluster_membership_setting4 <- rep(NA, prod(grid))
# 
for(i in 1:length(unique_cluster_intensity)){
  temp_arr_ind <- which(z == unique_cluster_intensity[i], arr.ind=T)
  temp_arr_ind <- data.frame(id=paste0(temp_arr_ind[,1],":",temp_arr_ind[,2]))
  # table join to figure out its number
  temp_arr_number <- (temp_arr_ind %>% inner_join(pp_region_list_setting4[[1]]$nbrt_regionid_number, by=c("id"="id")))$number
  true_cluster_membership_setting4[temp_arr_number] <- i
}
#####################
saveRDS(pp_region_list_setting4, "pp_region_list_setting4.rds")
saveRDS(true_cluster_membership_setting4, "true_cluster_membership_setting4.rds")
#####################