#### ddCRP with MRF
# library(MHadaptive)
# library(purrr)
# library(spdep)
# library(expp)
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("MHadaptive","purrr","dplyr","fossil")
ipak(packages)
# need package: spdep and expp
# https://cran.r-project.org/web/packages/spdep/vignettes/nb_igraph.html
##' model: lambda(s) = lambda0(s)exp(X * beta)
##' data: (x, y), length N 
##' grid: n

## variable selection, spike-slab for beta

##' function for Collapsed sampler;
##' input: X_grid (covariates value for every gird, n by p matrix);
##' input: n, number of grid;
##' input: index_grid, index showing which grid the data belongs to, N by 1 vector,
##'        integer valued from 1 to n;
##' input: hyperparameters, a, b, alpha;
##' input: A, area of every grid (equal area now);
##' input: distance, distance matrix, n by n matrix
##' input: decay function, f(d).
##' input: initNCluster, initial number of cluster
##' input: niterations, posterior sample size;
##' output: sample, matrix, with column: lambda0 (nCluster), index_cluster (n)
## ddCRP without covariates: ddMRFCRP
##' input: eta corresponds to the inverse temperature parameter in the Gibbs distribution
pp_Gibbs_ddMRFCRP <- function(n, index_grid, a, b, alpha, A, eta, neighbor_list,
                           initNCluster, niterations) {
  N <- length(index_grid)
  count_grid0 <- as.matrix(table(index_grid)) # may have missing grid
  count_grid <- rep(0, n)
  count_grid[as.numeric(rownames(count_grid0))] <- count_grid0
  
  ## initial posterior sample
  ## index_cluster: cluster index for each grid, corresponds to "z" in draft
  ## from 1:initNCluster, each cluster has elements (grids), no wasting number for cluster labeling.
  index_cluster <- c(sample(1:initNCluster, size = initNCluster, replace = FALSE),
                     sample(1:initNCluster, size = n-initNCluster, replace = TRUE))
  lambda <- rgamma(initNCluster, shape = a, rate = b)
  nCluster <- initNCluster
  
  
  History <- vector("list", niterations)
  
  ## start Gibb's sampling
  for (iter in 1:niterations) {
    
    ## permutate n grids
    np <- sample(1:n, size = n)
    ## random scan gibbs sampler
    ## update index_cluster
    for (i in np) {
      count_cluster <- table(index_cluster)
      ## determine if i-th grid is a singleton
      if (count_cluster[index_cluster[i]] > 1) {
        ## if not a singleton, then have nCluster + 1 choice
        count_cluster[index_cluster[i]] <- count_cluster[index_cluster[i]] - 1
        ## goes to an existing cluster: z_i = 1:nCluster
        logclusterProb <- sapply(1:nCluster, function(x) {
          # log(sum(f_d[i, -i] * (index_cluster[-i] == x)))
          # count the number of i's neighbors in this cluster
           eta * sum(index_cluster[neighbor_list[[i]]] == x) + 
            count_grid[i] * log(lambda[x]) - 
            A*lambda[x]
        })
        ## goes to a new cluster: z_i = nCluster+1 # this part remains unchanged
        logclusterProb[nCluster+1] <- log(alpha) + a * log(b) + 
          log(gamma(min(count_grid[i] + a, 150))) - 
          (count_grid[i] + a) * log(b+A) - 
          log(gamma(a))
        ## get the posterior sample for Z_i
        clusterProb <- exp(logclusterProb - max(logclusterProb))
        
        index_i <- sample(1:(nCluster+1), size = 1,
                          prob = clusterProb/sum(clusterProb))
        ## if the i-th grid really goes to a new cluster
        if (index_i > nCluster) {
          lambda_new <- rep(0, nCluster + 1)
          lambda_new[1:nCluster] <- lambda
          lambda_new[nCluster+1] <- rgamma(1, shape = a, rate = b)
          lambda <- lambda_new
          
          index_cluster[i] <- index_i
          nCluster <- nCluster + 1
        } else { ## if i-th grid goes to an existing cluster
          index_cluster[i] <- index_i
        }
      } else { ## if grid is a singleton, then has nCluster choices
        ## move all the cluster index that are greater than index_cluster[i] 
        ## foward by one to fill in the blank, also change the count_cluster 
        ## and lambda0 coresspondingly
        index_cluster[index_cluster > index_cluster[i]] <- 
          index_cluster[index_cluster > index_cluster[i]] - 1
        lambda[index_cluster[i] : (nCluster-1)] <- lambda[(index_cluster[i]+1):nCluster]
        count_cluster <- table(index_cluster)
        count_cluster[index_cluster[i]] <- count_cluster[index_cluster[i]] - 1
        ## only nCluster-1 clusters remained after removing i-th grid
        lambda <- lambda[1:(nCluster-1)]
        ## goes to existing cluster : 1:(nCluster-1)
        logclusterProb <- sapply(1:(nCluster-1), function(x) { 
          eta * sum(index_cluster[neighbor_list[[i]]] == x) + 
            count_grid[i] * log(lambda[x]) - 
            A*lambda[x]
        })
        ## goes to new cluster: nCluster
        logclusterProb[nCluster] <- log(alpha) + a * log(b) + 
          log(gamma(min(count_grid[i] + a, 150))) - 
          (count_grid[i] + a) * log(b+A) - 
          log(gamma(a))
        ## get the posterior sample for Z_i
        clusterProb <- exp(logclusterProb - max(logclusterProb))
        index_i <- sample(1:nCluster, size = 1, 
                          prob = clusterProb/sum(clusterProb))
        index_cluster[i] <- index_i
        if (index_i < nCluster) {
          nCluster <- nCluster-1
        } else {
          lambda_new <- rep(0, nCluster)
          lambda_new[1:(nCluster-1)] <- lambda
          lambda_new[nCluster] <- rgamma(1, shape = a, rate = b)
          lambda <- lambda_new
        }
      }
    }
    
    ## update lambda
    for (i in 1:nCluster) {
      lambda[i] <- rgamma(1, shape = a+sum(count_grid[index_cluster == i]), 
                          rate = b + A * sum(index_cluster == i))
    }
    
    
    
    History[[iter]] <- list(lambda = lambda, index_cluster = index_cluster)
    cat("eta =", eta, " iteration:", iter,"\n")
  }
  
  list(Iterates = History)
}

## exponential decay function
f_exp <- function(d, h, neighbor) {
  n <- nrow(d)
  d[d <= neighbor] <- 0
  x <- exp(-d*h)
  return(x/sum(x)*n*(n-1))
}

##' given grid, compute the n by n distance matrix, with diagonal elements being 10sqrt(2n)
dist_mat <- function(grid) {
  n <- prod(grid)
  distance <- matrix(10*n*sqrt(2), n, n)
  for (i in 1:(n-1)) {
    rowindex_o <- ceiling(i/grid[2])
    colindex_o <- i - (rowindex_o-1)*grid[2]
    for (j in (i+1):n) {
      rowindex <- ceiling(j/grid[2])
      colindex <- j - (rowindex-1)*grid[2]
      distance[i, j] <- sqrt((rowindex - rowindex_o)^2 + (colindex - colindex_o)^2)
      distance[j, i] <- distance[i, j]
    }
  }
  
  return(distance)
}

###################################################
RI <- function(z_true, samples) {
  n <- length(z_true)
  
  z0_pair <- combn(z_true, m = 2)
  
  ri <- rep(0, length(samples))
  for (i in 1:length(ri)) {
    z_pair <- combn(samples[[i]]$index_cluster, m = 2)
    a <- sum((z0_pair[1, ] == z0_pair[2, ]) * (z_pair[1, ] == z_pair[2, ]))
    b <- sum((z0_pair[1, ] != z0_pair[2, ]) * (z_pair[1, ] != z_pair[2, ]))
    ri[i] <- (a+b) / choose(n, 2)
  }
  return(ri)
}

#####################################################
dahl <- function(samples) {
  n <- length(samples[[1]]$index_cluster)
  niter <- length(samples)
  z_sample <- matrix(0, nrow = niter, ncol = n)
  for (i in 1:niter) {
    z_sample[i, ] <- samples[[i]]$index_cluster
  }
  BList <- purrr::map(1:niter, ~outer(z_sample[.x,], z_sample[.x,], "=="))
  BBar <- Reduce("+", BList) / niter
  SSE <- map_dbl(BList, ~sum((.x - BBar)^2))
  return(list(min.sse = min(SSE), cLS = which.min(SSE), 
              cluster = as.numeric(z_sample[which.min(SSE),])))
}


###############################################################
## sample after burnin
LPML <- function(samples, index_grid, A = 1) {
  lambda1 <- rep(0, length = length(index_grid))
  lambda2 <- rep(0, length = length(samples[[1]]$index_cluster))
  
  for (j in 1:length(index_grid)) {
    lambda_sample <- rep(0, length = length(samples))
    for (b in 1:length(samples)) {
      lambda_sample[b] <- samples[[b]]$lambda[samples[[b]]$index_cluster[index_grid[j]]]
    }
    lambda1[j] <- 1/mean(1/lambda_sample)
  }
  
  for (i in 1:length(samples[[1]]$index_cluster)) {
    lambda_sample <- rep(0, length = length(samples))
    for (b in 1:length(samples)) {
      lambda_sample[b] <- samples[[b]]$lambda[samples[[b]]$index_cluster[i]]
    }
    lambda2[i] <- mean(lambda_sample)
  }
  
  LPML <- sum(log(lambda1)) - A * sum(lambda2)
  return(LPML)
}

#################################################################
## samples after burnin
DIC <- function(samples, index_grid, grid = c(20, 20)) {
  
  n <- length(samples[[1]]$index_cluster)
  
  ## for each poster sample draw
  lambda_grid <- matrix(0, nrow = n, ncol = length(samples))
  dev_sample <- rep(0, length(samples))
  for (i in 1:length(samples)) {
    for (j in 1:n) {
      lambda_grid[j, i] <- samples[[i]]$lambda[samples[[i]]$index_cluster[j]]
    }
    dev_sample[i] <- -2 * (sum(log(lambda_grid[index_grid, i])) - sum(lambda_grid[, i]))
  }
  
  ## first get posterior mean of parameters
  beta_sample <- matrix(0, nrow = length(samples), ncol = n)
  lambda_sample <- matrix(0, nrow = length(samples), ncol = n)
  for (i in 1:length(samples)) {
    lambda_sample[i, ] <- samples[[i]]$lambda[samples[[i]]$index_cluster]
  }
  lambda_pos <- colMeans(lambda_sample)
  dev_mean <- -2 * (sum(log(lambda_pos[index_grid])) - sum(lambda_pos))
  
  pd <- mean(dev_sample) - dev_mean
  return(list(DIC = dev_mean + 2*pd, pD = pd))
}

BIC <- function(samples, Dahl_index, index_grid, A = 1) {
  N <- length(index_grid)
  
  
  lambda_data <- samples[[Dahl_index]]$lambda[samples[[Dahl_index]]$index_cluster[index_grid]]
  lambda_area <- samples[[Dahl_index]]$lambda[samples[[Dahl_index]]$index_cluster]
  
  loglikelihood <- sum(log(lambda_data)) - A * sum(lambda_area)
  
  n_par <- length(samples[[Dahl_index]]$lambda)
  
  return(-2 * loglikelihood + n_par * log(N))
}
###############################
### function: pp2regionid
###############################
pp2region <- function(pp,nrow,ncol,neigh_type="rook",...){
  ## input: pp with x and y
  ## output: pp region id:, e.g., "1:1"; and region number; and neighorlist
  # region number
  nbrt <- cell2nb(nrow=nrow, ncol = ncol, type=neigh_type, legacy = TRUE)
  # neighbordataframe
  nbrt_df <- neighborsDataFrame(nbrt)
  # need to map the id "a:b" to a number
  nbrt_regionid_number <- data.frame(id=attr(nbrt, "region.id"),
                                      number=1:length(attr(nbrt, "region.id")))
  #
  pp_region_id <- c()
  pp_region_number <- c()
  for(i in 1:length(pp$x)){
    x_id <- ceiling(pp$x[i])
    y_id <- ceiling(pp$y[i])
    pp_region_id <- c(pp_region_id, paste0(x_id,":",y_id))
    pp_region_number <- c(pp_region_number, nbrt_regionid_number$number[which(nbrt_regionid_number$id==pp_region_id[i])])
  }
  # first figure out the numbers for id
  # second figure out the numbers for id_neigh
  nbrt_df_number <- nbrt_df %>% inner_join(nbrt_regionid_number, by = c("id"="id")) %>% inner_join(nbrt_regionid_number, by = c("id_neigh"="id"))
  #
  colnames(nbrt_df_number) <- c("id", "id_neigh", "number", "number_neigh")
  # now generate a list, for each grid (represented by a number), 
  # figure out the numbers of its neighbors
  neighbor_list <- vector(mode="list", length = length(unique(nbrt_df_number$number)))
  ####
  for(i in 1:length(neighbor_list)){
    neighbor_list[[i]] <- filter(nbrt_df_number, number==i)$number_neigh
  }
  
  #
  return(list(pp_region_id = pp_region_id,pp_region_number = pp_region_number, neighbor_list = neighbor_list,
              nbrt_regionid_number = nbrt_regionid_number))
}
## test
# samples <- pp_Gibbs_ddCRP_cov(X_grid = simudata$X_grid, n = n, index_grid = index_grid, 
#                           a = a, b = b, alpha = alpha, A = A, distance = distance, 
#                           f_decay, initNCluster = initNCluster, niterations = niterations)

## focus on 35*50 grid
# setwd("/Users/fan/Documents/dd_spp_grf")
# curry_data <- read.csv("./data/curry.csv")
# curry_data <- subset(curry_data, y >= -47.5 & y <= 302.5)
# curry_data$x <- (curry_data$x+250)/10 # rescale the x coordinate value
# curry_data$y <- (curry_data$y+47.5)/10 # rescale the y coordinate value
# grid <- c(35, 50)
# n <- prod(grid)
# 
# colindex <- ceiling(curry_data$x)
# rowindex <- ceiling(curry_data$y)
# index_grid <- (rowindex - 1) * grid[2] + colindex # data
# a <- 1
# b <- 1
# alpha <- 1
# initNCluster <- 5
# niterations <- 50000
# burnin <- 20000
# thin <- 10
# 
# h <- seq(0, 0.1, by = 0.01)
# neighbor <- 5
# distance <- dist_mat(grid) # 1750 by 1750 matrix, because we have 1750 grids
# 
# results <- vector("list", length = length(h))
# for (i in 1:length(h)) {
#   fd <- f_exp(distance, h[i], neighbor)
#   sample <- pp_Gibbs_ddCRP(n = n, index_grid = index_grid,
#                            a = a, b = b, alpha = alpha, A = 1,
#                            f_d = fd, initNCluster = initNCluster,
#                            niterations = niterations, h = h[i])
#   
#   post_sample <- sample$Iterates[seq(from = (burnin+1), to = niterations, by = thin)]
#   Dahl_index <- dahl(post_sample)$cLS
#   z_fitted <- post_sample[[Dahl_index]]$index_cluster
#   ri <- RI(z_fitted, sample$Iterates)
#   bic <- BIC(samples = post_sample,
#              Dahl_index = Dahl_index, index_grid = index_grid, A = 1)
#   
#   lpml <- LPML(samples = post_sample,
#                index_grid = index_grid, A = 1)
#   dic <- DIC(samples = post_sample,
#              index_grid = index_grid, grid = grid)$DIC
#   results[[i]] <- list(sample = sample, simudata = simudata, ri = ri, Dahl_index = Dahl_index,
#                        bic = bic, lpml = lpml, dic = dic)
#   
# }
# 
# save(results, file = "~/Desktop/results_neighbor5.rdata")
##################################################################