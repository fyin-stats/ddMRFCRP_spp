########################
########### summary of NBA data analysis results
#########################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("dplyr","fossil","xtable","ggplot2","fields","lattice","gridExtra")
ipak(packages)
######################
######################
######################
######################
grid <- c(50, 35)
tool_matrix <- matrix(1, nrow = grid[1], grid[2])
#
X <- data.frame(cbind(which(tool_matrix==1, arr.ind = T),1:prod(grid)))
colnames(X) <- c("x", "y", "linear_index")

############# MFM-NHPP for NBA players
NBAdata_results_summary_list_all <- readRDS("./ppMFM_NBAdata_results_list_new.rds") 
MFM_NBAdata_results_summary_list_all <- NBAdata_results_summary_list_all
######################
player_names <- c("curry", "jordan", "durant", "james", "harden", "Embiid", "Antetokounmpo",
                  "DeRozan", "george", "Lillard", "griffin", "towns", "westbrook", "howard", 
                  "paul", "Aldridge","Porzingis", "irving", "Butler", "thompson")
player_positions <- c("PG", "C", "SF", "SF", "SG", "C", "PF", "SG", "SF", "PG",
                      "PF", "C", "SG", "C", "PG", "PF", "PF", "PG", "SF", "SG")
################
# p_list <- vector(mode="list", length = length(NBAdata_results_summary_list_all))
# ################
# for(i in 1:length(NBAdata_results_summary_list_all)){
#   inten0fit_dat <- data.frame(y = X$y, x = X$x, value = c(NBAdata_results_summary_list_all[[i]]$est_lambda), 
#                               player = player_names[i])
#   #############
#   breaks <- seq(0, max(c(inten0fit_dat$value))+1, length.out = 100)
#   cols <- terrain.colors(100)[length(terrain.colors(100)):1]
#   p_list[[i]] <- levelplot(value~x*y,inten0fit_dat, 
#                            col.regions = cols, at = breaks, 
#                            xlab = "X", ylab = "Y", main = player_names[i],
#                            scales = list(x = list(at = seq(0, 50, length.out = 11)), 
#                                          y = list(at = seq(0, 35, length.out = 8))))
# }
#
# results for MRFCRP
player_est_lambda_list <- vector(mode="list", length = length(NBAdata_results_summary_list_all))
# 
for(i in 1:length(NBAdata_results_summary_list_all)){
  player_est_lambda_list[[i]] <- list(est_lambda = NBAdata_results_summary_list_all[[i]]$est_lambda,
                                      player_name = player_names[i],
                                      player_position = player_positions[i])
  # player_est_lambda_list[[i]]$player_name <- player_names[i]
}
saveRDS(player_est_lambda_list, "MFM_player_est_lambda_list.rds")
################### MRF-DPM-NHPP, NBA data analysis
# setwd("/Users/fan/Documents/dd_spp_grf/")
NBAdata_results_summary_list <- readRDS("./NBAdata_results_summary_list.rds")
NBAdata_results_summary_list_new <- readRDS("./NBAdata_results_summary_list_new.rds")
NBAdata_results_summary_list_new2 <- readRDS("./NBAdata_results_summary_list_new2.rds")
NBAdata_results_summary_list_new3 <- readRDS("./NBAdata_results_summary_list_new3.rds")
NBAdata_results_summary_list_new4 <- readRDS("./NBAdata_results_summary_list_new4.rds")
NBAdata_results_summary_list_new5 <- readRDS("./NBAdata_results_summary_list_new5.rds")
##########
NBAdata_results_summary_list_all <- c(NBAdata_results_summary_list,
                                      NBAdata_results_summary_list_new,
                                      NBAdata_results_summary_list_new2,
                                      NBAdata_results_summary_list_new3,
                                      NBAdata_results_summary_list_new4,
                                      NBAdata_results_summary_list_new5)
######################
grid <- c(50, 35)
player_names <- c("curry", "jordan", "durant", "james", "harden", "Embiid", "Antetokounmpo",
                  "DeRozan", "george", "Lillard", "griffin", "towns", "westbrook", "howard", 
                  "paul", "Aldridge","Porzingis", "irving", "Butler", "thompson")
################
p_list <- vector(mode="list", length = length(NBAdata_results_summary_list_all))
################
for(i in 1:length(NBAdata_results_summary_list_all)){
  inten0fit_dat <- data.frame(y = X$y, x = X$x, value = c(NBAdata_results_summary_list_all[[i]]$bic_est_lambda), 
                              player = player_names[i])
  #############
  breaks <- seq(0, max(c(inten0fit_dat$value))+1, length.out = 100)
  cols <- terrain.colors(100)[length(terrain.colors(100)):1]
  p_list[[i]] <- levelplot(value~x*y,inten0fit_dat, 
                         col.regions = cols, at = breaks, 
                         xlab = "X", ylab = "Y", main = player_names[i],
                         scales = list(x = list(at = seq(0, 50, length.out = 11)), 
                                       y = list(at = seq(0, 35, length.out = 8))))
}
#
# results for MRFCRP
player_est_lambda_list <- vector(mode="list", length = length(NBAdata_results_summary_list_all))
# 
for(i in 1:length(NBAdata_results_summary_list_all)){
  player_est_lambda_list[[i]] <- NBAdata_results_summary_list_all[[i]]
  player_est_lambda_list[[i]]$player_name <- player_names[i]
  player_est_lambda_list[[i]]$player_position <- player_positions[i]
}
saveRDS(player_est_lambda_list, "MRF_player_est_lambda_list.rds")

########## Table for players
# 
eta_vec <- seq(0,7,by=0.5)
player_names <- c("curry", "jordan", "durant", "james", "harden", "Embiid", "Antetokounmpo",
                  "DeRozan", "george", "Lillard", "griffin", "towns", "westbrook", "howard", 
                  "paul", "Aldridge","Porzingis", "irving", "Butler", "thompson")
player_full_names <- c("Stephen Curry", "DeAndre Jordan", "Kevin Durant", "LeBron James", "James Harden",
                       "Joel Embiid", "Giannis Antetokounmpo", "DeMar DeRozan", "Paul George", "Damian Lillard",
                       "Blake Griffin", "Karl-Anthony Towns", "Russell Westbrook", "Dwight Howard", "Chris Paul",
                       "LaMarcus Aldridge", "Kristaps Porziņģis", "Kyrie Irving", "Jimmy Butler", "Klay Thompson")
player_table <- data.frame(player_names = player_full_names,
                           position = player_positions,
                           bic_Khat = do.call("c", lapply(NBAdata_results_summary_list_all, function(x) x$bic_Khat)),
                           bic_eta = eta_vec[do.call("c", lapply(NBAdata_results_summary_list_all, function(x) x$minbic_model_id))],
                           dic_Khat = do.call("c", lapply(NBAdata_results_summary_list_all, function(x) x$dic_Khat)),
                           dic_eta = eta_vec[do.call("c", lapply(NBAdata_results_summary_list_all, function(x) x$mindic_model_id))],
                           lpml_Khat = do.call("c", lapply(NBAdata_results_summary_list_all, function(x) x$lpml_Khat)),
                           lpml_eta = eta_vec[do.call("c", lapply(NBAdata_results_summary_list_all, function(x) x$minlpml_model_id))])
# reorder the table by column position
#
xtable(player_table %>% dplyr::arrange(position), digits = 1)

# BIC only + MFM
player_table <- data.frame(player_names = player_full_names,
                           position = player_positions,
                           bic_Khat = do.call("c", lapply(NBAdata_results_summary_list_all, function(x) x$bic_Khat)),
                           bic_eta = eta_vec[do.call("c", lapply(NBAdata_results_summary_list_all, function(x) x$minbic_model_id))],
                           MFM_Khat = do.call("c", lapply(MFM_NBAdata_results_summary_list_all, function(x) x$K_hat)))
# reorder the table by column position
#
xtable(player_table %>% dplyr::arrange(position), digits = 1)