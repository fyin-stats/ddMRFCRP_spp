################# create neighors for grid 
library(MHadaptive)
library(purrr)
library(spdep)
library(expp)
library(dplyr)
#################
#################

# create nb object using cell2nb function
nb7rt <- cell2nb(nrow=7, ncol = 7, type="rook")
# generate neighbor dataframe
nb7rt_df <- neighborsDataFrame(nb7rt)
# need to map the id "a:b" to a number
nb7rt_regionid_number <- data.frame(id=attr(nb7rt, "region.id"),
                                    number=1:length(attr(nb7rt, "region.id")))
##################
# first figure out the numbers for id
# second figure out the numbers for id_neigh
nb7rt_df_number <- nb7rt_df %>% inner_join(nb7rt_regionid_number, by = c("id"="id")) %>% inner_join(nb7rt_regionid_number, by = c("id_neigh"="id"))
#
colnames(nb7rt_df_number) <- c("id", "id_neigh", "number", "number_neigh")
##################
# now generate a list, for each grid (represented by a number), 
# figure out the numbers of its neighbors
neighbor_list <- vector(mode="list", length = length(unique(nb7rt_df_number$number)))
####
for(i in 1:length(neighbor_list)){
  neighbor_list[[i]] <- filter(nb7rt_df_number, number==i)$number_neigh
}

















