ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("dplyr","fossil","xtable","ggplot2","fields","lattice","gridExtra","jpeg","RCurl")
ipak(packages)
##########################
##########################
courtImg.URL <- "https://thedatagame.files.wordpress.com/2016/03/nba_court.jpg"

##########################
court <- rasterGrob(readJPEG(getURLContent(courtImg.URL)),
                    width=unit(1,"npc"), height=unit(1,"npc"))
########################
setwd("/Users/fan/Documents/dd_spp_grf")
# data_filenames <- list.files("./data/")
# data_filenames <- c("curry.csv", "Antetokounmpo.csv", "james.csv", "jordan.csv", "harden.csv")
# results_summary_list <- vector(mode="list", length = length(data_filenames))
data_filenames <- c("curry.csv","durant.csv","harden.csv","DeRozan.csv" )
# DeRozan.csv
#
data_list <- vector(mode="list", length=length(data_filenames))
#########################
for(i in 1:length(data_list)){
  temp_data <- read.csv(paste0("./data/",data_filenames[i]))
  # read.csv("./data/curry.csv")
  temp_data <- subset(temp_data, y >= -47.5 & y <= 302.5)
  # temp_data$x <- (temp_data$x+250)/10 # rescale the x coordinate value
  # temp_data$y <- (temp_data$y+47.5)/10 # rescale the y coordinate value
  data_list[[i]] <- temp_data
  # grid <- c(50, 35)
  # n <- prod(grid)
  # pp_region <- pp2region(temp_data, nrow=grid[1], ncol=grid[2])
}
########################
# mark_curry <- as.factor(shotdata[[5]]$shot_made_flag)

# levels(mark_curry) <- c("missed", "made")

p1 <- ggplot(data_list[[1]], aes(x=x, y=y)) + coord_fixed(ratio = 47/50) +
  
  annotation_custom(court, -250, 250, -50, 420) +
  
  geom_point(alpha = 0.5, size = 1, color="blue") + xlim(-250, 250) +
  
  ylim(-50, 420) +
  
  annotate("text", x = 160, y = 400, label="Curry", size=8)+
  
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),
        
        panel.grid.major = element_blank(), legend.justification = c(0, 1),
        
        legend.position = c(0.05, 0.94), legend.title = element_blank()) +
  
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        
        plot.margin = unit(c(0, -0.5, 0, -0.5), "cm"))


####
p2 <- ggplot(data_list[[2]], aes(x=x, y=y)) + coord_fixed(ratio = 47/50) +
  
  annotation_custom(court, -250, 250, -50, 420) +
  
  geom_point(alpha = 0.5, size = 1, color="blue") + xlim(-250, 250) +
  
  ylim(-50, 420) +
  
  annotate("text", x = 160, y = 400, label="Durant", size=8)+
  
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),
        
        panel.grid.major = element_blank(), legend.justification = c(0, 1),
        
        legend.position = c(0.05, 0.94), legend.title = element_blank()) +
  
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        
        plot.margin = unit(c(0, -0.5, 0, -0.5), "cm"))

####
# harden
p3 <- ggplot(data_list[[3]], aes(x=x, y=y)) + coord_fixed(ratio = 47/50) +
  
  annotation_custom(court, -250, 250, -50, 420) +
  
  geom_point(alpha = 0.5, size = 1, color="blue") + xlim(-250, 250) +
  
  ylim(-50, 420) +
  
  annotate("text", x = 160, y = 400, label="Harden", size = 8)+
  
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),
        
        panel.grid.major = element_blank(), legend.justification = c(0, 1),
        
        legend.position = c(0.05, 0.94), legend.title = element_blank()) +
  
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        
        plot.margin = unit(c(0, -0.5, 0, -0.5), "cm"))

####
# James
p4 <- ggplot(data_list[[4]], aes(x=x, y=y)) + coord_fixed(ratio = 47/50) +
  
  annotation_custom(court, -250, 250, -50, 420) +
  
  geom_point(alpha = 0.5, size = 1, color="blue") + xlim(-250, 250) +
  
  ylim(-50, 420) +
  
  annotate("text", x = 160, y = 400, label="DeRozan", size = 8)+
  
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),
        
        panel.grid.major = element_blank(), legend.justification = c(0, 1),
        
        legend.position = c(0.05, 0.94), legend.title = element_blank()) +
  
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        
        plot.margin = unit(c(0, -0.5, 0, -0.5), "cm"))


p <- grid.arrange(p1, p2, p3, p4, nrow=1)

ggsave("EDAplot.pdf", plot = p, width = 50/10, height = 47/10)