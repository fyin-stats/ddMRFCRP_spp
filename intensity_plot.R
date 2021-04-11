library(ggplot2)
library(lattice)
setwd("~/Dropbox/GAM_Point Process/earthquake_data")

load("100_100Dahl.Rdata")
Dahl.lambda <- Dahl$lambda[Dahl$z]

x <- y <- seq(0, 1, length.out = 100)
data <- expand.grid(x = x, y = y)
data$z <- Dahl.lambda

pdf("100_100_Dahl_intensity.pdf")
levelplot(z~x*y,data, col.regions = terrain.colors(100)[length(terrain.colors(100)):1], 
          xlab = "X", ylab = "Y",
          scales = list(x = list(at = seq(0, 1, length.out = 5)), y = list(at = seq(0, 1, length.out = 5))))
dev.off()

load("100_100results.Rdata")
lambda.matrix <- sapply(5001:20000, function(i) res$lambda[[i]][res$z[i, ]])
lambda.mean <- rowMeans(lambda.matrix)

data1 <- expand.grid(x = x, y = y)
data1$z <- lambda.mean
pdf("100_100_mean_intensity.pdf")
levelplot(z~x*y,data1, col.regions = terrain.colors(100)[length(terrain.colors(100)):1], 
          xlab = "X", ylab = "Y",
          scales = list(x = list(at = seq(0, 1, length.out = 5)), y = list(at = seq(0, 1, length.out = 5))))
dev.off()


#ggplot(data, aes(x = x, y = y)) + geom_tile(aes(fill = z)) + 
#  scale_fill_gradientn(colours = terrain.colors(256))

#levelplot(data, col.regions = terrain.colors(100)[length(terrain.colors(100)):1])
#heatmap(data1, Rowv = NA, Colv = NA, col = terrain.colors(256))