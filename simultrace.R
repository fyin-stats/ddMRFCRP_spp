
#load(FinalResult.RData)
load("FinalResult.RData")
library(tidyverse)
totalrep = 50
fit.rep = list()
tbar <- txtProgressBar(min = 0, max = 50, style = 3)
for (i in 1:totalrep) {
  set.seed(i)
  fit.rep[[i]] = CDMFM(data = A, data1 = AAA, lambda1 = 1, neighbour = 20,
                       niterations = 1000, mu_0 = 5, mu_0off = -5,
                       t_0 = 2, alpha = 1, beta = 1, GAMMA = 1, LAMBDA = 1,
                       initNClusters = 9)
  #result.rep = getDahl(fit.rep, 500)
  #RI.single[i]=rand.index(asm,result.rep$zout)
  setTxtProgressBar(tbar, i)
}



# 
# sim <- function(seed) {
#   set.seed(seed)
#   temp <- CDMFM(data = A, data1 = AAA, lambda1 = 1, neighbour = 20,
#                niterations = 1000, mu_0 = 5, mu_0off = -5, t_0 = 2,
#                alpha = 1, beta = 1, GAMMA=1, LAMBDA = 1, initNClusters = 9)
#   return(temp)
# }
# 
# library(furrr)
# plan(multisession(workers = 10))
# fit.rep <- future_map(1:totalrep, ~sim(.x), .progress = TRUE)
# 



RI.rep=matrix(0,1000,totalrep)
for(j in 1:totalrep){
  for(i in 1:1000){
    RI.rep[i,j]=fossil::rand.index(asm,fit.rep[[j]]$Iterates[[i]]$zout)
  }
}

df <- data.frame(method=c(rep("MFM",1000)),
                 iteration=c(1:1000),
                 RI=apply(RI.rep,1,mean))
RIplot=ggplot(data=df, aes(x=iteration, y=RI, group=method)) +
  geom_line(aes(color=method))+
  labs(x="Iteration", y = "Rand Index")+ theme(legend.position = c(0.8, 0.2))


df2 <- data.frame(method=c(rep("REP",1000*50)),
                  iteration=c(rep(1:1000,50)),
                  RI=c(RI.rep))
RIplot <- ggplot(data=df2, aes(x=iteration, y=RI, group=method)) +
  geom_line(aes(color=method),colour="darkgrey")+
  labs(x="Iteration", y = "Rand Index")+ 
  geom_line(data = df,aes(x = iteration,y = RI,group = method),colour = "red")+
  theme(legend.position = c(0.8, 0.2),
        text = element_text(size = 16)) + theme_bw()

ggsave(plot = RIplot, width = 6, height = 3,
       filename = "~/Desktop/simtrace.pdf")

saveRDS(fit.rep, file = "./fit_rep.rds")
