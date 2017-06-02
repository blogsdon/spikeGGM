#simulate test data
library(dplyr)
n = 300
p = 100
theta = 2/p
set.seed(47)
LAM = (rnorm(p^2)*rbinom(p^2,1,theta)) %>%
  matrix(p,p)
diag(LAM) <- 1
OM <- LAM%*%t(LAM)
E <- diag(.5,p)
Y=MASS::mvrnorm(n=n,mu=rep(0,p),solve(OM))
#pheatmap::pheatmap(scale(Y))
#pheatmap::pheatmap(cor((Y)))
#pheatmap::pheatmap(cov2cor(solve(OM)))
set.seed(1)
test=spikeGGM::spikeGGM(scale(Y),n_orderings=500,cleanSolution=TRUE)
test2=FastGGM::FastGGM(scale(Y))

fastggm1 = test2$p_precision < 0.05/choose(p,2)
#test = pchisq(test^2,1,lower.tail=F)
test5 = cbind(c(test),c(test2$p_precision)) %>%
  apply(1,min)

plot(-log10(test),-log10(test2$p_precision))
table(test<0.05/choose(p,2),OM!=0)/2
table(fastggm1,OM!=0)/2
table(test5<0.05/choose(p,2),OM!=0)/2
