#simulate test data
library(dplyr)
n = 101
p = 100
theta = 2/p
set.seed(47)
LAM = (rnorm(p^2)*rbinom(p^2,1,theta)) %>%
  matrix(p,p)
diag(LAM) <- 1
OM <- LAM%*%t(LAM)
E <- diag(.5,p)
Y <- MASS::mvrnorm(n=n,mu=rep(0,p),solve(OM))
#pheatmap::pheatmap(scale(Y))
#pheatmap::pheatmap(cor((Y)))
#pheatmap::pheatmap(cov2cor(solve(OM)))
res <- spikeGGM::combinedGGMTopology(Y,
                                     adjmethod='fdr',
                                     nrestarts=10)
table(res$lasso$invcovEst!=0,OM!=0)/2
table(res$spike$invcovEst!=0,OM!=0)/2
table(res$combined$invcovEst!=0,OM!=0)/2

#pairs(cbind(c(res$lasso$covEst),c(res$spike$covEst),c(res$combined$covEst),c(cov(Y)),c(solve(OM))))

cor(cbind(c(res$lasso$covEst),c(res$spike$covEst),c(res$combined$covEst),c(cov(Y)),c(solve(OM))))

#pairs(cbind(c(res$lasso$invcovEst),c(res$spike$invcovEst),c(res$combined$invcovEst),c(solve(cov(Y))),c((OM))))

cor(cbind(c(res$lasso$invcovEst),c(res$spike$invcovEst),c(res$combined$invcovEst),c(solve(cov(Y))),c((OM))))
