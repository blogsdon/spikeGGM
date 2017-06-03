getGGM = function(S,edgeList){
  library(dplyr)
  res <- glasso::glasso(S,rho=0,zero=edgeList)
  list(covEst = res$w,
       invcovEst = res$wi) %>%
    return()
}
