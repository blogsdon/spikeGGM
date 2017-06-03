combinedGGMTopology = function(Y,
                               adjmethod='fdr',
                               seed=47,
                               nrestarts=10){
  library(dplyr)
  S <- cov(Y)
  Y <- scale(Y)
  p <- ncol(Y)
  set.seed(seed)

  cat('computing spike ggm...\n')
  spikeResult <- spikeGGM::spikeGGM(Y,
                          n_orderings=nrestarts,
                          cleanSolution=TRUE)

  cat('computing lasso ggm...\n')
  lassoResult <- FastGGM::FastGGM(Y)

  cat('combing and adjusting...\n')
  fdrAdjust1 = spikeResult[which(upper.tri(spikeResult))] %>%
    p.adjust(method=adjmethod)
  fdrAdjust2 = lassoResult$p_precision[which(upper.tri(lassoResult$p_precision))] %>%
    p.adjust(method=adjmethod)

  spikeMat = spikeResult
  lassoMat = lassoResult$p_precision

  spikeMat[upper.tri(spikeMat)] <- fdrAdjust1
  spikeMat <- t(spikeMat)
  spikeMat[upper.tri(spikeMat)] <- fdrAdjust1
  diag(spikeMat) <- 0

  lassoMat[upper.tri(lassoMat)] <- fdrAdjust2
  lassoMat <- t(lassoMat)
  lassoMat[upper.tri(lassoMat)] <- fdrAdjust2
  diag(lassoMat) <- 0

  combinedVec <- cbind(c(lassoMat),c(spikeMat)) %>%
    apply(1,min)

  combinedMat <- matrix(combinedVec,p,p)

  cat('fitting ggm given topologies...\n')
  #S <- cov(Y)
  lassoSol <- spikeGGM::getGGM(S,which(lassoMat>=0.05,T))
  spikeSol <- spikeGGM::getGGM(S,which(spikeMat>=0.05,T))
  combSol <- spikeGGM::getGGM(S,which(combinedMat>=0.05,T))

  list(lasso = lassoSol,
       spike = spikeSol,
       combined = combSol) %>%
    return

}
