###function to infer a sparse GGM using vbsr
spikeGGM = function(x,
                    ...){
  neighborhoodSelection = function(i,
                                   x,
                                   ...){

    #run individual vbsr's
    vbsrResult = vbsr::vbsr(y=x[,i],
                            X=x[,-i],
                            ...)

    #store and return appropriate results
    vec <- rep(0,ncol(x))
    vec[-i] <- vbsrResult$z
    return(vec)
  }

  ind = 1:ncol(x)

  zmatrix = sapply(ind,
                  neighborhoodSelection,
                  x,
                  ...)

  zmatrix <- zmatrix/2 + t(zmatrix)/2

  return(pchisq(zmatrix^2,1,lower.tail=F))

}
