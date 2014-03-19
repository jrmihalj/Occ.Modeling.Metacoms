get_Rhat <- function (D, family = NA, scaling = 1.5) 
{
  if (attributes(D)$nChains < 2) {
    stop("At least two chains are required")
  }
  if (!is.na(family)) {
    D <- get_family(D, family = family)
  }
  psi.dot <- ddply(D, .(Parameter, Chain), summarize, psi.dot = mean(value), 
                   .parallel = attributes(D)$parallel)
  psi.j <- ddply(D, .(Parameter), summarize, psi.j = mean(value), 
                 .parallel = attributes(D)$parallel)
  b.df <- merge(psi.dot, psi.j)
  attr(b.df, "nIterations") <- attributes(D)$nIterations
  b.df <- cbind(b.df, nIterations = attributes(D)$nIterations)
  B <- ddply(b.df, .(Parameter), summarize, B = var(psi.j - 
                                                      psi.dot) * nIterations, .parallel = attributes(D)$parallel)
  s2j <- ddply(D, .(Parameter, Chain), summarize, s2j = var(value), 
               .parallel = attributes(D)$parallel)
  W <- ddply(s2j, .(Parameter), summarize, W = mean(s2j), .parallel = attributes(D)$parallel)
  BW <- merge(B, W)
  BW <- cbind(BW, nIterations = attributes(D)$nIterations)
  BW <- ddply(BW, .(Parameter), transform, wa = ((((nIterations - 
                                                      1)/nIterations) * W) + ((1/nIterations) * B)), .parallel = attributes(D)$parallel)
  BW <- ddply(BW, .(Parameter), transform, Rhat = sqrt(wa/W), 
              .parallel = attributes(D)$parallel)
  return(unique(BW$Rhat))
}
