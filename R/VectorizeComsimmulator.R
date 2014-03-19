y <- array(0, dim=c(75, 12))
samp <- function(x) x<-sample(c(0,1), 1)
y <- apply(y, c(1,2), samp)

nr <- nrow(y)
nc <- ncol(y)
rs <- rowSums(y)
p <- colSums(y)
out <- matrix(0, nrow = nr, ncol = nc)

for (i in 1:nr) {
  out[i, sample.int(nc, rs[i], prob = p)] <- 1
}
