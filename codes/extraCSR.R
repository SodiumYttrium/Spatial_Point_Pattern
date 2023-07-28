# H function on vesicles (irregular)

vs <- vesicles
n = 37

pd <- pairdist(vs, periodic=FALSE, method="C", squared=FALSE, metric=NULL)
Hhat = function(t){sum(pd[lower.tri(pd)] <= t) / (n*(n-1)/2)}

dist <- seq(from = 0, to = 905, by = 5)

Hhat_apply <- lapply(dist, Hhat)

pdlist <- vector("list", length = length(dist))
pdmat <- matrix(nrow = 99, ncol = length(dist))
for (i in 1:99) {
  pts <- runifpoint(n, win=window(vs))
  pds <- pairdist(pts, periodic=FALSE, method="C", squared=FALSE, metric=NULL)
  pdsvec <- pds[lower.tri(pds)]
  for (t in dist) {
    pdmat[i, which(dist == t)] <- sum(pdsvec <= t)
  }
}

Hbar <- apply(pdmat, 2, mean) / (n*(n-1)/2)

ulist <- RowMax(t(pdmat))
llist <-apply(pdmat, 2, min, simplify = TRUE)
U <- lapply(ulist, function(t) t/((n*(n-1))/2))
L <- lapply(llist, function(t) t/((n*(n-1))/2))

plot(Hbar, Hhat_apply, col='green', type='l', xlab='H(t)')
lines(Hbar, U, col='red')
lines(Hbar, L, col='cyan')
legend(x=0,y=1,c("data","upper envelope", "lower envelope"),cex=.8,col=c("green","red", "cyan"),pch=c(1,2,3))
abline(0, 1, col='blue')