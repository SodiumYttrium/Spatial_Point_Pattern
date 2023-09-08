library(spatstat)
# library(attention)

jpp <- japanesepines

# plot(jpp)
n <- 65
pd <- pairdist(jpp, periodic=FALSE, method="C", squared=FALSE, metric=NULL)

# Inter-event distances

Hhat = function(t){sum(pd[lower.tri(pd)] <= t) / (n*(n-1)/2)}
H = function(t){
  if (0 <= t && t <= 1) {
    pi*t^2 - 8/3*t^3 + t^4/2
  } else if (1 < t && t <= 1.42) {
    1/3 - 2*t^2 - t^4/2 + 4*(t^2 - 1)^(1/2)*(2*t^2 + 1)/3 +2*t^2*asin(2*t^(-2) - 1)
  }
}

# pts <- rpoispp(lambda = n, win = window(jpp), nsim = 65)
# plot(pts, main = "Simulated Homogeneous Point Process", pch = 16)

dist <- seq(from = 0, to = 1.42, by = 0.01)

H_apply <- lapply(dist, H)
Hhat_apply <- lapply(dist, Hhat)

pdlist <- vector("list", length = length(dist))
pdmat <- matrix(nrow = 99, ncol = length(dist))
for (i in 1:99) {
  pts <- runifpoint(n, win=window(jpp))
  pds <- pairdist(pts, periodic=FALSE, method="C", squared=FALSE, metric=NULL)
  pdsvec <- pds[lower.tri(pds)]
  for (t in dist) {
    pdmat[i, which(dist == t)] <- sum(pdsvec <= t)
  }
}
ulist <- apply(pdmat, 2, max, simplify = TRUE) # RowMax(t(pdmat))
llist <-apply(pdmat, 2, min, simplify = TRUE)
U <- lapply(ulist, function(t) t/((n*(n-1))/2))
L <- lapply(llist, function(t) t/((n*(n-1))/2))


plot(H_apply, Hhat_apply, col='green', type='l', xlab='H(t)')
lines(H_apply, U, col='red')
lines(H_apply, L, col='cyan')
legend(x=0,y=1,c("data","upper envelope", "lower envelope"),cex=.8,col=c("green","red", "cyan"),pch=c(1,2,3))
abline(0, 1, col='blue')

# Nearest neighbour distance

Y <- vector("list", length = 65)
for (i in 1:65) {
  Y[[i]] <- min(pd[i, -i])
}
A <- 1 # I assume the area is 1
G  = function(y){1 - (1 - pi*y^2/A)^(n-1)}
Ghat = function(y){sum(Y <= y)/n}

dist <- seq(from = 0, to = 0.54, by = 0.01)
G_apply <- lapply(dist, G)
Ghat_apply <- lapply(dist, Ghat)

gpdmat <- matrix(nrow = 99, ncol = length(dist))

for (i in 1:99) {
  gpdsvec <- vector("list", length = 65)
  gpts <- runifpoint(n, win=window(jpp))
  gpds <- pairdist(gpts, periodic=FALSE, method="C", squared=FALSE, metric=NULL)
  diag(gpds) <- diag(gpds) + 1000
  # gpdsvec[[i]] <- apply(gpds[row(gpds)!=col(gpds)], 2, min, simplify = TRUE)
  gpdsvec[[i]] <- apply(gpds, 2, min, simplify = TRUE)
  # add big num to diagonal, take min
  for (y in dist) {
    gpdmat[i, which(dist == y)] <- sum(unlist(gpdsvec) <= y)
  }
}
gulist <- apply(gpdmat, 2, max, simplify = TRUE) # RowMax(t(gpdmat))
gllist <-apply(gpdmat, 2, min, simplify = TRUE)
gU <- lapply(gulist, function(t) t/n)
gL <- lapply(gllist, function(t) t/n)

plot(G_apply, Ghat_apply, col='green', type='l', xlim=c(0,1))
lines(G_apply, gU, col='red')
lines(G_apply, gL, col='cyan')
abline(0, 1, col='blue')

# Point to nearest event distances

# plot(Fest(jpp))
F <- function(z){1-exp(-pi*n*z^2/A)}

x_coords <- seq(from = 0, to = 1, by = 1/15)
y_coords <- seq(from = 0, to = 1, by = 1/15)
jpp_grid <- expand.grid(x = x_coords, y = y_coords)
jpp_grid <- as.ppp(jpp_grid, W = owin(c(0, 1), c(0, 1)))

dpdist <- nncross(jpp_grid, jpp, what='dist')

dist <- seq(from = 0, to = 0.2, by = 0.005)

Fhat <- function(z){sum(dpdist <= z)/jpp_grid$n}

F_apply <- lapply(dist, F)
Fhat_apply <- lapply(dist, Fhat)

fpdmat <- matrix(nrow = 99, ncol = length(dist))
for (i in 1:99) {
  fpts <- runifpoint(n, win=window(jpp))
  fpds <- nncross(jpp_grid, fpts, what='dist')
  for (t in dist) {
    fpdmat[i, which(dist == t)] <- sum(fpds <= t)
  }
}

fulist <- apply(fpdmat, 2, max, simplify = TRUE)
fllist <-apply(fpdmat, 2, min, simplify = TRUE)
fU <- lapply(fulist, function(t) t/jpp_grid$n)
fL <- lapply(fllist, function(t) t/jpp_grid$n)

plot(F_apply, Fhat_apply, col='green', type='l', xlim=c(0,1))
lines(F_apply, fU, col='red')
lines(F_apply, fL, col='cyan')
abline(0, 1, col='blue')

# Quadrat counts

quadrats <- quadratcount(jpp, nx = 3, ny = 3)
quadrat.test(jpp, nx=3, ny=3)