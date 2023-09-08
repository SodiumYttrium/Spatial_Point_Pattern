library(spatstat)

# H function on vesicles (irregular)

vs <- vesicles
n = 37
data <- vs

dist <- seq(from = 0, to = 905, by = 5)

Hfunc <- function(data, dist){
  n <- data$n
  pd <- pairdist(data, periodic=FALSE, method="C", squared=FALSE, metric=NULL)
  Hhat = function(t){sum(pd[lower.tri(pd)] <= t) / (n*(n-1)/2)}
  
  Hhat_apply <- lapply(dist, Hhat)
  
  pdlist <- vector("list", length = length(dist))
  pdmat <- matrix(nrow = 99, ncol = length(dist))
  for (i in 1:99) {
    pts <- runifpoint(n, win=window(data))
    pds <- pairdist(pts, periodic=FALSE, method="C", squared=FALSE, metric=NULL)
    pdsvec <- pds[lower.tri(pds)]
    for (t in dist) {
      pdmat[i, which(dist == t)] <- sum(pdsvec <= t)
    }
  }
  
  Hbar <- apply(pdmat, 2, mean) / (n*(n-1)/2)
  
  ulist <- apply(pdmat, 2, max, simplify = TRUE)
  llist <-apply(pdmat, 2, min, simplify = TRUE)
  U <- lapply(ulist, function(t) t/((n*(n-1))/2))
  L <- lapply(llist, function(t) t/((n*(n-1))/2))
  
  plot(Hbar, Hhat_apply, col='green', type='l', xlab='H(t)')
  lines(Hbar, U, col='red')
  lines(Hbar, L, col='cyan')
  legend(x=0,y=1,c("data","upper envelope", "lower envelope"),cex=.8,col=c("green","red", "cyan"),pch=c(1,2,3))
}

# H on redwood3

rw <- redwood3
n <- 62
rwpd <- pairdist(rw, periodic=FALSE, method="C", squared=FALSE, metric=NULL)

Hhat = function(t){sum(rwpd[lower.tri(rwpd)] <= t) / (n*(n-1)/2)}
H = function(t){
  if (0 <= t && t <= 1) {
    pi*t^2 - 8/3*t^3 + t^4/2
  } else if (1 < t && t <= 1.42) {
    1/3 - 2*t^2 - t^4/2 + 4*(t^2 - 1)^(1/2)*(2*t^2 + 1)/3 +2*t^2*asin(2*t^(-2) - 1)
  }
}

dist <- seq(from = 0, to = 1.42, by = 0.01)

H_apply <- lapply(dist, H)
Hhat_apply <- lapply(dist, Hhat)

pdlist <- vector("list", length = length(dist))
pdmat <- matrix(nrow = 99, ncol = length(dist))
for (i in 1:99) {
  pts <- runifpoint(n, win=window(rw))
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

# H on cells
m <- 42
cpd <- pairdist(cells, periodic=FALSE, method="C", squared=FALSE, metric=NULL)

Hhat = function(t){sum(cpd[lower.tri(cpd)] <= t) / (m*(m-1)/2)}
H = function(t){
  if (0 <= t && t <= 1) {
    pi*t^2 - 8/3*t^3 + t^4/2
  } else if (1 < t && t <= 1.42) {
    1/3 - 2*t^2 - t^4/2 + 4*(t^2 - 1)^(1/2)*(2*t^2 + 1)/3 +2*t^2*asin(2*t^(-2) - 1)
  }
}
H_apply <- lapply(dist, H)
Hhat_apply <- lapply(dist, Hhat)

pdlist <- vector("list", length = length(dist))
pdmat <- matrix(nrow = 99, ncol = length(dist))
for (i in 1:99) {
  pts <- runifpoint(m, win=window(cells))
  pds <- pairdist(pts, periodic=FALSE, method="C", squared=FALSE, metric=NULL)
  pdsvec <- pds[lower.tri(pds)]
  for (t in dist) {
    pdmat[i, which(dist == t)] <- sum(pdsvec <= t)
  }
}
ulist <- apply(pdmat, 2, max, simplify = TRUE)
llist <-apply(pdmat, 2, min, simplify = TRUE)
U <- lapply(ulist, function(t) t/((m*(m-1))/2))
L <- lapply(llist, function(t) t/((m*(m-1))/2))

plot(H_apply, Hhat_apply, col='green', type='l', xlab='H(t)')
lines(H_apply, U, col='red')
lines(H_apply, L, col='cyan')
legend(x=0,y=1,c("data","upper envelope", "lower envelope"),cex=.8,col=c("green","red", "cyan"),pch=c(1,2,3))

# H on bei
subbei <- subset.ppp(bei, x<=100 & y>=400)
subbei <- ppp(subbei$x, subbei$y, window=owin(xrange=c(0,100),yrange=c(400,500)))
dist <- seq(from = 0, to = 125, by = 1)
Hfunc(bei, seq(from = 0, to = 1100, by = 10))

# G function
Gfunc <- function(data, dist){
  n <- data$n
  Y <- vector("list", length = n)
  pd <- pairdist(data, periodic=FALSE, method="C", squared=FALSE, metric=NULL)
  for (i in 1:n) {
    Y[[i]] <- min(pd[i, -i])
  }
  A <- area(window(data))
  G  = function(y){1 - (1 - pi*y^2/A)^(n-1)}
  Ghat = function(y){sum(Y <= y)/n}
  G_apply <- lapply(dist, G)
  Ghat_apply <- lapply(dist, Ghat)
  
  gpdmat <- matrix(nrow = 99, ncol = length(dist))
  
  for (i in 1:99) {
    gpdsvec <- vector("list", length = n)
    gpts <- runifpoint(n, win=window(data))
    gpds <- pairdist(gpts, periodic=FALSE, method="C", squared=FALSE, metric=NULL)
    diag(gpds) <- diag(gpds) + 1000
    gpdsvec[[i]] <- apply(gpds, 2, min, simplify = TRUE)
    for (y in dist) {
      gpdmat[i, which(dist == y)] <- sum(unlist(gpdsvec) <= y)
    }
  }
  gulist <- apply(gpdmat, 2, max, simplify = TRUE)
  gllist <-apply(gpdmat, 2, min, simplify = TRUE)
  gU <- lapply(gulist, function(t) t/n)
  gL <- lapply(gllist, function(t) t/n)
  
  plot(G_apply, Ghat_apply, col='green', type='l', xlim=c(0,1))
  lines(G_apply, gU, col='red')
  lines(G_apply, gL, col='cyan')
  legend(x=0,y=1,c("data","upper envelope", "lower envelope"),cex=.8,col=c("green","red", "cyan"),pch=c(1,2,3))
}

# G for redwood3
Y <- vector("list", length = n)
for (i in 1:n) {
  Y[[i]] <- min(rwpd[i, -i])
}
A <- 1 # I assume the area is 1
G  = function(y){1 - (1 - pi*y^2/A)^(n-1)}
Ghat = function(y){sum(Y <= y)/n}

dist <- seq(from = 0, to = 0.54, by = 0.01)
G_apply <- lapply(dist, G)
Ghat_apply <- lapply(dist, Ghat)

gpdmat <- matrix(nrow = 99, ncol = length(dist))

for (i in 1:99) {
  gpdsvec <- vector("list", length = n)
  gpts <- runifpoint(n, win=window(rw))
  gpds <- pairdist(gpts, periodic=FALSE, method="C", squared=FALSE, metric=NULL)
  diag(gpds) <- diag(gpds) + 1000
  gpdsvec[[i]] <- apply(gpds, 2, min, simplify = TRUE)
  for (y in dist) {
    gpdmat[i, which(dist == y)] <- sum(unlist(gpdsvec) <= y)
  }
}
gulist <- apply(gpdmat, 2, max, simplify = TRUE)
gllist <-apply(gpdmat, 2, min, simplify = TRUE)
gU <- lapply(gulist, function(t) t/n)
gL <- lapply(gllist, function(t) t/n)

plot(G_apply, Ghat_apply, col='green', type='l', xlim=c(0,1))
lines(G_apply, gU, col='red')
lines(G_apply, gL, col='cyan')
legend(x=0,y=1,c("data","upper envelope", "lower envelope"),cex=.8,col=c("green","red", "cyan"),pch=c(1,2,3))

# G on cells
Y <- vector("list", length = m)
for (i in 1:m) {
  Y[[i]] <- min(cpd[i, -i])
}
A <- 1 # I assume the area is 1
G  = function(y){1 - (1 - pi*y^2/A)^(m-1)}
Ghat = function(y){sum(Y <= y)/m}

dist <- seq(from = 0, to = 0.54, by = 0.01)
G_apply <- lapply(dist, G)
Ghat_apply <- lapply(dist, Ghat)

gpdmat <- matrix(nrow = 99, ncol = length(dist))

for (i in 1:99) {
  gpdsvec <- vector("list", length = m)
  gpts <- runifpoint(m, win=window(cells))
  gpds <- pairdist(gpts, periodic=FALSE, method="C", squared=FALSE, metric=NULL)
  diag(gpds) <- diag(gpds) + 1000
  gpdsvec[[i]] <- apply(gpds, 2, min, simplify = TRUE)
  for (y in dist) {
    gpdmat[i, which(dist == y)] <- sum(unlist(gpdsvec) <= y)
  }
}
gulist <- apply(gpdmat, 2, max, simplify = TRUE)
gllist <-apply(gpdmat, 2, min, simplify = TRUE)
gU <- lapply(gulist, function(t) t/m)
gL <- lapply(gllist, function(t) t/m)

plot(G_apply, Ghat_apply, col='green', type='l', xlim=c(0,1))
lines(G_apply, gU, col='red')
lines(G_apply, gL, col='cyan')
legend(x=0,y=1,c("data","upper envelope", "lower envelope"),cex=.8,col=c("green","red", "cyan"),pch=c(1,2,3))

# F function
Ffunc <- function(data, x, y, dist){
  n <- data$n
  A <- area(window(data))
  F <- function(z){1-exp(-pi*n*z^2/A)}
  
  x_coords <- seq(from = window(data)$window$xrange[1], to = window(data)$window$xrange[2],
                  by = (window(data)$window$xrange[2] - window(data)$window$xrange[1])/(x-1))
  y_coords <- seq(from = window(data)$window$yrange[1], to = window(data)$window$yrange[2],
                  by = (window(data)$window$yrange[2] - window(data)$window$yrange[1])/(y-1))
  d_grid <- expand.grid(x = x_coords, y = y_coords)
  d_grid <- as.ppp(d_grid, W = window(data))
  
  dpdist <- nncross(d_grid, data, what='dist')
  
  Fhat <- function(z){sum(dpdist <= z)/d_grid$n}
  
  F_apply <- lapply(dist, F)
  Fhat_apply <- lapply(dist, Fhat)
  
  fpdmat <- matrix(nrow = 99, ncol = length(dist))
  for (i in 1:99) {
    fpts <- runifpoint(n, win=window(data))
    fpds <- nncross(d_grid, fpts, what='dist')
    for (t in dist) {
      fpdmat[i, which(dist == t)] <- sum(fpds <= t)
    }
  }
  
  fulist <- apply(fpdmat, 2, max, simplify = TRUE)
  fllist <-apply(fpdmat, 2, min, simplify = TRUE)
  fU <- lapply(fulist, function(t) t/d_grid$n)
  fL <- lapply(fllist, function(t) t/d_grid$n)
  
  plot(F_apply, Fhat_apply, col='green', type='l', xlim=c(0,1), ylim=c(0,1))
  lines(F_apply, fU, col='red')
  lines(F_apply, fL, col='cyan')
  legend(x=0,y=1,c("data","upper envelope", "lower envelope"),cex=.8,col=c("green","red", "cyan"),pch=c(1,2,3))
}

# F on redwood3
F <- function(z){1-exp(-pi*n*z^2/A)}

x_coords <- seq(from = 0, to = 1, by = 1/15)
y_coords <- seq(from = -1, to = 0, by = 1/15)
rw_grid <- expand.grid(x = x_coords, y = y_coords)
rw_grid <- as.ppp(rw_grid, W = window(rw))

dpdist <- nncross(rw_grid, rw, what='dist')

dist <- seq(from = 0, to = 1.15, by = 0.01)

Fhat <- function(z){sum(dpdist <= z)/rw_grid$n}

F_apply <- lapply(dist, F)
Fhat_apply <- lapply(dist, Fhat)

fpdmat <- matrix(nrow = 99, ncol = length(dist))
for (i in 1:99) {
  fpts <- runifpoint(n, win=window(rw))
  fpds <- nncross(rw_grid, fpts, what='dist')
  for (t in dist) {
    fpdmat[i, which(dist == t)] <- sum(fpds <= t)
  }
}

fulist <- apply(fpdmat, 2, max, simplify = TRUE)
fllist <-apply(fpdmat, 2, min, simplify = TRUE)
fU <- lapply(fulist, function(t) t/rw_grid$n)
fL <- lapply(fllist, function(t) t/rw_grid$n)

plot(F_apply, Fhat_apply, col='green', type='l', xlim=c(0,1), ylim=c(0,1))
lines(F_apply, fU, col='red')
lines(F_apply, fL, col='cyan')
legend(x=0,y=1,c("data","upper envelope", "lower envelope"),cex=.8,col=c("green","red", "cyan"),pch=c(1,2,3))

# F on cells
F <- function(z){1-exp(-pi*m*z^2/A)}

x_coords <- seq(from = 0, to = 1, by = 1/13)
y_coords <- seq(from = 0, to = 1, by = 1/13)
c_grid <- expand.grid(x = x_coords, y = y_coords)
c_grid <- as.ppp(c_grid, W = window(cells))

dpdist <- nncross(c_grid, cells, what='dist')

dist <- seq(from = 0, to = 0.25, by = 0.005)

Fhat <- function(z){sum(dpdist <= z)/c_grid$n}

F_apply <- lapply(dist, F)
Fhat_apply <- lapply(dist, Fhat)

fpdmat <- matrix(nrow = 99, ncol = length(dist))
for (i in 1:99) {
  fpts <- runifpoint(m, win=window(cells))
  fpds <- nncross(c_grid, fpts, what='dist')
  for (t in dist) {
    fpdmat[i, which(dist == t)] <- sum(fpds <= t)
  }
}

fulist <- apply(fpdmat, 2, max, simplify = TRUE)
fllist <-apply(fpdmat, 2, min, simplify = TRUE)
fU <- lapply(fulist, function(t) t/c_grid$n)
fL <- lapply(fllist, function(t) t/c_grid$n)

plot(F_apply, Fhat_apply, col='green', type='l', xlim=c(0,1), ylim=c(0,1))
lines(F_apply, fU, col='red')
lines(F_apply, fL, col='cyan')
legend(x=0,y=1,c("data","upper envelope", "lower envelope"),cex=.8,col=c("green","red", "cyan"),pch=c(1,2,3))

# F on bei
# x = 120, y = 60 ?
# dist = 0 to 120 by 10
Ffunc(bei, 120, 60, seq(from = 0, to = 119, by = 1))

# Quadrat on bei
quadratcount(bei, nx = 38, ny = 19)
quadrat.test(bei, nx=38, ny=19, alternative='clustered')