library(spatstat)
library(data.table)

W <- owin(c(0,1000),c(0,500))

df <- fread("BCI_3specs.csv")
bci <- ppp(df$x, df$y, window=W, marks=df$Latin)

AD <- unmark(subset.ppp(bci, marks(bci) == "Acalypha diversifolia"))
CI <- unmark(subset.ppp(bci, marks(bci) == "Cecropia insignis"))
CS <- unmark(subset.ppp(bci, marks(bci) == "Cupania seemannii"))


# Model fitting
model <- ppm(AD, trend=~ polynom(x, y, 2), interaction=Poisson(), correction="isotropic")
unname(exp(coef(model)))

bei.m <- ppm(bei ~ elev + grad, data = bei.extra, covariates=list(elev=elevim, grad=gradim))
asympcov <- vcov(bei.m)
sqrt(diag(asympcov))
rho <- predict.ppm(bei.m)
Kbei <- Kinhom(bei,rho)



# berman.test.ppm(bei.m, bei.extra$elev)


mfit <- function(data){
  model <- ppm(data, trend=~ polynom(x, y, 2), interaction=Poisson(), correction="isotropic")
  K <- Kest(data, var.approx = TRUE)
  PCF <- pcfinhom(data)
  Kmod <- Kmodel.ppm(model)
  mat <- matrix(nrow = 99, ncol = length(PCF$r))
  for (i in 1:99) {
    sims <- simulate(model)
    p <- pcfinhom(sims$`Simulation 1`)$iso
    for (t in PCF$r) {
      mat[i, which(PCF$r == t)] <- p[which(PCF$r == t)]
    }
  }
  
  ulist <- apply(mat, 2, max, simplify = TRUE)
  llist <-apply(mat, 2, min, simplify = TRUE)
  plot(PCF$r, PCF$iso, ylim=c(0,5), type='l', col='blue')
  lines(PCF$r, ulist, ylim=c(0,5), type='l', col='red')
  lines(PCF$r, llist, ylim=c(0,5), type='l', col='green')
}


# Near distance aggregation, not a lot of pairs of events at large dist -> var is big
# Poisson not good fit, would have to explore other models
# pcf = K' / 2pir, should be 1 if stationary

# lines(K$r, Kmod(K$r), type='l', xlim = c(0,20), ylim = c(0, 4000))
