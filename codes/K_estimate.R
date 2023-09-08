library(spatstat)
Kfunc <- function(data){
  K <- Kest(data, var.approx = TRUE)
  plot(K$r,K$iso, type='l')
  lines(K$r,K$theo, type='l', col='blue', lty=2)
  lines(K$r,K$theo+2*sqrt(K$ls), col='red')
  lines(K$r,K$theo-2*sqrt(K$ls), col='green')
  legend(x='topleft',c("data","CSR","CSR +2 standard error", "CSR -2 standard error"),cex=.8,
         col=c("black","blue", "red", "green"),pch=c(1,2,3,4))
}

subbei <- subset.ppp(bei, x<=100 & y>=400)
subbei <- ppp(subbei$x, subbei$y, window=owin(xrange=c(0,100),yrange=c(400,500)))

# bci
library(data.table)
# library(ecespa)
W <- owin(c(0,1000),c(0,500))

# bci <- scanpp("BCI_3specs.txt", window=W, header=TRUE)
df <- fread("BCI_3specs.csv")
# ad <- subset(df, Latin="Acalypha diversifolia")
# ci <- subset(df, Lain="Cecropia insignis")
# cs <- subset(df, Lain="Cupania seemannii")
ad <- function(X) {marks(X) == "Acalypha diversifolia"}
ci <- function(X) {marks(X) == "Cecropia insignis"}
cs <- function(X) {marks(X) == "Cupania seemannii"}
bci <- ppp(df$x, df$y, window=W, marks=df$Latin)
plot(bci)
K1 <- Kmulti(bci, ad, ci, correction="isotropic")
K2 <- Kmulti(bci, ci, cs, correction="isotropic")
K3 <- Kmulti(bci, ad, cs, correction="isotropic")

# var1 <- Kmulti.ls(bci, I=ad(bci), J=ci(bci), r = NULL, corre = "isotropic")

AD <- subset.ppp(bci, marks(bci) == "Acalypha diversifolia")
CI <- subset.ppp(bci, marks(bci) == "Cecropia insignis")
CS <- subset.ppp(bci, marks(bci) == "Cupania seemannii")

b <- function(t, A, P){pi*t^2 / A *(1 - pi*t^2 / A) + A^(-2)*(1.0716*P*t^2 + 2.2375*t^4)}
a_1 <- function(t, A, P){A^(-2)*(0.21*P*t^3 + 1.3*t^4)}
a_2 <- function(t, A, P){A^(-3)*(0.24*P*t^5 + 2.62*t^6)}

v <- function(n, t, A, P){(2*b(t, A, P) - a_1(t, A, P) + (n-2)*a_2(t, A, P))/(n*(n-1))}

var <- function(t, data1, data2){
  n1 <- data1$n
  n2 <- data2$n
  A <- area(window(data1))
  P <- perimeter(window(data1))
  c <- n2/(n1+n2)
  A^2/(n1*n2)*(b(t, A, P) - 2*c*(1-c)*a_1(t, A, P) + ((n1-1)*c^2 + (n2-1)*(1-c)^2)*a_2(t, A, P))
}

plt <- function(K, d1, d2){
  plot(K$r,K$iso, type='l', col='cyan')
  lines(K$r,K$theo, type='l', col='blue', lty=2)
  lines(K$r,K$theo+2*sqrt(var(K$r, d1, d2)), col='red')
  lines(K$r,K$theo-2*sqrt(var(K$r, d1, d2)), col='green')
  legend(x='topleft',c("data","CSR","CSR +2 standard error", "CSR -2 standard error"),cex=.8,
         col=c("cyan","blue", "red", "green"),pch=c(1,2,3,4))
}
