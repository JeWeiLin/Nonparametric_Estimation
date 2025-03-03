##########Multivariate kernel regression
##########Exercise 1.
ker.est <- function(X, y, h, x0) {
  k0 <- function(u) {
    return(dnorm(u)) 
  }
  kernel_func <- function(x) {
    return(apply(x, 1, function(d) prod(sapply(d, k0))))
  }
  
  weights <- kernel_func((x0 - X) / h)
  
  numerator <- sum(y * weights)
  denominator <- sum(weights)
  
  f_hat <- numerator / denominator
  return(f_hat)
}


X <- matrix(rnorm(20), ncol = 2)  
y <- rnorm(10)                   
h <- c(1, 1)                     
x0 <- c(0, 0)               

f_x0 <- ker.est(X, y, h, x0)


set.seed(1)
f <- function(x1,x2){ 
  dnorm(x1-0.5, sd=0.2)*dnorm(x2-0.5, sd=0.2) 
}
n <- 1000
X <- matrix(runif(n*2), n,2)
y <- f(X[,1],X[,2]) + rnorm(n,sd=0.4)

h0=0.05
xlist <- (1:20)/21
ylist <- (1:10)/11
n1 <- length(xlist)
n2 <- length(ylist)
zm <- matrix(0, n1, n2)
for (i in 1:n1){  
  for (j in 1:n2){  
    zm[i,j] <- f(xlist[i], ylist[j])  
  } 
}
f.persp <- persp(xlist, ylist, zm, theta=20)

for (i in 1:n1){
  xi <- xlist[i]
  fhat.xi <- rep(0, n2)
  for (j in 1:n2){
    fhat.xi[j] <- ker.est(X, y, rep(h0,2), c(xi,ylist[j]))
  }
  lines(trans3d(xi, ylist, fhat.xi, pmat=f.persp), col=2)
}
dif1 <- function(u, v){ (f(u,v) - ker.est(X, y, rep(h0,2), c(u,v)))^2 }

tem1 <- function(u){
  tem2 <- function(v){ dif1(u,v)}
  vtem2 <- Vectorize(tem2)
  return(integrate(vtem2, 0, 1)$value)
}
vtem1 <- Vectorize(tem1)
integrate(vtem1, 0,1)$value
#h0 = 0.05 => ISE:  0.009570552




##########Exercise 2.
set.seed(1)
f <- function(x) { dnorm(x - 0.5, sd = 0.2) } # True regression function
n <- 1000
X <- matrix(runif(n))       
y <- f(X) + rnorm(n, sd = 0.4)  


ker.est.bias_corrected <- function(X, y, h, x0) {
  
  k0 <- function(u) { dnorm(u) }
  
  kernel_func <- function(x) {
    return(apply(x, 1, function(d) prod(sapply(d, k0))))
  }
  
  weights <- kernel_func((x0 - X) / h)
  
  numerator <- sum(y * weights)
  denominator <- sum(weights)
  
  f_hat <- numerator / denominator
  
  kernel_weights <- sapply((x0 - X) / h, k0)
  numerator_corrected <- sum(y * kernel_weights * ((x0 - X) / h)^2)
  f_hat_corrected <- f_hat - h^2 * numerator_corrected / (2 * denominator)
  
  return(f_hat_corrected)
}

h <- 0.05


x_grid <- seq(0, 1, length.out = 100)
f_true <- sapply(x_grid, f)
f_estimated <- sapply(x_grid, function(x0) ker.est.bias_corrected(X, y, h, x0))


ISE <- sum((f_estimated - f_true)^2) * (1 / length(x_grid))
ISE
