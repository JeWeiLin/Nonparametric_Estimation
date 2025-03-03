##########Density estimation based on basis function approximation
##########Exercise 1.
###generate data of size 1000 (stored in x) from density f
set.seed(1)
mu1 <- 0.2
mu2 <- 0.7
n <- 10000
m <- n * 10
z <- rnorm(m, mean = mu1, sd = 0.1)
x <- z[(z > 0) & ( z< 1)]
x <- x[1:n]

z <- rnorm(m, mean = mu2, sd = 0.2)
x2 <- z[(z > 0) & (z < 1)]
x2 <- x2[1:n]
z <- sample(0:1, size = n, replace=T)
x[z==1] <- x2[z==1]

#### compute the matrix whose (i,j)th element is the integral of B_iB_j
knotlist <- (1:8) / 9 
nb <- length(knotlist) + 4
M <- matrix(0, nb, nb)
for (i in 1:nb){
  for (j in i:nb){
    tem <- function(u){
      bx <- bs(u, knots = knotlist, Boundary.knots = c(0,1), intercept=T)
      return( bx[,i]*bx[,j])
    }
    M[i,j] <- integrate(tem, 0, 1)$value
    if (j > i) { M[j,i] <- M[i,j] }
  }
}
#### compute fhat,  the estimator of f using method of moments
moments <- apply(bs(x, knots = knotlist, Boundary.knots=c(0,1), intercept=T), 2, mean)
ahat <- solve( M, moments)
fhat <- function(u){
  ans <- bs(u, knots = knotlist, Boundary.knots = c(0,1), intercept=T) %*% ahat
  return( as.numeric(ans) )
}
##### compare fhat with the true density f
k0 = pnorm(1, mean = mu1, sd = 0.1) - pnorm(0, mean = mu1, sd = 0.1)
k1 = pnorm(1, mean = mu2, sd = 0.2) - pnorm(0, mean = mu2, sd = 0.2)
f <- function(x){
  ans <- 0.5 * dnorm(x, mean = mu1, sd = 0.1) / k0 + 0.5 * dnorm(x, mean = mu2, sd = 0.2) / k1
  ans[x>1]=0
  ans[x<0]=0
  return(ans)
}
curve(f,0,1)
curve(fhat,0,1, add=T, col=2)
## compute ISE
tem <- function(u){ (fhat(u)-f(u))^2 }
integrate(tem, 0, 1)

#(a) The sample size n increases to 5000 or 10000.
###ISE = 0.003107452 (n = 5000)
###ISE = 0.002577162 (n = 10000)

#(b) The knots are replaced with 1/9, 2/9, . . ., 8/9.
###ISE = 0.01607938 (n = 1000)
###ISE = 0.001435022 (n = 5000)


#(c) The knots are replaced with the knots in Part (b) and the sample size n increases to 10000.
###ISE = 0.001098023 (n = 10000)



##########Exercise 2.
library(splines)
set.seed(1)
mu1 <- 0.2
mu2 <- 0.7
n <- 10000
knotlist <- (1:4) / 5
nb <- length(knotlist) + 4  

m <- n * 10
z <- rnorm(m, mean = mu1, sd = 0.1)
x <- z[(z > 0) & (z < 1)][1:n]
z <- rnorm(m, mean = mu2, sd = 0.2)
x2 <- z[(z > 0) & (z < 1)][1:n]
z <- sample(0:1, size = n, replace = TRUE)
x[z == 1] <- x2[z == 1]


bx <- bs(x, knots = knotlist, Boundary.knots = c(0, 1), intercept = TRUE)

log_likelihood <- function(a) {
  fa_denominator <- integrate(function(u) {exp(bs(u, knots = knotlist, Boundary.knots = c(0,1), intercept = TRUE) %*% a)}
    , 0, 1)$value
  
  lambda_a <- log(fa_denominator)
  log_fa <- bx %*% (a - lambda_a)
  return(-sum(log_fa))  
}

initial <- rep(0, nb)
mle_result <- optim(initial, log_likelihood)
a_hat <- mle_result$par

f_hat <- function(u) {
  exp(bs(u, knots = knotlist, Boundary.knots = c(0, 1), intercept = TRUE) %*% a_hat) /
    integrate(function(u) exp(bs(u, knots = knotlist, Boundary.knots = c(0, 1), intercept = TRUE) %*% a_hat), 0, 1)$value
}


k0 <- pnorm(1, mean = mu1, sd = 0.1) - pnorm(0, mean = mu1, sd = 0.1)
k1 <- pnorm(1, mean = mu2, sd = 0.2) - pnorm(0, mean = mu2, sd = 0.2)

f <- function(x) {
  ans <- 0.5 * dnorm(x, mean = mu1, sd = 0.1) / k0 + 0.5 * dnorm(x, mean = mu2, sd = 0.2) / k1
  ans[x > 1] <- 0
  ans[x < 0] <- 0
  return(ans)
}
ISE <- integrate(function(u) (f_hat(u) - f(u))^2, 0, 1)$value
ISE
###ISE = 0.004256832 (n = 10000) knots are 1/9, 2/9, . . ., 8/9

###ISE = 0.001375642 (n = 10000) knots are 1/5, 2/5, . . ., 4/5
###ISE = 0.005324884 (n = 1000) knots are 1/5, 2/5, . . ., 4/5
curve(f, 0, 1, col = "black")
curve(f_hat, 0, 1, col = "red", add = TRUE)








##########Exercise 3
set.seed(1)
n <- 1000
m <- n * 10
z <- rnorm(m, mean = mu1, sd = 0.1)
x <- z[(z > 0) & (z < 1)][1:n]
z <- rnorm(m, mean = mu2, sd = 0.2)
x2 <- z[(z > 0) & (z < 1)][1:n]
z <- sample(0:1, size = n, replace = TRUE)
x[z == 1] <- x2[z == 1]


looml_estimate <- function(x, knotlist) {

nb <- length(knotlist) + 4
M <- matrix(0, nb, nb)

for (i in 1:nb){
  for (j in i:nb){
    tem <- function(u){
      bx <- bs(u, knots = knotlist, Boundary.knots = c(0,1), intercept=T)
      return( bx[,i]*bx[,j])
    }
    M[i,j] <- integrate(tem, 0, 1)$value
    if (j > i) { M[j,i] <- M[i,j] }
  }
}

moments <- apply(bs(x, knots = knotlist, Boundary.knots=c(0,1), intercept=T), 2, mean)
ahat <- solve( M, moments)
fhat <- function(u){
  ans <- bs(u, knots = knotlist, Boundary.knots = c(0,1), intercept=T) %*% ahat
  return( as.numeric(ans) )
}

log_ml<- 0

for (i in 1:length(x)) {
  x_looml <- x[-i] 
  moments_looml <- apply(bs(x_looml, knots = knotlist, Boundary.knots = c(0, 1), intercept = TRUE), 2, mean)
  ahat_looml <- solve(M, moments_looml)
  
  fhat_looml <- function(u) {
    ans <- bs(u, knots = knotlist, Boundary.knots = c(0, 1), intercept = TRUE) %*% ahat_looml
    return(as.numeric(ans))
  }
  log_ml <- log_ml + log(fhat_looml(x[i]))  
}
return(log_ml)
}





likelihood_cv_result <- rep(0, 30)
for (i in 1:30) {
    
   n_knots_4 <- (1:4) / 5
   n_knots_8 <- (1:8) / 9
    
   looml_estimate
    
    
    if (log_lik_4 > log_lik_8) {
      results[sim] <- 4
    } else {
      results[sim] <- 8
    }
}
  