###########Kernel regression
###########Exercise 2.(a)
m <- function(x) {
  sin(20 * x)
}



E_I <- function(x0, h) {
  integrand <- function(u) {
    (m(x0 - h * u) - m(x0)) * dnorm(u)
    
  }
  integrate(integrand, lower = -Inf, upper = Inf)$value
}

x0 <- 0.1
h <- c(0.01, 0.005, 0.001, 0.0005)
E_I_results <- sapply(h, function(h) E_I(x0, h))


E_I_results


###########Exercise 2.(b)
n <- 1000
N <- 10000

f <- function(x) {
  if (x >= -1 && x <= 1) {
    return(1 / 2)
  } else {
    return(0)
  }
}

I_x0 <- function(x0, X, Y, h) {
  (1 / (n * h)) * sum((Y - m(x0)) * k((x0 - X) / h) * sapply(x0 - (X / h), f_uniform))
}

x0 <- 0.1
h_values <- c(0.01, 0.0005)


results <- rep(0, length(h_values))
I_results <- matrix(0, nrow = N, ncol = length(h_values))

for (i in 1:N) {
  
  X <- runif(n, -1, 1)
  epsilon <- rnorm(n, mean = 0, sd = 0.01)
  Y <- m(X) + epsilon
  
  I_results[i, ] <- sapply(h_values, function(h) I_x0(x0, X, Y, h))
}
E_I_approx <- colMeans(I_results)

