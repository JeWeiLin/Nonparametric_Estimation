###########Kernel regression
###########Exercise 1.
set.seed(12345)
# gaussian_kernel <- function(x) {
#   return(1 / sqrt(2 * pi ) * exp(-x^2 / 2))
# }
kernel_estimator <- function(x, y, x_eval, bandwidth) {
  
  x_diff <- x_eval - x
  
  kernel <- dnorm(x_diff /bandwidth)
  
  density <- sum(y * kernel ) / sum(kernel)
  return(density)
}

n <- 1000
x <- rnorm(n)
y <- rnorm(n)

bandwidth <- 0.1

x_eval <- 0.1

kernel_estimate <- kernel_estimator(x, y, x_eval, bandwidth)
kernel_estimate



####plot
x_value <- seq(-2.5, 2.5, length=100)

y_0 <- sin(x_value) + rnorm(100, sd =0.1)

plot(x_value, y_0, xlab = "x", ylab = "y")

kernel_estimate <- NULL
x_eval <- seq(-2.5, 2.5, length=100)

for(i in 1:100){
  kernel_estimate[i] <- kernel_estimator(x_value, y_0 , x_eval[i], bandwidth)
}

plot(x_value, y_0, xlab = "x", ylab = "y")
lines(x_value, sin(x_value) )
lines(x_value, kernel_estimate, xlab = "x", ylab = "y", col = "red",lty=4)



###########Exercise 3.
n <- 1000
x <- runif(n)
y <- rnorm(n,sd =0.1) + sin(6*x)

loocv <- function(x, y, bandwidths) {
  rss <- NULL
  
  for (i in seq_along(bandwidths)) {
    h <- bandwidths[i]
    
    m_hat <- NULL
    
    for (j in seq_along(x)) {
      x_train <- x[-j]
      x_test <- x[j]
      
      m_hat <-c(m_hat, kernel_estimator(x[-j], y[-j] , x_test,  h))
    }
    
    rss <- c (rss, sum((y - m_hat) ^ 2))
  }
  index_h <- which.min(rss)
  
  return(list(best_bandwidth = bandwidths[index_h], 
              RSSCV = kernel_estimator(x, y, x_eval, bandwidths[index_h])))
}

h <- c(0.01, 0.05, 0.1, 0.4, 0.002, 0.025, 0.2)

loocv(x, y, bandwidths = h )

#####

loocv_2 <- function(x, y, bandwidths) {
  rss <- NULL
  
  for (i in seq_along(bandwidths)) {
    h <- bandwidths[i]
    
    m_hat <- NULL
    
    for (j in seq_along(x)) {
      
      m_hat <-c(m_hat, kernel_estimator(x, y , x[j],  h))
    }
    
    rss <- c (rss, sum((y - m_hat) ^ 2))
  }
  index_h <- which.min(rss)
  
  return(list(best_bandwidth = bandwidths[index_h], 
              RSSCV = kernel_estimator(x, y, x_eval, bandwidths[index_h])))
}

h <- c(0.01, 0.05, 0.1, 0.4, 0.002, 0.025, 0.2)

loocv_2(x, y, bandwidths = h )
