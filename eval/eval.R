###########Evalution of a nonparametric function estimator
###########Exercise 1.
N <- 100   
n <- 1000 
sigma <- 0.05  
h <- c(0.005, 0.01, 0.1) 

x_eval <- seq(0, 1, 1000)  

m_hat <- function(x, y, x_eval, h) {
  m_hat <- rep(0, length(x_eval))
  
  for (i in seq_along(x_eval)) {
    x_diff <- x_eval[i] - x
    kernel <- dnorm(x_diff / h)
    m_hat[i] <- sum(y * kernel) / sum(kernel)
  }
  return(m_hat)
}

m_hat_1 <- function(x_0){
  m_hat(x, y, x_0, h)
}
g <- Vectorize(m_hat_1)


ise <- rep(0, N)

imse_f <- function(h){
  for (i in 1:N) {
    x <- runif(n)
    epsilon <- rnorm(n, 0, sigma)
    y <- sin(20 * x) + epsilon
    
    m_hat_1 <- function(x_0){
      m_hat(x, y, x_0, h)
    }
    g <- Vectorize(m_hat_1)
    
    integrand <- function(u) {
      (g(u) - sin(20 * u)) ^ 2
    }
    ise[i] <- integrate(integrand, 0, 1)$value
  }
  return(mean(ise))
}
imse_f(0.005)
# 0.0003374082
imse_f(0.01)
# 0.0007617356
imse_f(0.1)
# 0.3460142

###########Exercise 2.(a)

n <- 1000
sigma <- 0.05
h <- c(0.1, 0.01)

kernel_estimator <- function(x, y, x_eval, h) {
  
  x_diff <- x_eval - x
  
  kernel <- dnorm(x_diff /h)
  
  density <- sum(y * kernel ) / sum(kernel)
  return(density)
}

rsscv_f <- function(h) {
  rsscv_values <- rep(0, 10)
  
  for (i in 1:10) {

    x <- runif(n)  
    epsilon <- rnorm(n, 0, sigma)
    y <- sin(20 * x) + epsilon  
    
    m_hat <-rep(0, 10)
    
    for (j in seq_along(x)) {
      x_train <- x[-j]
      x_test <- x[j]
      y_train <- y[-j]
      
      m_hat[j] <- kernel_estimator(x_train, y_train, x_test,  h)
    }

    rsscv <- sum((y - m_hat) ^2)  
    rsscv_values[i] <- rsscv / 1000 - sigma^2  
  }
  return(rsscv_values)
}

rsscv_results <- sapply(h, rsscv_f)
rsscv_results


###########Exercise 2.(b)

rsscv_f_beta <- function(h) {
  rsscv_values <- rep(0, 10)
  
  for (i in 1:10) {
    x <- rbeta(n, 2, 2)  
    epsilon <- rnorm(n, 0, sigma)
    y <- sin(20 * x) + epsilon
    
    m_hat <- rep(0, n)
    
    for (j in seq_along(x)) {
      x_train <- x[-j]
      x_test <- x[j]
      y_train <- y[-j]
      
      m_hat[j] <- kernel_estimator(x_train, y_train, x_test, h)
    }
    
    rsscv <- sum((y - m_hat) ^ 2)
    rsscv_values[i] <- rsscv / n - sigma^2
  }
  
  return(rsscv_values)
}
rsscv_results_beta <- sapply(h, rsscv_f_beta)
rsscv_results_beta


###########Exercise 2.(c)

ise <- rep(0, N)

imse_f_beta <- function(h) {
  for (i in 1:N) {
    x <- rbeta(n, 2, 2)  
    epsilon <- rnorm(n, 0, sigma)
    y <- sin(20 * x) + epsilon
    
    m_hat_1 <- function(x_0) {
      m_hat(x, y, x_0, h)
    }
    g <- Vectorize(m_hat_1)
    
    integrand <- function(u) {
      (g(u) - sin(20 * u)) ^ 2 * dbeta(u, 2, 2)  
    }
    
    ise[i] <- integrate(integrand, 0, 1)$value
  }
  return(mean(ise))
}
imse_beta_estimate <- imse_f_beta(h = 0.1)



