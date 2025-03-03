##########Goodness of fit tests involving nonparametric function estimation
##########Exercise 1.
#####(a) chi-square approximation of the distribution of nlog(W) under H0
library(splines)
n <- 150
set.seed(1)
x <- runif(n)
y <- sin(2 * x) + runif(n, -0.1, 0.1) * 10


spline_term <- bs(x, knots = 0.5, degree = 3) 
spline_fit <- lm(y ~ spline_term)
lm_fit <- lm(y ~ x)

lm_residual <- summary(lm_fit)$resid
spline_residual <- summary(spline_fit)$resid

lm_rss <- sum(lm_residual ^ 2)
spline_rss <- sum(spline_residual ^ 2)

W_statistic <-  lm_rss / spline_rss 
test_statistic <- n * log(W_statistic)

df_difference <- df.residual(lm_fit) - df.residual(spline_fit)
p_value <- 1 - pchisq(test_statistic, df = df_difference)
p_value
###[1] 0.02878809


#####(b) bootstrap with 200 bootstrap trials.
set.seed(2) 
bootstrap <- 200
bootstrap_reesult <- numeric(bootstrap)

for (i in 1:bootstrap) {
 
  bootstrap_sample <- fitted(lm_fit) +  sample(residuals(lm_fit), replace = TRUE)
  
  bootstrap_lm_fit <- lm(bootstrap_sample ~ x)
  bootstrap_spline_fit <- lm(bootstrap_sample ~ spline_term)
  
  bootstrap_lm_residual <- summary(bootstrap_lm_fit)$resid
  bootstrap_spline_residual <- summary(bootstrap_spline_fit)$resid
  
  bootstrap_lm_rss <- sum(bootstrap_lm_residual ^ 2)
  bootstrap_spline_rss <- sum(bootstrap_spline_residual ^ 2)
  
  bootstrap_W_statistic <-  bootstrap_lm_rss / bootstrap_spline_rss 
  bootstrap_test_statistic <- n * log(bootstrap_W_statistic)

  bootstrap_reesult[i] <- bootstrap_test_statistic > test_statistic
}

bootstrap_p_value <- mean(bootstrap_reesult)
bootstrap_p_value
###[1] 0.015



##########Exercise 2.
p_value_result <- rep(0, 500)
for(i in 1:500){
  set.seed(i)
  n <- 5000
  x <- runif(n)
  y <-  1 +x+runif(n, -0.1, 0.1)*5
  
  spline_term <- bs(x, knots = 0.5, degree = 3) 
  spline_fit <- lm(y ~ spline_term)
  lm_fit <- lm(y ~ x)
  
  lm_residual <- summary(lm_fit)$resid
  spline_residual <- summary(spline_fit)$resid
  
  lm_rss <- sum(lm_residual ^ 2)
  spline_rss <- sum(spline_residual ^ 2)
  
  W_statistic <-  lm_rss / spline_rss 
  test_statistic <- n * log(W_statistic)
  
  df_difference <- df.residual(lm_fit) - df.residual(spline_fit)
  p_value_result[i] <- 1 - pchisq(test_statistic, df = df_difference)

}
hist(p_value_result)
ks.test(p_value_result, punif)


###Do Part (a) again with n <- 150 replaced by n <- 5000.
# Asymptotic one-sample Kolmogorov-Smirnov test
# data:  p_value_result
# D = 0.023597, p-value = 0.9435
# alternative hypothesis: two-sided


##########Exercise 3.
set.seed(1) 
n <- 1000
x <- matrix(runif(n*2), n, 2)
y <- 1 + x[,1] * sin(3 * x[,2]) + rnorm(n, sd=0.1)


##########tensor basis  
bspline <- function(x, k = 3) {
  bs_result <- bs(x, degree = 3, knots = seq(1 / (k + 1), k / (k + 1), length.out = k), intercept = TRUE)
  return(bs_result)
}

bx.tensor <- function(x, bx.uni) {
  if (is.null(dim(x))) {
    return(bx.uni(x))
  }
  n <- dim(x)[1]
  d <- dim(x)[2]
  mat.list <- vector("list", d)
  for (i in 1:d) {
    mat.list[[i]] <- bx.uni(x[, i])
  }
  n.basis1 <- dim(mat.list[[1]])[2]###單個維度中基底函數的數量
  n.basisd <- n.basis1 ^ d###每個維度的基底函數數量的乘積
  v <- vector("list", d)
  for (i in 1:d) {
    v[[i]] <- 1:n.basis1
  }
  ind.mat <- as.matrix(expand.grid(v))
  bx <- matrix(0, n, n.basisd)
  for (j in 1:n.basisd) {
    ind <- ind.mat[j, ]
    b.prod <- rep(1, n)
    for (k in 1:d) {
      b.prod <- b.prod * mat.list[[k]][, ind[k]]
    }
    bx[, j] <- b.prod
  }
  return(bx)
}

get.fhat <- function(data.x, data.y, bx.uni){
  bx.data <- bx.tensor(data.x, bx.uni)
  coef.hat <- lm(data.y ~ bx.data - 1)$coef
  
  fhat <- function(x){
    bx.x <-  bx.tensor(x, bx.uni)
    return( as.numeric( bx.x %*% coef.hat )  )
  }
  return(fhat)
}

fhat_tensor <- get.fhat(x, y, bspline)
tensor_residual <- sum((y - fhat_tensor(x)) ^ 2)


##########additive
additive <- function(data.y, data.x, k) {
  n <- nrow(data.x)  
  d <- ncol(data.x)  
  
  additive_matrix <- bs(data.x[,1], degree = 3, knots = seq(1/(k+1), k/(k+1), length.out = k), intercept = TRUE)
  
  for (i in 2:d) {
    basis <- bs(data.x[,i], degree = 3, knots = seq(1/(k+1), k/(k+1), length.out = k), intercept = TRUE)[,-1]
    
    additive_matrix <- cbind(additive_matrix, basis)
  }
  
  lm_y <- lm(data.y ~ additive_matrix - 1)  
  
  #return(summary(lm_y))
  fhat <- function(u){
    additive_matrix <- bs(u[,1], degree = 3, knots = seq(1/(k+1), k/(k+1), length.out = k), intercept = TRUE)
    
    for (i in 2:d) {
      basis <- bs(u[,i], degree = 3, knots = seq(1/(k+1), k/(k+1), length.out = k), intercept = TRUE)[,-1]
      
      additive_matrix <- cbind(additive_matrix, basis)
    }
    return(as.numeric(additive_matrix %*% lm_y$coef))
  }
  return(fhat)
}
fhat_additive <- additive(y, x, 3)
additive_residual <- sum((y - fhat_additive(x)) ^ 2)

additive_residual / n 
tensor_residual / n 
W_statistic <-  additive_residual / tensor_residual 
test_statistic <- n * log(W_statistic)

# (a) set.seed(1) n <- 1000
# x <- matrix(runif(n*2), n,2)
# y <-  1+x[,1]+sin(3*x[,2])+rnorm(n, sd=0.1)

# (b) set.seed(1) n <- 1000
# x <- matrix(runif(n*2), n,2)
# y <-  1+x[,1]*sin(3*x[,2])+rnorm(n, sd=0.1)

# 在H0(可加性模型)成立下，W_statistic 相較於相乘模型(即虛無假設不成立情況下)，
# 檢定統計量(W_statistic)會較小，即此檢定可以區別兩模型的差別。


# set.seed(1) 
# n <- 1000
# x <- matrix(runif(n*2), n, 2)
# y <- 1 + x[,1] + sin(3 * x[,2]) + rnorm(n, sd=0.1)
# 
# fhat_additive <- additive(y, x, 3)
# additive_residual <- sum((y - fhat_additive(x)) ^ 2)
# 
# W_statistic <-  additive_residual / tensor_residual 
# test_statistic <- n * log(W_statistic)
# test_statistic



n <- 1000
x <- matrix(runif(n*2), n, 2)
y <- 1 + x[,1] + sin(3 * x[,2]) + rnorm(n, sd=0.1)

fhat_additive <- additive(y, x, 3)
fhat_tensor <- get.fhat(x, y, bspline)

tensor_residual <- sum((y - fhat_tensor(x)) ^ 2)
additive_residual <- sum((y - fhat_additive(x)) ^ 2)

W_statistic <-  additive_residual / tensor_residual 
test_statistic <- n * log(W_statistic)



bootstrap_compute_h0 <- function( bootstrap = 20) {
  
bootstrap_statistic_reesult <- rep(0, bootstrap)
  
for (i in 1:bootstrap) {
  bootstrap_sample <- fhat_additive(x) +  sample(((y - fhat_additive(x)) ^ 2), replace = TRUE)
  
  additive_result <- additive(bootstrap_sample, x, 3)
  additive_residual <- sum((bootstrap_sample - additive_result(x)) ^ 2)
  
  tensor_result <- get.fhat(x, bootstrap_sample, bspline)
  tensor_residual <- sum((bootstrap_sample - tensor_result(x)) ^ 2)
  
  W_statistic <- additive_residual / tensor_residual
  bootstrap_statistic_reesult[i] <- n * log(W_statistic)
}
  return(  1 - (sum(test_statistic > bootstrap_statistic_reesult ) / bootstrap)) 

}

bootstrap_statistic_reesult_h0 <- rep(0, 100)

for (i in 1:100 ) {
  bootstrap_statistic_reesult_h0[i] <- bootstrap_compute_h0()
}
sum(bootstrap_statistic_reesult_h0 > 0.05)
#######################


n <- 1000
x <- matrix(runif(n*2), n, 2)
y <- 1 + x[,1] * sin(3 * x[,2]) + rnorm(n, sd=0.1)

fhat_additive <- additive(y, x, 3)
fhat_tensor <- get.fhat(x, y, bspline)

tensor_residual <- sum((y - fhat_tensor(x)) ^ 2)
additive_residual <- sum((y - fhat_additive(x)) ^ 2)

W_statistic <-  additive_residual / tensor_residual 
test_statistic <- n * log(W_statistic)


bootstrap_compute_h1 <- function( bootstrap = 20) {
  
  bootstrap_statistic_result <- rep(0, bootstrap)
  
  for (i in 1:bootstrap) {
    bootstrap_sample <- fhat_additive(x) +  sample(((y - fhat_additive(x)) ^ 2), replace = TRUE)
    
    additive_result <- additive(bootstrap_sample, x, 3)
    additive_residual <- sum((bootstrap_sample - additive_result(x)) ^ 2)
    
    tensor_result <- get.fhat(x, bootstrap_sample, bspline)
    tensor_residual <- sum((bootstrap_sample - tensor_result(x)) ^ 2)
    
    W_statistic <- additive_residual / tensor_residual
    bootstrap_statistic_result[i] <- n * log(W_statistic)
  }
  return(  1 - (sum(test_statistic > bootstrap_statistic_result ) / bootstrap)) 
  
}

bootstrap_statistic_reesult_h1 <- rep(0, 100)

for (i in 1:100 ) {
  bootstrap_statistic_reesult_h1[i] <- bootstrap_compute_h1()
}

