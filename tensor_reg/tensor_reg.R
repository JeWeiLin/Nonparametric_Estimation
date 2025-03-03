############Regression using tensor basis functions
############Exercise 1.
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
  coef.hat <- lm(data.y~ bx.data - 1)$coef
  print(length(coef.hat))
  fhat <- function(x){
    bx.x <-  bx.tensor(x, bx.uni)
    return( as.numeric( bx.x %*% coef.hat )  )
  }
  return(fhat)
}



set.seed(1)
f <- function(x){ dnorm(x[,1]-0.5, sd=0.2)*dnorm(x[,2]-0.5, sd=0.2) }
n <- 1000
X <- matrix(runif(n*2), n,2)
y <- f(X)  + rnorm(n,sd=0.4)
fhat <- get.fhat(X, y, bspline)



n.mc <- 10000
u <- matrix( runif(n.mc*2), n.mc, 2)
mean( (fhat(u) - f(u))^2 )





###########
############Exercise 2.
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





##############################
############Exercise 3.
set.seed(1)
f <- function(x){
  dnorm(x[,1]-0.5, sd=0.2) + dnorm(x[,2]-0.5, sd=0.2) + dnorm(x[,3]-0.5, sd=0.2)
}
n <- 100000
X <- matrix(runif(n*3), n,3)
y <- f(X) + rnorm(n,sd=0.4)

###addictive(ISE)
fhat_result <- additive(y, X, 3)
n.mc <- 10000
u <- matrix( runif(n.mc*3), n.mc, 3)
mean( (fhat_result(u) - f(u))^2 )



###compare ISE
tensor_result <- get.fhat(X, y, bspline)
u <- matrix( runif(n.mc*3), n.mc, 3)
mean( (tensor_result(u) - f(u))^2 )
