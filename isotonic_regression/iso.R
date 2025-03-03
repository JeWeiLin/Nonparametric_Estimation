##########Constrained curve fitting based on splines
##########Exercise 1.
library(quadprog)
dbs <- function(x, knotlist, bknots, m=4, der=0){
  J <- m+length(knotlist)
  n <- length(x)
  knots.all <- c(rep(bknots[1], m), knotlist, rep(bknots[2], m))
  dbx <- matrix(0, n, J)
  for (j in 1:J){
    k <- knots.all[j:(j+m)]
    dbx[,j] <- splineDesign(k, x, ord=m, derivs=rep(der,n), outer.ok=TRUE)
  }
  return(dbx)
}

get.fit <- function(x, y, knotlist, deg = 2, bknots = c(0, 1)) {
  
  B <- bs(x, knots = knotlist, degree = deg, Boundary.knots = bknots, intercept = TRUE)
  
  Dmat <- t(B) %*% B
  dvec <- as.numeric(y %*% B)
  
  xi <- c(bknots[1], knotlist, bknots[2])
 
  Amat <- -t(dbs(xi,  knotlist,  bknots, m = deg + 1, der=1))
  bvec <- rep(0, length(xi)) 
  
  a.hat <- solve.QP(Dmat, dvec, Amat, bvec)$solution
  
  fhat <- function(u) {
    bu <- bs(u, knots = knotlist, degree = deg, Boundary.knots = bknots, intercept = TRUE)
    ans <- bu %*% a.hat
    return(ans[,1])
  }
  return(fhat)
}


set.seed(1)
n <- 1000
x <- seq(0, 1, length=n)
f <- function(x){ ans <- exp(-x); ans[x>=0.5] <- exp(-0.5); return(ans) }
y <- f(x)+rnorm(n, sd=0.05)


fhat <- get.fit(x, y, (1:7) / 8 )
plot(x,y)
curve(fhat, 0, 1, add=TRUE, col=2)
curve(f, 0,1, add=TRUE, col=3)

#compute ISE
g <- function(x){ (fhat(x)-f(x))^2 }
integrate(g,0,1)$value




##########Exercise 2.
set.seed(1)
n <- 1000
x <- seq(-3, 3, length = n)
f <- function(x) { dnorm(x, 0, 1) }
y <- f(x) + rnorm(n, sd = 0.05)

knots <-  (1:7) * 6 / 8 - 3
deg <- 2
B <- bs(x, knots = knots, degree = deg, Boundary.knots = c(-3, 3), intercept = TRUE)

Dmat <- t(B) %*% B
dvec <- t(y) %*% B

xi <- c(-3, knots, 3)  

der_bspline_matrix <- t(dbs(xi, knotlist = knots, bknots = c(-3, 3), m = deg + 1, der = 1))

left_constraints <- der_bspline_matrix[, xi <= 0]    # increasing on [-3, 0]
right_constraints <- -der_bspline_matrix[, xi >= 0]  # decreasing on [0, 3]

Amat <- cbind(left_constraints, right_constraints)

bvec <- rep(0, ncol(Amat))

a_hat <- solve.QP(Dmat, dvec, Amat, bvec)$solution

fhat <- function(u) {
  bu <- bs(u, knots = knots, degree = deg, Boundary.knots = c(-3, 3), intercept = TRUE)
  as.numeric(bu %*% a_hat)
}

ISE <- integrate(function(u) (f(u) - fhat(u))^2, lower = -3, upper = 3)$value
ISE

plot(x, y)
curve(fhat, -3, 3, add=TRUE, col=2, lwd =2)
curve(f, -3, 3, add=TRUE, col=3, lwd =2)

