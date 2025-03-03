###########B-splines
###########Exercise 2.
library(splines)

m <- 4 
k <- 3  
knotlist <- (1:3) / 4  
x <- (1:1000) / 1001 

bx <- bs(x, knots = knotlist, degree = m - 1, intercept = TRUE)
base_func_terms <- cbind(1, x, x^2, x^3) 

xi_terms <- sapply(knotlist, function(xi) {(x - xi)^3 * (x > xi)})

terms <- cbind(base_func_terms, xi_terms)

lm_list <- lapply(1:ncol(bx), function(i) {lm(bx[, i] ~ terms - 1)})

for(i in 1:7){
  print(summary(lm_list[[i]])$r.squared)
}


###########Exercise 3.(a)
library(splines)

f <- function(x) x * sin(20 * x)

n <- 1000
x <- seq(0, 1, length = n)
y <- f(x)

ise_f <- function(y, y_hat, x) {
  ise <- function(x) {
    (y(x) - y_hat(x))^2
  }
  result <- integrate(ise, 0, 1)
  return(result$value)
}


bspline_f <- function(k) {
  knotlist <- seq(1, k) / (k + 1)
  bx <- bs(x, knots = knotlist, degree = 3, intercept = TRUE)
  fit <- lm(y ~ bx - 1)
  
  y_pred <- function(x) {
    basis_func <- predict(bx, x)
    basis_func %*% coef(fit)
  }
  ise <- ise_f(f, y_pred, x)
  return(ise)
}


k_values <- 1:7
ise_values <- sapply(k_values, bspline_f)
ise_values
best_k <- k_values[which.min(ise_values)]
best_k

###########Exercise 3.(b)
ise_poly <- function(m, x, y, f) {
  
  x_poly_term <- sapply(0:m, function(p) x^p)  
  
  fit <- lm(y ~ x_poly_term - 1) 
  
  f_hat <- function(xi) {
    x_poly_hat <- sapply(0:m, function(p) xi^p)
    x_poly_hat %*% coef(fit)  
  }
  curve(f_hat)
  ise <- integrate(function(x) (f(x) - f_hat(x))^2, 0, 1)$value
  
  return(ise)
}

m_degree <- 15 
ise_poly_value <- rep(0, m_degree)

for (m in 1:m_degree) {
  ise_poly_value[m] <- ise_poly(m, x, y, f)
  
  if (ise_values_poly[m] < min(ise_values)) {
    m0 <-m
    break  
  }
}
m0
