###########B-splines
###########Exercise 4.(a)
f0 <- function(x){
  ans <- x*sin(20*x)
  ans[x<0] <- 0
  return(ans)
}

f <- function(x){ f0(2*(x-0.5))}
curve(f,0,1)


n <- 1000
x <- seq(0, 1, length=n)
y <- f(x) + rnorm(n, sd=0.02)

knotlist <- (1:13) /14

base_func_terms <- cbind(1, x, x^2, x^3)

xi_terms <- sapply(knotlist, function(xi) {(x - xi)^3 * (x > xi)})

terms <- cbind(base_func_terms, xi_terms)

colnames(terms) <- 1:17

fit <- lm(y ~ terms - 1)

pvalue_threshold <- 0.05  
terms_pvalues <- summary(fit)$coefficients[, 4]  


while (any(terms_pvalues[-(1:4)] > pvalue_threshold)) {  

  insignificant_knots_index <- which.max(terms_pvalues[-(1:4)]) + 4  
  
  terms <- terms[, -insignificant_knots_index]
  
  fit <- lm(y ~ terms - 1)
  
  terms_pvalues <- summary(fit)$coefficients[, 4]
}

summary(fit)


###########Exercise 4.(b)
f0 <- function(x){
  ans <- x*sin(20*x)
  ans[x<0] <- 0
  return(ans)
}

f <- function(x){ f0(2*(x-0.5))}


knotlist <-  (6:13) / 14
n <- 1000
x <- seq(0, 1,length=n)
y <- f(x)

bx <- bs(x, knots = knotlist, degree = 3, intercept = TRUE)
f_fit <- lm(y ~ bx - 1)

y_pred <- function(x) {
  bx <- bs(x, knots = knotlist, degree = 3, intercept = TRUE)
  return((bx %*% f_fit$coefficients)[,1])
}

ise_f <- function(f, f_hat) {
  result <- integrate(ise, 0, 1)
  return(result$value)
}
ise <- function(x) {
  (f(x) - y_pred(x))^2
}
integrate(ise, 0, 0.99)
curve(ise)
curve(y_pred)
curve(f, add=T)
ise_f(f, y_pred)

ise_poly(12, x, y, f)
curve()