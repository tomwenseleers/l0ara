glm_irls = function(X, y, weights=rep(1,nrow(X)), family=poisson(log), maxit=25, tol=1e-16) {
    if (!is(family, "family")) family = family()
    variance = family$variance
    linkinv = family$linkinv
    mu.eta = family$mu.eta
    etastart = NULL
    
    nobs = nrow(X)    # needed by the initialize expression below
    nvars = ncol(X)   # needed by the initialize expression below
    eval(family$initialize) # initializes n and fitted values mustart
    eta = family$linkfun(mustart) # we then initialize eta with this
    dev.resids = family$dev.resids
    dev = sum(dev.resids(y, linkinv(eta), weights))
    devold = 0
    beta_old = rep(1, nvars)
    
    for(j in 1:maxit)
    {
      mu = linkinv(eta) 
      varg = variance(mu)
      gprime = mu.eta(eta)
      z = eta + (y - mu) / gprime # potentially -offset if you would have an offset argument as well
      W = weights * as.vector(gprime^2 / varg)
      beta = solve(crossprod(X,W*X), crossprod(X,W*z), tol=2*.Machine$double.eps)
      eta = X %*% beta # potentially +offset if you would have an offset argument as well
      dev = sum(dev.resids(y, mu, weights))
      if (abs(dev - devold) / (0.1 + abs(dev)) < tol) break
      # if (norm(beta-beta_old, "2") < tol) break
      devold = dev
      beta_old = beta
    }
    list(coefficients=beta, iterations=j)
}




# ## Dobson (1990) Page 93: Randomized Controlled Trial :
# y <- counts <- c(18,17,15,20,10,20,25,13,12)
# outcome <- gl(3,1,9)
# treatment <- gl(3,3)
# X <- model.matrix(counts ~ outcome + treatment)
# print(d.AD <- data.frame(treatment, outcome, counts))
# library(microbenchmark)
# microbenchmark(glmfit.D93 <- glm.fit(x=X, y=y, family = poisson(log))) # 540 µs
# # glm.D93.ngl <- glm.cons(counts ~ outcome + treatment, family = poisson(log),
# #    method="glm.fit.cons")
# # summary(glm.D93)
# # summary(glm.D93.ngl)
# coef(glm.D93)
# microbenchmark(glm_irls(X=X, y=y, family=poisson(log), maxit=10, tol=1E-10)) # 558 µs
# coef(glm_irls(X=X, y=y, family=poisson(log), maxit=10, tol=1E-10))
# 
# microbenchmark(l0araR(X=X, y=y, family=poisson(log),
#                       lam=1E-10, maxit=10, tol=1E-5)$coefficients) # 532 µs with solve 
# fit = l0araR(X=X, y=y, family=poisson(log), lam=1E-12, maxit=30, tol=1E-5)
# fit
# fit$coefficients


