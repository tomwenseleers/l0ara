# # pure R implementation of l0ara, allowing for any family & link function & with nonnegativity constraints
# # also based on IRLS algo (using expected Hessian) instead of Newton-Raphson (using observed Hessian) as in l0araC,
# # with initialization of GLM as in R's GLM code
# 
# library(pracma)
# pinv.solve = function (X, y) pinv(X) %*% y # Moore-Penrose / SVD-based solve for p>n case
# 

#' @export
l0araR = function(X, y, weights=rep(1,nrow(X)), family, lam=0, nonnegative=FALSE, maxit=30, tol=1E-8, thresh=1E-3) {
  if (!is(family, "family")) family = family()
  variance = family$variance
  linkinv = family$linkinv
  mu.eta = family$mu.eta
  etastart = NULL

  nobs = nrow(X)    # nobs & nvars are needed by the initialize expression below
  nvars = ncol(X)
  eval(family$initialize) # initializes n and mustart
  eta = family$linkfun(mustart) # we initialize eta based on this
  dev.resids = family$dev.resids
  dev = sum(dev.resids(y, linkinv(eta), weights))
  devold = 0

  iter = 1
  # ybar = mean(y) # this was initialization used by l0ara, but I use initialization of GLMs as in R's glm code
  # beta = rep(0, nvars)
  # beta[1] = family$linkfun(ybar) # assumes 1st column of X is intercept term

  Xt = X
  if (nobs > nvars) { M = .sparseDiagonal(nvars,1)
                      # diag(1,nvars)
  } else { M = .sparseDiagonal(nobs,1)
            # M = diag(1,nobs)
    }
  beta = rep(1,nvars) # make sparse instead?
  while (iter < maxit) {
    mu = linkinv(eta)
    varg = variance(mu)
    gprime = mu.eta(eta)
    z = eta + (y - mu) / gprime # potentially -offset
    W = weights * as.vector(gprime^2 / varg)
    nonzero = (beta!=0)
    beta_old = beta

    # TODO: take into account zero & converged coefficients, add block fit, add coordinate descent option, allow sparse X
    # add nonnegativity or box constraints also via nnls or nnnpls or osqp instead of via clipping only

    # lm_adridge = function (X, y, weights, start=rep(1, ncol(X))) # TODO put function below in separate function lm_adridge
    # cf Liu & Li 2016 & Liu et al 2017 & https://github.com/tomwenseleers/l0ara/blob/master/src/l0araC.cpp

    if ((nobs > nvars)&lam>0) {
      # beta[nonzero] = solve(Matrix::t(Xt[,nonzero,drop=F]) %*% (X[,nonzero,drop=F] * W) + diag(lam,sum(nonzero)),
      #                      Matrix::t(Xt[,nonzero,drop=F]) %*% (z * W), tol=2*.Machine$double.eps)
      # PS: faster than row augmentation method:
      # beta = matrix(qr.solve(rbind(X*sqrt(W),sqrt(lam)*diag(1,ncol(X))),c(y*sqrt(W),rep(0,ncol(X)))),ncol=1)
      # or
      beta[nonzero] = solve_sparse_lsconjgrad( list(rbind(X*sqrt(weights), .sparseDiagonal(ncol(X), x=sqrt(lam))), c(y*sqrt(weights), rep(0,ncol(X))), rep(0,ncol(X))) )$res
    }
    if ((nobs <= nvars)&lam>0) { # fastest option also for n=p case
      # beta[nonzero] = Matrix::t(Xt[,nonzero,drop=F]) %*% solve((X[,nonzero,drop=F] * W) %*% Matrix::t(Xt[,nonzero,drop=F]) + lam*M,
      #                                                 z * W, tol=2*.Machine$double.eps) # 25 ms for n=p=500
      # faster for large matrices: 
      beta[nonzero] = Matrix::t(Xt[,nonzero,drop=F]) %*% solve_sparse_lsconjgrad(list((X[,nonzero,drop=F] * W) %*% Matrix::t(Xt[,nonzero,drop=F]) + lam*M,z * W,rep(0,nrow(X))))$res
      # or try with row-augment matrix?
    }
    # as possible alternatives for solve:
    # RcppEigen::fastLmPure(mm, y, 2L) # LDLt Cholesky decomposition with rank detection
    # for very large problems, conjugate gradient / iterative solvers might be worth checking,
    # eg cPCG package (pcgsolve), cf https://github.com/styvon/cPCG
    # library(cPCG)
    # microbenchmark(pcgsolve((X[,nonzero,drop=F] * W) %*% t(Xt[,nonzero,drop=F]) + lam*M, z * W,"ICC")) # 100 ms for n=p=500
    # cPCG:::pcgsolve_sparseOMP((X[,nonzero,drop=F] * W) %*% t(Xt[,nonzero,drop=F]) + lam*M, z * W, "ICC", nThreads=2)

    # Rlinsolve package, cf
    # microbenchmark(lsolve.bicgstab((X[,nonzero,drop=F] * W) %*% t(Xt[,nonzero,drop=F]) + lam*M,
    #        z * W, xinit=NA, reltol = 1e-05, maxiter = 1000,
    #        preconditioner = diag(ncol(A)), verbose = FALSE), times=10) # 293 ms for n=p=500

    if (lam==0) beta[nonzero] = solve(crossprod(X[,nonzero,drop=F],W*X[,nonzero,drop=F]),
                                      crossprod(X[,nonzero,drop=F],W*z), tol=2*.Machine$double.eps) # use pinv.solve() here for p>=n case or when lam is too small above?

    if (thresh!=0) beta[abs(beta) < thresh] = 0

    if (nonnegative) beta[beta<0] = 0 # clipping to enforce nonnegativity

    eta = X %*% beta # potentially +offset
    dev = sum(dev.resids(y, mu, weights))

    if (lam>0) Xt = matrix(beta * beta, nrow=nobs, ncol=nvars, byrow=TRUE) * X
    iter = iter + 1

    if (iter >= maxit) { warning("Algorithm did not converge. Increase maxit.") }
    # print(sum(nonzero))
    if (sum(nonzero)==0) break
    if (norm(beta-beta_old, "2") < tol) break
    # if (abs(dev - devold) / (0.1 + abs(dev)) < tol) break # iterations converge when |dev - dev_{old}|/(|dev| + 0.1) < tol

    devold = dev
  }

  if (thresh!=0) beta[abs(beta) < thresh] = 0

  return(list(coefficients = beta, iter = iter))
}
# 
# 
# # TEST ON SPIKE TRAIN CONVOLUTED WITH GAUSSIAN PSF AND WITH POISSON NOISE
# library(L0glm)
# library(Matrix)
# # Simulate some data
# n <- 500
# p <- 500
# s <- 0.1 # sparsity as proportion of p that have nonzero coefficients
# k <- round(p*s) # nr of nonzero covariates
# sim <- simulate_spike_train(n = n, p = p, k = k,
#                             mean_beta = 1000, sd_logbeta = 1,
#                             family = "poisson", seed = 123, Plot = TRUE)
# X <- sim$X
# y <- sim$y
# weights <- rep(1,nrow(X))
# family <- poisson(identity)
# lam <- 1
# nonnegative <- TRUE
# maxit <- 100
# tol <- 1E-8
# weighted.rmse <- function(actual, predicted, weight){
#   sqrt(sum((predicted-actual)^2*weight)/sum(weight))
# }
# weightedrmse_betas <- function(betas) apply(betas, 2, function (fitted_coefs) {
#   weighted.rmse(actual=sim$y_true,
#                 predicted=sim$X %*% fitted_coefs,
#                 weight=1/sim$y_true) } )
# library(microbenchmark)
# library(l0ara)
# library(L0glm)
# library(nnls)
# microbenchmark("wnnls" = { wnnls_fit <- nnls(A=sim$X*sqrt(1/(sim$y+1)), b=sim$y*sqrt(1/(sim$y+1))) },
#                "l0araR_pois" = { l0araR_pois_fit <- l0araR(X=sim$X, y=sim$y, family=poisson(identity), 
#                                                            lam=1, nonnegative=TRUE, maxit=25, tol=1E-8) },
#                "l0araR_pois_wnnls_prescreen" = { wnnls_fit <- nnls(A=sim$X*sqrt(1/(sim$y+1)), b=sim$y*sqrt(1/(sim$y+1)))
#                wnnls <- wnnls_fit$x
#                l0araR_pois_wnnls_coefs <- rep(0,ncol(sim$X))
#                l0araR_pois_wnnls_coefs[wnnls!=0] <- coef(l0araR(X=sim$X[,wnnls!=0],
#                                                                 y=sim$y,
#                                                                 family=poisson(identity),
#                                                                 lam=1, nonnegative=TRUE, maxit=25, tol=1E-8)) },
#                "l0araR_wgaus" = { l0araR_wgaus_fit <- l0araR(X=sim$X, 
#                                                              y=sim$y, 
#                                                              weights=1/(sim$y+1),
#                                                              family=gaussian(identity), 
#                                                              lam=1, nonnegative=TRUE, maxit=25, tol=1E-8) },
#                "l0araR_wgaus_wnnls_prescreen" = { wnnls_fit <- nnls(A=sim$X*sqrt(1/(sim$y+1)), b=sim$y*sqrt(1/(sim$y+1)))
#                wnnls <- wnnls_fit$x
#                l0araR_wgaus_wnnls_coefs <- rep(0,ncol(sim$X))
#                l0araR_wgaus_wnnls_coefs[wnnls!=0] <- coef(l0araR(X=sim$X[,wnnls!=0],
#                                                                  y=sim$y,
#                                                                  family=gaussian(identity),
#                                                                  lam=1, nonnegative=TRUE, maxit=25, tol=1E-8)) },
#                "l0ara_wgaus" = { l0ara_wgaus_fit <- l0ara(x=sim$X*sqrt(1/(sim$y+1)), y=sim$y*sqrt(1/(sim$y+1)), family="gaussian", 
#                                                           lam=1, standardize=FALSE, maxit=25, eps=1E-8)},
#                "L0glm_pois" = { L0glm_pois_fit <- L0glm.fit(X=sim$X, y=sim$y,
#                                                             family = poisson(identity),
#                                                             lambda = 1, nonnegative = TRUE, normalize = FALSE,
#                                                             control.l0 = list(maxit = 25, rel.tol = 1e-04, 
#                                                                               delta = 1e-05, gamma = 2, warn = FALSE), 
#                                                             control.iwls = list(maxit = 1, rel.tol = 1e-04, thresh = 1e-03, warn = FALSE), 
#                                                             control.fit = list(maxit = 1, block.size = NULL, tol = 1e-07)) },
#                "L0glm_wgaus" = { L0glm_wgaus_fit <- L0glm.fit(X=sim$X, y=sim$y,
#                                                               weights=1/(sim$y+1),
#                                                               family = gaussian(identity),
#                                                               lambda = 1, nonnegative = TRUE, normalize = FALSE,
#                                                               control.l0 = list(maxit = 25, rel.tol = 1e-04, 
#                                                                                 delta = 1e-05, gamma = 2, warn = FALSE), 
#                                                               control.iwls = list(maxit = 1, rel.tol = 1e-04, thresh = 1e-03, warn = FALSE), 
#                                                               control.fit = list(maxit = 1, block.size = NULL, tol = 1e-07)) }, 
#                times=3)
# 
# 
# # Unit: milliseconds
# #   expr                             min        lq      mean    median        uq       max neval
# # wnnls                         142.4643  142.4643  142.4643  142.4643  142.4643  142.4643     1
# # l0araR_pois                   648.3939  648.3939  648.3939  648.3939  648.3939  648.3939     1
# # l0araR_pois_wnnls_prescreen   303.8300  303.8300  303.8300  303.8300  303.8300  303.8300     1
# # l0araR_wgaus                  750.4784  750.4784  750.4784  750.4784  750.4784  750.4784     1
# # l0araR_wgaus_wnnls_prescreen  346.9583  346.9583  346.9583  346.9583  346.9583  346.9583     1
# # l0ara_wgaus                  1726.4469 1726.4469 1726.4469 1726.4469 1726.4469 1726.4469     1
# # L0glm_pois                   1303.9442 1303.9442 1303.9442 1303.9442 1303.9442 1303.9442     1
# # L0glm_wgaus                  1523.9141 1523.9141 1523.9141 1523.9141 1523.9141 1523.9141     1
# 
# l0araR_pois_fit # 518 ms  / 536 ms     # 612 ms, 26 iters for n=p=500; 18 ms, 33 iters for n=p=100; 9 ms, 18 iters for n=100, p=101
# l0araR_wgaus_fit # 612 ms, 26 iters for n=p=500; 18 ms, 33 iters for n=p=100; 9 ms, 18 iters for n=100, p=101 
# # 2 iterations - not correct...
# l0ara_wgaus_fit # 1.7 s, 52 iters for n=p=500; 21 ms, 55 iters for n=p=100; 10 ms, 22 iters for n=100, p=101
# L0glm_pois_fit # 1.3 s, 22 iters for n=p=500; 18 ms, 21 iters for n=p=100; 21 ms, 15 iters for n=100, p=101
# L0glm_wgaus_fit # conv after 10 iters
# 
# # some benchmarks
# betas = data.frame(beta_true=sim$beta_true,
#                    wnnls=wnnls_fit$x, 
#                    l0araR_pois=coef(l0araR_pois_fit),
#                    l0araR_pois_wnnls_prescreen=l0araR_pois_wnnls_coefs,
#                    l0araR_wgaus=coef(l0araR_wgaus_fit),
#                    l0araR_wgaus_wnnls_prescreen=l0araR_wgaus_wnnls_coefs,
#                    l0ara_wgaus=coef(l0ara_wgaus_fit),
#                    L0glm_pois=coef(L0glm_pois_fit),
#                    L0glm_wgaus=coef(L0glm_wgaus_fit))
# betas[betas<0] = 0
# TP=colSums(betas[sim$beta_true!=0,]>0) # TPs
# FP=colSums(betas[sim$beta_true==0,]>0) # FPs
# FN=colSums(betas[sim$beta_true>0,]==0) # FNs
# TN=colSums(betas[sim$beta_true==0,]==0) # TNs
# ACC = (TP+TN)/p # ACCURACY
# SENS = TP/(TP + FN) # sensitivity
# SPEC = TN/(TN + FP) # specificity
# RELABSBIAS = colMeans(100*abs(betas[sim$beta_true!=0,]-sim$beta_true[sim$beta_true!=0])/sim$beta_true[sim$beta_true!=0])
# ABSBIAS = colMeans(abs(betas[sim$beta_true!=0,]-sim$beta_true[sim$beta_true!=0]))
# 
# data.frame(wrmse=weightedrmse_betas(betas), TP=TP, FP=FP, FN=FN, TN=TN, 
#            ACC=ACC, SENS=SENS, SPEC=SPEC, RELABSBIAS=RELABSBIAS, ABSBIAS=ABSBIAS)
# #                                     wrmse TP FP FN  TN   ACC SENS      SPEC RELABSBIAS  ABSBIAS
# # beta_true                    0.000000e+00 50  0  0 450 1.000 1.00 1.0000000    0.00000   0.0000
# # wnnls                        6.279938e-09 47 48  3 402 0.898 0.94 0.8933333   21.53539 278.5526
# # l0araR_pois                  5.465149e-09 44  9  6 441 0.970 0.88 0.9800000   16.74723 154.0157
# # l0araR_pois_wnnls_prescreen  5.265800e-09 45  9  5 441 0.972 0.90 0.9800000   14.89734 185.5639
# # l0araR_wgaus # BEST          5.188420e-09 46  7  4 443 0.978 0.92 0.9844444   12.81361 139.3991
# # l0araR_wgaus_wnnls_prescreen 7.000543e-09 47 30  3 420 0.934 0.94 0.9333333   21.98483 314.8148
# # l0ara_wgaus                  5.242743e-09 45  9  5 441 0.972 0.90 0.9800000   15.00487 143.7398
# # L0glm_pois                   5.737391e-09 45  9  5 441 0.972 0.90 0.9800000   15.78284 153.2482
# # L0glm_wgaus                  5.555031e-09 44  8  6 442 0.972 0.88 0.9822222   16.70547 144.0974
# 
# 
# # l0araR_pois  0.03892734 18   3  2 177 0.975 0.90 0.9833333   15.78212  156.5168
# # l0araR_pois  0.03892734 18   3  2 177 0.975 0.90 0.9833333   15.78212  156.5168
# 
# # plot solution quality
# par(mfrow=c(2,1))
# plot(x = sim$x, y = sim$y, ylim=c(-max(sim$y),max(sim$y)), type="l", ylab="Signal", xlab="Time", main="Real spike train (red) & l0araR estimated spike train (blue)")
# lines(x = sim$x, y = -sim$y)
# lines(x = sim$x_beta[sim$beta_true>0], y=sim$beta_true[sim$beta_true>0], col="red", type="h")     
# lines(x = sim$x_beta[coef(l0araR_fit)>0], y=-coef(l0araR_fit)[coef(l0araR_fit)>0], col="blue", type="h")     
# plot(x = sim$x, y = sim$y, ylim=c(-max(sim$y),max(sim$y)), type="l", ylab="Signal", xlab="Time", main="Real spike train (red) & l0ara estimated spike train (blue)")
# lines(x = sim$x, y = -sim$y)
# lines(x = sim$x_beta[sim$beta_true>0], y=sim$beta_true[sim$beta_true>0], col="red", type="h")     
# lines(x = sim$x_beta[coef(l0ara_fit)>0], y=-coef(l0ara_fit)[coef(l0ara_fit)>0], col="blue", type="h")     
# 
# 
# 
# # CHECK TO SEE THAT WITH LAMBDA=0 WE GET THE REGULAR GLM SOLUTION
# ## Dobson (1990) Page 93: Randomized Controlled Trial :
# y <- counts <- c(18,17,15,20,10,20,25,13,12)
# outcome <- gl(3,1,9)
# treatment <- gl(3,3)
# X <- model.matrix(counts ~ outcome + treatment)
# print(data <- data.frame(treatment, outcome, counts))
# library(microbenchmark)
# microbenchmark("glm" = { glmfit <- glm.fit(x=X, y=y, family = poisson(log), control=glm.control(epsilon = 1e-8, maxit = 25)) },
#                "l0araR" = { L0araRfit <- l0araR(X=X, y=y, family=poisson(log),
#                                                 lam=0, nonnegative=FALSE, maxit=25, tol=1E-8, thresh=0) }, times=500) # comparable timings
# # PS iterations converge when |dev - dev_{old}|/(|dev| + 0.1) < ε
# 
# cbind(coef(glmfit),coef(L0araRfit))
# 
# # glm.D93.ngl <- glm.cons(counts ~ outcome + treatment, family = poisson(log),
# #    method="glm.fit.cons")
# # summary(glm.D93)
# # summary(glm.D93.ngl)
# 
# 
# 
# # DEBUGGING / TEST AREA
# 
# # lmridge_solve = function (X, Y, lambda) {
# #   nobs = nrow(X)
# #   nvars = ncol(X)
# #   if (nobs > nvars) {
# #     return(solve(crossprod(X) + diag(lambda, ncol(X)), crossprod(X, Y))[, 1]) # cf Liu & Li 2016 algo 1, see Liu et al 2017 & https://github.com/tomwenseleers/l0ara/blob/master/src/l0araC.cpp for application to GLMs
# #     l
# #   } else {
# #     return(t(X) %*% solve(tcrossprod(X) + diag(lambda, nrow(X)), Y)[, 1])
# #   } # cf Liu & Li 2016 algo 2, see Liu et al 2017 & https://github.com/tomwenseleers/l0ara/blob/master/src/l0araC.cpp for application to GLMs
# # }
# # 
# # lmridge_qrsolve = function (X, Y, lambda) {
# #   nobs = nrow(X)
# #   nvars = ncol(X)
# #   if (nobs > nvars) {
# #     return(qr.solve(crossprod(X) + diag(lambda, ncol(X)), crossprod(X, Y))[, 1]) # cf Liu & Li 2016 algo 1, see Liu et al 2017 & https://github.com/tomwenseleers/l0ara/blob/master/src/l0araC.cpp for application to GLMs
# #     l
# #   } else {
# #     return(t(X) %*% qr.solve(tcrossprod(X) + diag(lambda, nrow(X)), Y)[, 1])
# #   } # cf Liu & Li 2016 algo 2, see Liu et al 2017 & https://github.com/tomwenseleers/l0ara/blob/master/src/l0araC.cpp for application to GLMs
# # }
# # 
# chol.solve = function (a, b) {
#   ch <- chol(crossprod(a)) # solve using Cholesky decomposition
#   backsolve(ch, forwardsolve(ch, crossprod(a, b), upper = TRUE, trans = TRUE))
# }
# # 
# # lmridge_cholsolve = function (X, Y, lambda) {
# #   nobs = nrow(X)
# #   nvars = ncol(X)
# #   if (nobs > nvars) {
# #     return(chol.solve(crossprod(X) + diag(lambda, ncol(X)), crossprod(X, Y))[, 1]) # cf Liu & Li 2016 algo 1, see Liu et al 2017 & https://github.com/tomwenseleers/l0ara/blob/master/src/l0araC.cpp for application to GLMs
# #     l
# #   } else {
# #     return(t(X) %*% chol.solve(tcrossprod(X) + diag(lambda, nrow(X)), Y)[, 1])
# #   } # cf Liu & Li 2016 algo 2, see Liu et al 2017 & https://github.com/tomwenseleers/l0ara/blob/master/src/l0araC.cpp for application to GLMs
# # }
# # 
# # 
# # # R equivalent of repmat (matlab)
# # repmat = function (X,m,n) {
# #   mx = dim(X)[1]
# #   nx = dim(X)[2]
# #   matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
# # }  
# # 

