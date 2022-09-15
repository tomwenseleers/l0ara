# TEST ON SPIKE TRAIN CONVOLUTED WITH GAUSSIAN PSF AND WITH POISSON NOISE
# remotes::install_github("zdebruine/RcppML", ref="main")
library(RcppML) # for fast nnls function
library(L0glm)
library(Matrix)
# Simulate some data
n <- 10000 # dimensions of chromatograms is 10000 timepoints x ca 550 mass channels
p <- n
s <- 0.01 # sparsity as proportion of p that have nonzero coefficients
k <- round(p*s) # nr of nonzero covariates
maxit = 500
sim <- simulate_spike_train(n = n, npeaks = k, 
                            peakhrange = c(10, 1000),
                            seed = 123, Plot = TRUE)
X <- sim$X # class(X) = "matrix" = coded as dense covariate matrix here
y <- sim$y
range(sim$a) # real coefficients 0.0000 936.0128
# approx 1/variance Poisson weights if identity link Poisson GLM is approximated by weighted least squares regression
weights = 1/(y+0.1) 
family <- poisson(identity)
lam <- 1E-5
nonnegative <- TRUE
maxit <- 100
tol <- 1E-8
weighted.rmse <- function(actual, predicted, weight){
  sqrt(sum((predicted-actual)^2*weight)/sum(weight))
}
weightedrmse_betas <- function(betas) apply(betas, 2, function (fitted_coefs) {
  weighted.rmse(actual=sim$y_true,
                predicted=sim$X %*% fitted_coefs,
                weight=1/sim$y_true) } )
library(microbenchmark)
library(l0ara)
library(L0glm)
library(nnls)


# for weighted least square regression we need to multiply X and y by sqrt(weights)
Xw = Matrix(X*sqrt(weights))
yw = y*sqrt(weights)
# class(Xw)
# dgeMatrix

# coded as sparse Matrix
X_sparse = Matrix(sim$X, sparse=TRUE)
Xw_sparse = Matrix(sim$X, sparse=TRUE)*sqrt(weights)
# class(Xw_sparse)
# dgC.Matrix


system.time(coefs <- wls_solve(Xw_sparse, yw, weights)$res) # 1.3s
range(coefs/sqrt(weights))
range(sim$a) # 0.0000 936.0128
range(coefs) # -4.848032 10.847179
plot(sim$a, coefs, type="p", pch=16) # coefficients not correctly scaled yet

lam = 1E-5
X_aug_sparse = rbind(X, .sparseDiagonal(ncol(X), x=sqrt(lam)))
weights_aug = c(weights, rep(1,p))
Xw_aug_sparse = X_aug_sparse*weights_aug
y_aug = c(y, rep(0,ncol(X)))
yw = y*sqrt(weights)
yw_aug = c(yw, rep(0,ncol(X)))

system.time(coefs <- wls_solve(X_aug_sparse, y_aug, weights_aug)$res) # 1.26s
plot(sim$a, coefs, type="p", pch=16) # coefs not correctly scaled
system.time(coefs <- solve_sparse_lsconjgrad( list(Xw_aug_sparse, 
                                          yw_aug, 
                                          rep(0,p)) )$res) # 4s, or as.numeric(y) as initialisation?
plot(sim$a, coefs, type="p", pch=16) # looks OK
system.time(coefs <- RcppML::nnls(a=as.matrix(crossprod(Xw)), 
                                  b=as.matrix(crossprod(Xw, yw)),
                                  fast_nnls = FALSE)) # 4s
plot(sim$a, coefs, type="p", pch=16) # looks OK


system.time(l0arafit <- l0ara(x=Xw, 
                              y=yw, 
      family="gaussian", lam=1, nonneg=TRUE, 
      standardize=FALSE, maxit=5, eps=1E-5)) 
# 22s with X either dense or sparse
plot(sim$a, coef(l0arafit), pch=16)
# plot solution quality
par(mfrow=c(1,1))
plot(x = sim$x, y = sim$y, ylim=c(-max(sim$y),max(sim$y)), type="l", ylab="Signal", xlab="Time", main="Real spike train (red) & l0ara estimated spike train (blue)")
lines(x = sim$x, y = -sim$y)
lines(x = sim$x[sim$a>0], y=sim$a[sim$a>0], col="red", type="h")
lines(x = sim$x[coef(l0arafit)>0], y=-coef(l0arafit)[coef(l0arafit)>0], col="blue", type="h")


# coding X as a sparse matrix (class dgCMatrix or dsCMatrix)
# using pure R implementation here
system.time(l0araRfit <- l0ara:::l0araR(X=Matrix(sim$X,sparse=TRUE), 
                                y=sim$y, weights=weights,
                 family=gaussian(identity), lam=1, 
                 nonnegative=TRUE, maxit=maxit, tol=1E-5)) # 5s
plot(sim$a, coef(l0araRfit), pch=16)
# plot solution quality
par(mfrow=c(1,1))
plot(x = sim$x, y = sim$y, ylim=c(-max(sim$y),max(sim$y)), type="l", ylab="Signal", xlab="Time", main="Real spike train (red) & l0araR estimated spike train (blue)")
lines(x = sim$x, y = -sim$y)
lines(x = sim$x_beta[sim$a>0], y=sim$a[sim$a>0], col="red", type="h")
lines(x = sim$x_beta[coef(l0araRfit)>0], y=-coef(l0araRfit)[coef(l0araRfit)>0], col="blue", type="h")



library(L0adri1.0)
system.time(l0adrifit <- l0adri(X=Matrix(sim$X,sparse=TRUE), y=sim$y,
                   weights=weights, 
                   family=gaussian(identity), 
                   lam=1, nonnegative=rep(TRUE,ncol(sim$X)), algo="clipping", miniter=1, maxit=50, tol=1E-6)) # 1s
plot(sim$a, l0adrifit$beta, pch=16)
# plot solution quality
par(mfrow=c(1,1))
plot(x = sim$x, y = sim$y, ylim=c(-max(sim$y),max(sim$y)), type="l", ylab="Signal", xlab="Time", main="Real spike train (red) & l0ara estimated spike train (blue)")
lines(x = sim$x, y = -sim$y)
lines(x = sim$x_beta[sim$a>0], y=sim$a[sim$a>0], col="red", type="h")
lines(x = sim$x_beta[coef(l0arafit)>0], y=-coef(l0arafit)[coef(l0arafit)>0], col="blue", type="h")


system.time(wnnls_fit <- nnls(A=sim$X*sqrt(1/(sim$y+1)), b=sim$y*sqrt(1/(sim$y+1)))) # 0.13s











microbenchmark("wnnls" = { wnnls_fit <- nnls(A=sim$X*sqrt(1/(sim$y+1)), b=sim$y*sqrt(1/(sim$y+1))) },
               "l0araR_pois" = { l0araR_pois_fit <- l0araR(X=sim$X, y=sim$y, family=poisson(identity),
                                                           lam=1, nonnegative=TRUE, maxit=25, tol=1E-8) },
               "l0araR_pois_wnnls_prescreen" = { wnnls_fit <- nnls(A=sim$X*sqrt(1/(sim$y+1)), b=sim$y*sqrt(1/(sim$y+1)))
               wnnls <- wnnls_fit$x
               l0araR_pois_wnnls_coefs <- rep(0,ncol(sim$X))
               l0araR_pois_wnnls_coefs[wnnls!=0] <- coef(l0araR(X=sim$X[,wnnls!=0],
                                                                y=sim$y,
                                                                family=poisson(identity),
                                                                lam=1, nonnegative=TRUE, maxit=25, tol=1E-8)) },
               "l0araR_wgaus" = { l0araR_wgaus_fit <- l0araR(X=sim$X,
                                                             y=sim$y,
                                                             weights=1/(sim$y+1),
                                                             family=gaussian(identity),
                                                             lam=1, nonnegative=TRUE, maxit=25, tol=1E-8) },
               "l0araR_wgaus_wnnls_prescreen" = { wnnls_fit <- nnls(A=sim$X*sqrt(1/(sim$y+1)), b=sim$y*sqrt(1/(sim$y+1)))
               wnnls <- wnnls_fit$x
               l0araR_wgaus_wnnls_coefs <- rep(0,ncol(sim$X))
               l0araR_wgaus_wnnls_coefs[wnnls!=0] <- coef(l0araR(X=sim$X[,wnnls!=0],
                                                                 y=sim$y,
                                                                 family=gaussian(identity),
                                                                 lam=1, nonnegative=TRUE, maxit=25, tol=1E-8)) },
               "l0ara_wgaus" = { l0ara_wgaus_fit <- l0ara(x=sim$X*sqrt(1/(sim$y+1)), y=sim$y*sqrt(1/(sim$y+1)), family="gaussian",
                                                          lam=1, standardize=FALSE, maxit=25, eps=1E-8)},
               "L0glm_pois" = { L0glm_pois_fit <- L0glm.fit(X=sim$X, y=sim$y,
                                                            family = poisson(identity),
                                                            lambda = 1, nonnegative = TRUE, normalize = FALSE,
                                                            control.l0 = list(maxit = 25, rel.tol = 1e-04,
                                                                              delta = 1e-05, gamma = 2, warn = FALSE),
                                                            control.iwls = list(maxit = 1, rel.tol = 1e-04, thresh = 1e-03, warn = FALSE),
                                                            control.fit = list(maxit = 1, block.size = NULL, tol = 1e-07)) },
               "L0glm_wgaus" = { L0glm_wgaus_fit <- L0glm.fit(X=sim$X, y=sim$y,
                                                              weights=1/(sim$y+1),
                                                              family = gaussian(identity),
                                                              lambda = 1, nonnegative = TRUE, normalize = FALSE,
                                                              control.l0 = list(maxit = 25, rel.tol = 1e-04,
                                                                                delta = 1e-05, gamma = 2, warn = FALSE),
                                                              control.iwls = list(maxit = 1, rel.tol = 1e-04, thresh = 1e-03, warn = FALSE),
                                                              control.fit = list(maxit = 1, block.size = NULL, tol = 1e-07)) },
               times=3)


# Unit: milliseconds
#   expr                             min        lq      mean    median        uq       max neval
# wnnls                         142.4643  142.4643  142.4643  142.4643  142.4643  142.4643     1
# l0araR_pois                   648.3939  648.3939  648.3939  648.3939  648.3939  648.3939     1
# l0araR_pois_wnnls_prescreen   303.8300  303.8300  303.8300  303.8300  303.8300  303.8300     1
# l0araR_wgaus                  750.4784  750.4784  750.4784  750.4784  750.4784  750.4784     1
# l0araR_wgaus_wnnls_prescreen  346.9583  346.9583  346.9583  346.9583  346.9583  346.9583     1
# l0ara_wgaus                  1726.4469 1726.4469 1726.4469 1726.4469 1726.4469 1726.4469     1
# L0glm_pois                   1303.9442 1303.9442 1303.9442 1303.9442 1303.9442 1303.9442     1
# L0glm_wgaus                  1523.9141 1523.9141 1523.9141 1523.9141 1523.9141 1523.9141     1

l0araR_pois_fit # 518 ms  / 536 ms     # 612 ms, 26 iters for n=p=500; 18 ms, 33 iters for n=p=100; 9 ms, 18 iters for n=100, p=101
l0araR_wgaus_fit # 612 ms, 26 iters for n=p=500; 18 ms, 33 iters for n=p=100; 9 ms, 18 iters for n=100, p=101
# 2 iterations - not correct...
l0ara_wgaus_fit # 1.7 s, 52 iters for n=p=500; 21 ms, 55 iters for n=p=100; 10 ms, 22 iters for n=100, p=101
L0glm_pois_fit # 1.3 s, 22 iters for n=p=500; 18 ms, 21 iters for n=p=100; 21 ms, 15 iters for n=100, p=101
L0glm_wgaus_fit # conv after 10 iters

# some benchmarks
betas = data.frame(a=sim$a,
                   wnnls=wnnls_fit$x,
                   l0araR_pois=coef(l0araR_pois_fit),
                   l0araR_pois_wnnls_prescreen=l0araR_pois_wnnls_coefs,
                   l0araR_wgaus=coef(l0araR_wgaus_fit),
                   l0araR_wgaus_wnnls_prescreen=l0araR_wgaus_wnnls_coefs,
                   l0ara_wgaus=coef(l0ara_wgaus_fit),
                   L0glm_pois=coef(L0glm_pois_fit),
                   L0glm_wgaus=coef(L0glm_wgaus_fit))
betas[betas<0] = 0
TP=colSums(betas[sim$a!=0,]>0) # TPs
FP=colSums(betas[sim$a==0,]>0) # FPs
FN=colSums(betas[sim$a>0,]==0) # FNs
TN=colSums(betas[sim$a==0,]==0) # TNs
ACC = (TP+TN)/p # ACCURACY
SENS = TP/(TP + FN) # sensitivity
SPEC = TN/(TN + FP) # specificity
RELABSBIAS = colMeans(100*abs(betas[sim$a!=0,]-sim$a[sim$a!=0])/sim$a[sim$a!=0])
ABSBIAS = colMeans(abs(betas[sim$a!=0,]-sim$a[sim$a!=0]))

data.frame(wrmse=weightedrmse_betas(betas), TP=TP, FP=FP, FN=FN, TN=TN,
           ACC=ACC, SENS=SENS, SPEC=SPEC, RELABSBIAS=RELABSBIAS, ABSBIAS=ABSBIAS)
#                                     wrmse TP FP FN  TN   ACC SENS      SPEC RELABSBIAS  ABSBIAS
# a                    0.000000e+00 50  0  0 450 1.000 1.00 1.0000000    0.00000   0.0000
# wnnls                        6.279938e-09 47 48  3 402 0.898 0.94 0.8933333   21.53539 278.5526
# l0araR_pois                  5.465149e-09 44  9  6 441 0.970 0.88 0.9800000   16.74723 154.0157
# l0araR_pois_wnnls_prescreen  5.265800e-09 45  9  5 441 0.972 0.90 0.9800000   14.89734 185.5639
# l0araR_wgaus # BEST          5.188420e-09 46  7  4 443 0.978 0.92 0.9844444   12.81361 139.3991
# l0araR_wgaus_wnnls_prescreen 7.000543e-09 47 30  3 420 0.934 0.94 0.9333333   21.98483 314.8148
# l0ara_wgaus                  5.242743e-09 45  9  5 441 0.972 0.90 0.9800000   15.00487 143.7398
# L0glm_pois                   5.737391e-09 45  9  5 441 0.972 0.90 0.9800000   15.78284 153.2482
# L0glm_wgaus                  5.555031e-09 44  8  6 442 0.972 0.88 0.9822222   16.70547 144.0974


# l0araR_pois  0.03892734 18   3  2 177 0.975 0.90 0.9833333   15.78212  156.5168
# l0araR_pois  0.03892734 18   3  2 177 0.975 0.90 0.9833333   15.78212  156.5168

# plot solution quality
par(mfrow=c(2,1))
plot(x = sim$x, y = sim$y, ylim=c(-max(sim$y),max(sim$y)), type="l", ylab="Signal", xlab="Time", main="Real spike train (red) & l0araR estimated spike train (blue)")
lines(x = sim$x, y = -sim$y)
lines(x = sim$x_beta[sim$a>0], y=sim$a[sim$a>0], col="red", type="h")
lines(x = sim$x_beta[coef(l0araR_fit)>0], y=-coef(l0araR_fit)[coef(l0araR_fit)>0], col="blue", type="h")
plot(x = sim$x, y = sim$y, ylim=c(-max(sim$y),max(sim$y)), type="l", ylab="Signal", xlab="Time", main="Real spike train (red) & l0ara estimated spike train (blue)")
lines(x = sim$x, y = -sim$y)
lines(x = sim$x_beta[sim$a>0], y=sim$a[sim$a>0], col="red", type="h")
lines(x = sim$x_beta[coef(l0ara_fit)>0], y=-coef(l0ara_fit)[coef(l0ara_fit)>0], col="blue", type="h")
