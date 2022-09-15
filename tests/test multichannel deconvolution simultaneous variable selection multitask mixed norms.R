# tests simultaneous variable selection using mixed norms

library(l0ara)

# simulate data
# GC-MS data is simulated as a linear combination of spectrum profiles S and
# elution profiles C.

n <- p <- 100 # take as many covariate p as the number of observations n
x <- 1:n
npeaks <- 10 # nr of peaks = nr of spectra
set.seed(1)
peakhrange.log <- c(2, 5) # peak height range
w1 = 1.25 # simulated peak widths bigaussian peak shape
w2 = 2.5

# Function to generate an asymmetrical Gaussian peak
bigauspeak = function(x, u, w1, w2, h=1) {
  out <- rep(0, length(x))
  out[x<=u] <- h*exp(-(x[x<=u]-u)^2/(2*w1^2))
  out[x>u] <- h*exp(-(x[x>u]-u)^2/(2*w2^2))
  return(out)
}
# Banded matrix with theoretical peak shape function used, this will serve as
# the covariate matrix from which we sample the true set of variables for selection
# TODO you can still add background covariates in this one as well (e.g. B spline basis)
X <- sapply(1:p, function (u) bigauspeak(x, u = u, w1 = w1, w2 = w2, h = 1))
# Generate matrix C containing the true set of covariates
u <- sample(x, npeaks, replace=FALSE) # simulated peak locations
# h <- 10^runif(npeaks, min = min(peakhrange), max = max(peakhrange)) # simulated peak heights (in total ion current)
h <- 10^seq(2, 4, length.out = npeaks)
C <- sweep(X[,u], 2, h, "*")
par(mfrow = c(2,1))
matplot(C+1, log = "y", type = "l")
image(C^0.2, col = topo.colors(255))

# Generate the spectrum profiles
# Load real-life profiles
file <- system.file("Simulated spectrum.rds", package = "l0ara") # TODO test this when building
file <- "./inst/Simulated spectrum.rds"
S <- readRDS(file)
m <- ncol(S)
# Randomly select some profiles
set.seed(1)
nonzero <- sample(1:nrow(S), npeaks) 
S <- S[nonzero,]

# Simulate response variable
Y0 <- C %*% S # noiseless simulated signal = linear convolution of spike train with peak shape function
Y <- apply(Y0 + 1, 2, function (col) {
  out <- rpois(n, col) # simulated signal with random poisson noise on it
  out[is.na(out)] <- rnorm(n = sum(is.na(out)), mean = col[is.na(out)],
                           sd = sqrt(col[is.na(out)])) # to deal with overflows
  return(out)
})

setwd("~/GitHub/l0ara/tests")
write.csv(X,"X.csv",row.names=F)
write.csv(Y,"Y.csv",row.names=F)
write.csv(Y0,"Y0.csv",row.names=F)
write.csv(C,"C.csv",row.names=F)
write.csv(S,"S.csv",row.names=F)


# Plot the data
const <- 1
par(mfrow = c(1,1))
matplot(C+const, type="l", log="y", lwd=3, lty=1, col="red")
matlines(Y+const, lty = 1, col = rgb(0,0,0,0.2))
abline(v=u, col="blue", lwd=log10(h)/3)


# Case I: Fit a multichannel glm (no penalty)
library(nnls)
system.time(fit1 <- group.pen.nnglm(X = X, Y = Y, family = poisson(identity),
                                    # Weights = matrix(1, nrow = nrow(Y), ncol = ncol(Y)),
                                    intercept = TRUE,
                                    group.pen.fun = function(mat) sqrt(rowSums(mat^2)), # default = L2 norm of the estimated coefficients per component
                                    lambdas = 0,
                                    maxit = 10))
par(mfrow = c(3,1))
image(y = 0:m, x = 0:n, Y0^0.2, col = topo.colors(255), main = "Ground truth")
abline(v=u, col="green")
image(y = 0:m, x = 0:n, Y^0.2, col = topo.colors(255), main = "Input data")
abline(v=u, col="green")
image(y = 0:m, x = 0:n, fit1$fitted.values^0.2, col = topo.colors(255), main = "Fitted data")
abline(v=u, col="green")


# Case II: Fit a multichannel L0 penalized glm
system.time(fit2 <- group.pen.nnglm(X = X, Y = Y, family = poisson(identity),
                                    # Weights = matrix(1, nrow = nrow(Y), ncol = ncol(Y)),
                                    intercept = TRUE,
                                    lambdas = rep(1E-2, m),
                                    # lambdas = 1E-5*(1+apply(Y, 2, max))^0.001, # lambda should be function of scaling of each mass channel, we specify small lambda because we initialise with nnls solution
                                    
                                    group.pen.fun = function(mat){
                                      out <- sqrt(rowSums(mat^2)) # row L2 norms are used for adaptive ridge penalty
                                      out <- out/max(out) # L infinity norm normalized L2 row norms to make them range between 0 and 1 - @CHRIS - maybe figure out why this normalisation works best - I also tried e.g. m*rownorms/sum(rownorms) but that didn't work so well
                                      return(out)
                                    },
                                    maxit = 30))

par(mfrow = c(3,1))
image(y = 0:m, x = 0:n, Y0^0.2, col = topo.colors(255), main = "Ground truth")
abline(v=u, col="green")
image(y = 0:m, x = 0:n, Y^0.2, col = topo.colors(255), main = "Input data")
abline(v=u, col="green")
image(y = 0:m, x = 0:n, fit2$fitted.values^0.2, col = topo.colors(255), main = "Fitted data")
abline(v=which(rowSums(fit2$coefficients)>1), col="green")












# tests using spams package, see ?spams.fistaFlat
library(spams)
?spams.fistaFlat
# Spams multi-task regression
# set.seed(0)
# m = 100;n = 200
# X = matrix(rnorm(m * n),nrow = m,ncol = n,byrow = FALSE)
# X = X - matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow = T)
# X = spams.normalize(X)
# Y = matrix(rnorm(100 * 100),nrow = 100,ncol = 100,byrow = FALSE)
# Y = Y - matrix(rep(colMeans(Y),nrow(Y)),nrow(Y),ncol(Y),byrow = T)
# Y = spams.normalize(Y)

W0 = matrix(c(0),nrow = ncol(X), ncol = ncol(Y))

print("FISTA nonnegative l1l2 penalized LS regression")
res = spams.fistaFlat(Y, X, W0, TRUE, numThreads = 1, verbose = TRUE, 
                      lambda1 = 1000, lambda2 = 0, lambda3 = 0,  
                      max_it = 100, L0 = 1, tol = 1e-5, intercept = FALSE, pos = TRUE,
                      compute_gram = FALSE, loss = "square", regul = "l1l2", transpose = FALSE,
                      ista = FALSE, subgrad = FALSE, 
                      a = 0.1, b = 100, size_group = 10)
image(res[[1]]^0.2)
image(S^0.2)
optim_info = res[[2]]
.printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info
        [1],optim_info[3],optim_info[4])


print("FISTA nonnegative l1linf penalized LS regression")
res = spams.fistaFlat(Y, X, W0, TRUE, numThreads = 1, verbose = TRUE,
                      lambda1 = 0, lambda2 = 0, lambda3 = 0,
                      max_it = 100, L0 = 0, tol = 1e-5, intercept = FALSE, pos = TRUE,
                      compute_gram = FALSE, loss = "square", regul = "l1linf", transpose = FALSE,
                      ista = FALSE, subgrad = FALSE, 
                      a = 0.1, b = 100, size_group = 10)
image(res[[1]])
optim_info = res[[2]]
.printf("mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n",optim_info
        [1],optim_info[3],optim_info[4])

