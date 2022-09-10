rm(list=ls())

library(glmnet)
library(flare)
library(MASS)
library(readstata13)
source("functions.txt")
dat=read.dta13("age5064female.dta")

set.seed(20201960)

D = dat$D
Z = dat$IV
Y=dat$pcs8_score

age=dat$age
black=dat$black
hispanic=dat$hispanic
Dage=D*age
Dage2=D*(age^2)
Dhispanic=D*hispanic
Dblack=D*black
which(names(dat)=='pcs8_score')
which(names(dat) =='IV')
which(names(dat)=='age')
which(names(dat)=='D')
which(names(dat)=='black')
which(names(dat)=='hispanic')
X=as.matrix(dat[,-c(4,48,2,47,45,43)])
Z1 = Z[!is.na(rowSums(X))]
D1=D[!is.na(rowSums(X))]
Y1 = Y[!is.na(rowSums(X))]
age1=age[!is.na(rowSums(X))]
black=black[!is.na(rowSums(X))]
hispanic=hispanic[!is.na(rowSums(X))]
Dage=Dage[!is.na(rowSums(X))]
Dage2=Dage2[!is.na(rowSums(X))]
Dhispanic=Dhispanic[!is.na(rowSums(X))]
Dblack=Dblack[!is.na(rowSums(X))]
X = X[!is.na(rowSums(X)), ]

which( apply(X, 2, sd) == 0)
XX=X
#generate polynomials
XX = as.data.frame(XX)
n=dim(XX)[1]
intercept = rep(1, n)
XX= cbind(intercept, XX)
which(apply(XX,2,sd)==0)
p0 = dim(XX)[2]
m = 0
W = c()
name.W = c()
name.X = names(XX)
for (i in 1 : p0){
  for (j in i : p0){
    m = m + 1
    Wtemp = XX[, i] * XX[, j]
    name.W[m] = paste(name.X[i], name.X[j], sep = "-")
    W = cbind(W, Wtemp)
  }
}   

W = data.frame(W)
names(W) = name.W
X1=as.matrix(W[,(-1)])
dim(X1)
which( apply(X1, 2, sd) == 0)

X=cbind(age1,hispanic,black,Dage, Dage2, Dhispanic, Dblack,X1)
XD = cbind(D1, X)
Yc = Y1 - mean(Y1)
Xc = t(t(XD) - colMeans(XD))
Ys = Yc / sd(Yc)
Xs = t(t(Xc) / apply(Xc, 2, sd))
Xs1 = Xs[, -1]
distr.quantile = 1.96
M = dim(Xs1)
n = M[1]; p = M[2]
s1.limit = ceiling(M[2] / 10)

fit.bin.1 = cv.glmnet(Xs1, Z1, alpha = 1, family = "binomial")
lambda1 = fit.bin.1$lambda.min
fit.bin = glmnet(Xs1, Z1, alpha = 1, family = "binomial", lambda = lambda1)
gamma.est.sparse = fit.bin$beta
intercept.est = fit.bin$a0

prob.est = exp(intercept.est + as.matrix(Xs1) %*% gamma.est.sparse) / (1 + exp(intercept.est + as.matrix(Xs1) %*% gamma.est.sparse))
weight0 = prob.est * (1 - prob.est)

Theta = matrix(1, p, p)
for (j in 1 : p){
  fit.nodewise.1 = glmnet(Xs1[, -j], Xs1[, j], alpha = 1, weights = weight0, family = "gaussian")
  temp1.nodewise = length(fit.nodewise.1$df[fit.nodewise.1$df <= s1.limit])
  lambda1.nodewise = fit.nodewise.1$lambda[temp1.nodewise]
  fit.nodewise = glmnet(Xs1[, -j], Xs1[, j], alpha = 1, weights = weight0, family = "gaussian", lambda = lambda1.nodewise)
  gamma.nodewise.est = fit.nodewise$beta
  if (j == 1) gamma.nodewise.row = c(1, -as.vector(gamma.nodewise.est))
  if (j == p) gamma.nodewise.row = c(-as.vector(gamma.nodewise.est), 1)
  if (j >= 2 & j < p) gamma.nodewise.row = c(-as.vector(gamma.nodewise.est)[1 : (j - 1)], 1, -as.vector(gamma.nodewise.est)[j : (p - 1)])
  tau = sum(weight0 * (Xs1[, j] - Xs1[, -j] %*% gamma.nodewise.est - fit.nodewise$a0)^2) / n + lambda1.nodewise * sum(abs(gamma.nodewise.est))
  Theta[j, ] = gamma.nodewise.row / tau
}

x1w = Xs1 * as.vector(Z1 - prob.est)
gamma.est.DBias = as.vector(gamma.est.sparse) + as.vector(Theta %*% matrix(colMeans(x1w), p, 1))
gamma.var = diag(Theta %*% (t(x1w) %*% x1w / n) %*% t(Theta))
gamma.stderr = sqrt(gamma.var / n)
# 
# #-------------------- Second stage -----------------

kappa0 = (1 - D1) * ((1 - Z1) - (1 - prob.est)) / (prob.est * (1 - prob.est))
kappa1 = D1 * (Z1 - prob.est) / (prob.est * (1 - prob.est))
kappa = kappa0 * (1 - prob.est) + kappa1 * prob.est

WSigma = 0
for (i in 1 : n){
  WSigma = WSigma + kappa[i] * matrix(Xs[i, ], p + 1, 1) %*% matrix(Xs[i, ], 1, p + 1)
}
WSigma = WSigma / n
s1=sum(gamma.est.sparse != 0)
lambda.second.range = seq(3 * sqrt(s1 * log(p) / n), 3 * sqrt(s1 * log(p) / n), 0.05)
lambda.clime =  sqrt(4 * log(p) / n)

beta.initial=rep(0,dim(Xs)[2])

fit.unweight.cv = cv.glmnet(Xs, Ys, alpha = 1, weights = rep(1, n), family = "gaussian", standardize = FALSE, intercept = FALSE) # new
lambda.unweight = fit.unweight.cv$lambda.min 
fit.unweight = glmnet(Xs, Ys, alpha = 1, weights = rep(1, n), family = "gaussian", lambda = lambda.unweight, standardize = FALSE, intercept = FALSE) 

L1bound = seq(1,sqrt(n),4)
# 
fit.second = glmnet.PGD.cv(Xs, Ys, kappa, lambda.second.range, L1bound, beta.initial, rate = 0.01)
beta.est.1 = fit.second$slope
steps.second = fit.second$step.convergence

WSigma.inv.out = sugm(WSigma, lambda = lambda.clime, method = "clime")
WSigma.inv = WSigma.inv.out$icov[[1]]
WSigma.df = WSigma.inv.out$df

Y.pred = mean(Ys) + Xs %*% beta.est.1
residual = Ys - Y.pred

W.residual = 0
N.residual = 0
for (i in 1 : n){
  W.residual = W.residual + kappa[i] * residual[i] * as.vector(Xs[i, ])
  N.residual = N.residual + residual[i] * Kappa.D1(Z[i], D[i], prob.est[i]) * prob.est[i] * (1 - prob.est[i]) * matrix(Xs[i, ], p + 1, 1) %*% matrix(Xs1[i, ], 1, p)
}
W.residual  = W.residual / n
N.residual  = N.residual / n

threshold.q = 0.9 
threshold.1 = quantile(abs(N.residual), threshold.q, names = FALSE) 
N.residual.threshold = N.residual
N.residual.threshold[abs(N.residual.threshold) < threshold.1] = 0 

Est.var.sandwich = 0
for (i in 1:n) {
  Est.var.temp = kappa[i] * residual[i] * as.vector(Xs[i,]) + as.vector(N.residual.threshold %*% Theta %*% x1w[i,])
  Est.var.sandwich = Est.var.sandwich + matrix(Est.var.temp, p + 1, 1) %*% matrix(Est.var.temp, 1, p + 1)
}
Est.var.sandwich = Est.var.sandwich / n
Est.var = diag(WSigma.inv %*% Est.var.sandwich %*% WSigma.inv)
Est.stdErr = sqrt(Est.var / n)

beta.est.DBias = beta.est.1 + as.vector(
  WSigma.inv %*% W.residual - WSigma.inv %*% N.residual.threshold %*% (as.vector(gamma.est.sparse) - gamma.est.DBias)
) 
lower = beta.est.DBias - distr.quantile * Est.stdErr
upper = beta.est.DBias + distr.quantile * Est.stdErr
t = beta.est.DBias / Est.stdErr

pvalue=	2 * (1 - pnorm(abs(t)))
result=cbind(beta.est.DBias, Est.stdErr,t,pvalue)
cutoff1= max(which(1 * (pvalue <= 0.10 * c(1 : (p + 1)) / (p + 1)) == 1))
cutoff2 = pvalue[cutoff1]
significant.Threshold = which(pvalue <= cutoff2)
s2=length(significant.Threshold)
#write.csv(result, file = "pcs8_s2_126_p996.csv")
