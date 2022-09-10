library(glmnet)
library(flare)
library(MASS)
source("functions.txt")
source("main.txt")

p = 100
n = 400
rho.X = 0.5
sigma.X = 1

gamma1 = 1
gamma0 = 0

alpha = 1
beta1.true = 2

r1 = main2(n, p, rho.X, sigma.X, gamma1, gamma0, lambda.NC = c(0.5, 0.5), alpha, beta1.true, 1000)

r2 = main2(n, p, rho.X, sigma.X, gamma1, gamma0, lambda.NC = c(0.5, 0.75), alpha, beta1.true, 1000)

r3 = main2(n, p, rho.X, sigma.X, gamma1, gamma0, lambda.NC = c(0.5, 1), alpha, beta1.true, 1000)

r4 = main2(n, p, rho.X, sigma.X, gamma1, gamma0, lambda.NC = c(0.5, 1.25), alpha, beta1.true, 1000)

write.csv(r1$CI, file = "CI-p100-n400-dgp2-0505.csv")
write.csv(r2$CI, file = "CI-p100-n400-dgp2-05075.csv")
write.csv(r1$CI, file = "CI-p100-n400-dgp2-0510.csv")
write.csv(r2$CI, file = "CI-p100-n400-dgp2-05125.csv")
