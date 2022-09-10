library(glmnet)
library(flare)
library(MASS)
source("functions.txt")
source("main.txt")

p = 100
n = 400
rho.X = 0.5
sigma.X = 0.5

s1 = 5
gamma1 = 0.5
gamma0 = 0

alpha = 1
s2 = 5
beta1.true = 1

r10 = main1(n, p, rho.X, sigma.X, s1, gamma1, gamma0, epsilon.rho = 0, alpha, s2, beta1.true, 1000)

r12 = main1(n, p, rho.X, sigma.X, s1, gamma1, gamma0, epsilon.rho = 0.2, alpha, s2, beta1.true, 1000)

r13 = main1(n, p, rho.X, sigma.X, s1, gamma1, gamma0, epsilon.rho = 0.3, alpha, s2, beta1.true, 1000)

r14 = main1(n, p, rho.X, sigma.X, s1, gamma1, gamma0, epsilon.rho = 0.4, alpha, s2, beta1.true, 1000)

write.csv(r10, file = "CI-p100-n400-dgp1-dep0.csv")
write.csv(r12, file = "CI-p100-n400-dgp1-dep2.csv")
write.csv(r13, file = "CI-p100-n400-dgp1-dep3.csv")
write.csv(r14, file = "CI-p100-n400-dgp1-dep4.csv")
