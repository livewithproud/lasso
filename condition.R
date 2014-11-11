#defination
n <- 80
p <- 100
iterate <- 100
d <- 1
# b <- matrix(c(3,1.5,0,0,2,0,0,0),ncol=1)
# b <- matrix(c(0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85),ncol=1)
# b <- matrix(c(5,0,0,0,0,0,0,0),ncol=1)
b <- matrix(c(rep(0,90),sample(10:30,10,replace = TRUE)/10),ncol = 1)
rou <- 0.5
sigma <- diag(p)
# sigma <- matrix(0, nrow = p, ncol = p)
# for (i in 1:p)
# {
# 	for(j in i:p)
# 	{
# 		temp <- abs(i - j)
# 		sigma[j,i] <- sigma[i,j] <-rou^temp 
# 	}
# }
# for (i in 1:p)
# {
# 	sigma[i,i] <- 2
# 	for(j in i:p)
# 	{
# 		sigma[j,i] <- sigma[i,j] <-0.5 
# 	}
# }

l <- simulation(n, p, iterate, d, b, sigma)

RMSq <- f$RSS[seq(20)]/(n + eps -  f$df[seq(20)])
Rq <-1- ((n-1)* f$RSS[seq(20)])/(f$TSS*(n-f$df[seq(20)])) 
Sq <- f$RSS[seq(20)]/((n - f$df[seq(20)])*(n - f$df[seq(20)] - 1))

AIC <- n * log(f$RSS) + 2 * f$df
BIC <- n * log(f$RSS) + log(n) * f$df

Cp <- f$RSS/(rev(f$RSS)[1]/n) - n + 2 * f$df