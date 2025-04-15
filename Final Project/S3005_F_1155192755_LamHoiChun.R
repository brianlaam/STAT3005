###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###
### Stat 3005 -- Final (2024)
###
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------

#--------------------------------------------
# Question 1
#--------------------------------------------
# Import Data
data = read.csv("C:\\Users\\brian\\Documents\\CUHK\\24-25 Term1\\STAT3005\\Final Project\\Q1.csv")

# Simulation parameters
Nrep = 10000  
alpha = 0.05
n = 30

#--------------------------------------------
# Question 1.1 Rank Sum Test
#--------------------------------------------
# Test Function
rankSum.test = function(x, y, alternative=c("two.sided", "less", "greater"), 
                        mu=0, exact=FALSE, correct=TRUE, print=FALSE){	
  alternative = match.arg(alternative) 
  xy = c(x,y+mu)
  nj = table(xy)
  A  = sum(nj^3-nj)
  n  = length(xy)
  ties = n>length(nj)
  n1 = length(x)
  n2 = n-n1
  R  = rank(xy)
  T  = sum(R[1:n1])/(n+1)
  muT = n1/2
  sdT = sqrt((n1*n2)/(12*(n+1)) * (1-A/(n^3-n)))
  if(correct&(!(exact&(!ties)))){
    if(alternative=="less"){ T=T+0.5/(n+1)}
    if(alternative=="greater"){ T=T-0.5/(n+1)}
    if(alternative=="two.sided"){ T=T-0.5*sign(T-muT)/(n+1)}
    
  }
  t = (T-muT)/sdT
  if(exact&(!ties)){
    nRep = 2^16
    out = rep(NA,nRep)
    for(iRep in 1:nRep){
      xy.sim = rt(n,1)
      R.sim = rank(xy.sim)
      T.sim = sum(R.sim[1:n1])/(n+1)
      t.sim = (T.sim-muT)/sdT
      out[iRep] = t.sim
    }
    p = switch(alternative, "less"=mean(out<=t),
               "greater"=mean(out>=t),
               "two.sided"=mean(out^2>=t^2))
  }else{
    p = switch(alternative, "less"=pnorm(t),
               "greater"=1-pnorm(t),
               "two.sided"=2*pnorm(-abs(t)))
  }
  if(print){
    cat(" Rank sum test\n-------------------------------------\n",
        paste0("H0: location(x)-location(y)=",mu,"\n"),
        paste0("H1: location(x)-location(y)",
               switch(alternative,"less"="<","greater"="<","two.sided"="!="),mu),
        "\n-------------------------------------\n")
  }
  list(statistic=t, p.value=p)
}

# (i) Size under H_0
set.seed(3005)
reject_count_H0 = 0
for(i in 1:Nrep){
  x = rnorm(n, mean=1, sd=1)
  y = rnorm(n, mean=1, sd=1)
  test_res = rankSum.test(x, y, alternative="two.sided")
  if(test_res$p.value < alpha) {
    reject_count_H0 = reject_count_H0 + 1
  }
}
size_estimate = reject_count_H0 / Nrep
cat("Estimated size under H0:", size_estimate, "\n")
# Estimated size under H0: = 0.0498 

# (ii) Power 
theta_x_values = c(1.0, 1.5, 2, 2.5)
power_estimates = numeric(length(theta_x_values))
set.seed(3005)
for(j in seq_along(theta_x_values)){
  theta_x = theta_x_values[j]
  reject_count_H1 = 0
  for(i in 1:Nrep){
    x = rnorm(n, mean=theta_x, sd=theta_x)
    y = rnorm(n, mean=1, sd=1)
    test_res = rankSum.test(x, y, alternative="two.sided")
    if(test_res$p.value < alpha){
      reject_count_H1 = reject_count_H1 + 1
    }
  }
  power_estimates[j] = reject_count_H1 / Nrep
}
list(theta = theta_x_values, Power = power_estimates)
# $theta
# [1] 1.0 1.5 2.0 2.5
# 
# $Power
# [1] 0.0498 0.3064 0.6425 0.8054


# (iii) Performing Test 
rankSum.test(data$x1, data$x2, alternative = "two.sided")$p.value
# Output = 0.01173492



#--------------------------------------------
# Question 1.2 Ansari-Bradley Test
#-------------------------------------------- 
# (i) Size under H_0
set.seed(3005)
reject_count_H0 = 0
for(i in 1:Nrep){
  x = rnorm(n, mean=1, sd=1)
  y = rnorm(n, mean=1, sd=1)
  test_res = ansari.test(x, y, alternative="two.sided")
  if(test_res$p.value < alpha) {
    reject_count_H0 = reject_count_H0 + 1
  }
}
size_estimate = reject_count_H0 / Nrep
cat("Estimated size under H0:", size_estimate, "\n")
# Estimated size under H0: 0.0488  

# (ii) Power 
theta_x_values = c(1.0, 1.5, 2, 2.5)
power_estimates = numeric(length(theta_x_values))
set.seed(3005)
for(j in seq_along(theta_x_values)){
  theta_x = theta_x_values[j]
  reject_count_H1 = 0
  for(i in 1:Nrep){
    x = rnorm(n, mean=theta_x, sd=theta_x)
    y = rnorm(n, mean=1, sd=1)
    test_res = ansari.test(x, y, alternative="two.sided")
    if(test_res$p.value < alpha){
      reject_count_H1 = reject_count_H1 + 1
    }
  }
  power_estimates[j] = reject_count_H1 / Nrep
}
list(theta = theta_x_values, Power = power_estimates)

# $theta
# [1] 1.0 1.5 2.0 2.5
# 
# $Power
# [1] 0.0488 0.3367 0.6839 0.8367

# (iii) Performing Test
ansari.test(data$x1, data$x2, alternative = "two.sided")$p.value
# Output = 0.1401889

#--------------------------------------------
# Question 1.3 Calvin's Test
#-------------------------------------------- 
# (i) Size under H_0
set.seed(3005)
reject_count_H0 = 0
for (i in 1:Nrep) {
  x = rnorm(n, mean=1, sd=1)
  y = rnorm(n, mean=1, sd=1)
  
  rs_test = rankSum.test(x, y, alternative = "two.sided")
  ab_test = ansari.test(x, y, alternative = "two.sided")
  
  # Combined test rejects if either rejects
  if (rs_test$p.value < alpha || ab_test$p.value < alpha) {
    reject_count_H0 = reject_count_H0 + 1
  }
}

size_estimate = reject_count_H0 / Nrep
cat("Estimated size under H0:", size_estimate, "\n")
# Estimated size under H0: 0.0965

# (ii) Power 
theta_x_values = c(1.0, 1.5, 2, 2.5)
power_estimates = numeric(length(theta_x_values))
set.seed(3005)
for(j in seq_along(theta_x_values)){
  theta_x = theta_x_values[j]
  reject_count_H1 = 0
  for(i in 1:Nrep){
    x = rnorm(n, mean=theta_x, sd=theta_x)
    y = rnorm(n, mean=1, sd=1)
    rs_test = rankSum.test(x, y, alternative = "two.sided")
    ab_test = ansari.test(x, y, alternative = "two.sided")
    if (rs_test$p.value < alpha || ab_test$p.value < alpha) {
      reject_count_H1 = reject_count_H1 + 1
    }
  }
  power_estimates[j] = reject_count_H1 / Nrep
}
list(theta = theta_x_values, Power = power_estimates)

# $theta
# [1] 1.0 1.5 2.0 2.5
# 
# $Power
# [1] 0.0965 0.5575 0.9266 0.9901

# (iii) Performing Test
min(rankSum.test(data$x1, data$x2, alternative = "two.sided")$p.value, 
    ansari.test(data$x1, data$x2, alternative = "two.sided")$p.value)
# Output = 0.01173492


#--------------------------------------------
# Question 1.4 Cramer-von Mises Test
#-------------------------------------------- 
# Test Function
cvm.test0 = function(x, F0, exact=TRUE){
  n = length(x)
  u = F0(sort(x))
  i = 1:n
  C = 1/12/n + sum((u-(2*i-1)/2/n)^2)
  if(exact){
    nRep = 10000
    out = rep(NA,nRep)
    for(iRep in 1:nRep) out[iRep] = cvm.test0(runif(n), F0=punif, exact=FALSE)$statistic
    p = mean(out>C)
  }else{
    p = NA
  }
  list(statistic=C, p.value=p)
}

# (i) Size under H_0
set.seed(3005)
Nrep = 500 # Lower the number of rep due to bad performance of the function
reject_count_H0 = 0
for(i in 1:Nrep){
  x = rnorm(n, mean=1, sd=1)
  y = rnorm(n, mean=1, sd=1)
  test_res = cvm.test0(x, ecdf(y))
  if(test_res$p.value < alpha) {
    reject_count_H0 = reject_count_H0 + 1
  }
  if(i %% 10 == 0) {
    cat("Simulation process:", i, "\n")
  }
}
size_estimate = reject_count_H0 / Nrep
cat("Estimated size under H0:", size_estimate, "\n")
# Estimated size under H0: 0.234  

# (ii) Power 
Nrep = 100 # Lower the number of rep due to bad performance of the function
theta_x_values = c(1.0, 1.5, 2, 2.5)
power_estimates = numeric(length(theta_x_values))
set.seed(3005)
for(j in seq_along(theta_x_values)){
  theta_x = theta_x_values[j]
  reject_count_H1 = 0
  for(i in 1:Nrep){
    x = rnorm(n, mean=theta_x, sd=theta_x)
    y = rnorm(n, mean=1, sd=1)
    test_res = cvm.test0(x, ecdf(y))
    if(test_res$p.value < alpha){
      reject_count_H1 = reject_count_H1 + 1
    }
    if(i %% 10 == 0) {
      cat("Simulation process:", i, "\n")
    }
  }
  power_estimates[j] = reject_count_H1 / Nrep
}
list(theta = theta_x_values, Power = power_estimates)

# $theta
# [1] 1.0 1.5 2.0 2.5
# 
# $Power
# [1] 0.24 0.64 0.97 1.00

# (iii) Performing Test
set.seed(3005)
cvm.test0(data$x1, ecdf(data$x2))$p.value
# Output: 3e-04

#--------------------------------------------
# Question 1.5 ks.test() function
#--------------------------------------------
# (i) Size under H_0
set.seed(3005)
reject_count_H0 = 0
for(i in 1:Nrep){
  x = rnorm(n, mean=1, sd=1)
  y = rnorm(n, mean=1, sd=1)
  test_res = ks.test(x, ecdf(y))
  if(test_res$p.value < alpha) {
    reject_count_H0 = reject_count_H0 + 1
  }
  if(i %% 1000 == 0) {
    cat("Simulation process:", i, "\n")
  }
}
size_estimate = reject_count_H0 / Nrep
cat("Estimated size under H0:", size_estimate, "\n")
# Estimated size under H0: 0.2437  

# (ii) Power
theta_x_values = c(1.0, 1.5, 2, 2.5)
power_estimates = numeric(length(theta_x_values))
set.seed(3005)
for(j in seq_along(theta_x_values)){
  theta_x = theta_x_values[j]
  reject_count_H1 = 0
  for(i in 1:Nrep){
    x = rnorm(n, mean=theta_x, sd=theta_x)
    y = rnorm(n, mean=1, sd=1)
    test_res = ks.test(x, ecdf(y))
    if(test_res$p.value < alpha){
      reject_count_H1 = reject_count_H1 + 1
    }
    if(i %% 1000 == 0) {
      cat("Simulation process:", i, "\n")
    }
  }
  power_estimates[j] = reject_count_H1 / Nrep
}
list(theta = theta_x_values, Power = power_estimates)

# $theta
# [1] 1.0 1.5 2.0 2.5
# 
# $Power
# [1] 0.2437 0.7003 0.9635 0.9959

# (iii) Performing Test
ks.test(data$x1, ecdf(data$x2))$p.value
# Output: 0.0001908552

#--------------------------------------------
# Question 1.6 2-sample Permutation test
#--------------------------------------------
# Test function
group.ptest = function(x, y, FUN=NULL, B=1e4, alternative=c("two.sided","less","greater"), plot=FALSE){
  # FUN         = function of test statistic 
  # mu0         = null value of mu
  # B           = number of Monte Carlo replications
  # alternative = direction of alternative hypothesis (e.g., "less" means H1: mu(x)<mu(y))
  # plot        = logical indicating whether to plot permutation distribution
  alternative = match.arg(alternative)
  if(is.null(FUN)) FUN = function(x,y) sum(rank(c(x,y))[1:length(x)])
  n1 = length(x)
  n2 = length(y)
  xy = c(x,y)
  n = n1+n2
  Tobs = FUN(x,y)
  Tb = rep(NA,B)
  for(b in 1:B){
    I = sample(1:n,n1,replace=FALSE)
    xb = xy[I]
    yb = xy[-I]
    Tb[b] = FUN(xb,yb)
  }
  if(plot){
    hist(Tb,nclass=30,freq=FALSE,xlim=range(c(Tb,Tobs)),xlab=expression(italic(T)),
         main="Permutation distribution",
         col="pink",border=FALSE,cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
    abline(v=Tobs,col="blue4",lty=2,lwd=2)
  }
  switch(alternative, "two.sided"=min(mean(Tb<=Tobs),mean(Tb>=Tobs))*2,
         "less"=mean(Tb<=Tobs),
         "greater"=mean(Tb>=Tobs))
}

# (i) Size under H_0
set.seed(3005)
Nrep = 500 # Lower the number of rep due to bad performance of the function
FUN_tTest = function(x,y) mean(x) - mean(y) # Permutation t-test use mean diff as test statistic
reject_count_H0 = 0
for(i in 1:Nrep){
  x = rnorm(n, mean=1, sd=1)
  y = rnorm(n, mean=1, sd=1)
  test_res = group.ptest(x, y, FUN=FUN_tTest)
  if(test_res < alpha) {
    reject_count_H0 = reject_count_H0 + 1
  }
  if(i %% 100 == 0) {
    cat("Simulation process:", i, "\n")
  }
}
size_estimate = reject_count_H0 / Nrep
cat("Estimated size under H0:", size_estimate, "\n")
# Estimated size under H0: 0.058


# (ii) Power
theta_x_values = c(1.0, 1.5, 2, 2.5)
power_estimates = numeric(length(theta_x_values))
FUN_tTest = function(x,y) mean(x) - mean(y) # Permutation t-test use mean diff as test statistic
set.seed(3005)
for(j in seq_along(theta_x_values)){
  theta_x = theta_x_values[j]
  reject_count_H1 = 0
  for(i in 1:Nrep){
    x = rnorm(n, mean=theta_x, sd=theta_x)
    y = rnorm(n, mean=1, sd=1)
    test_res = group.ptest(x, y, FUN=FUN_tTest)
    if(test_res < alpha){
      reject_count_H1 = reject_count_H1 + 1
    }
    if(i %% 100 == 0) {
      cat("Simulation process:", i, "\n")
    }
  }
  power_estimates[j] = reject_count_H1 / Nrep
}
list(theta = theta_x_values, Power = power_estimates)

# $theta
# [1] 1.0 1.5 2.0 2.5
# 
# $Power
# [1] 0.058 0.332 0.648 0.820

# (iii) Performing Test
set.seed(3005)
group.ptest(data$x1, data$x2,FUN=FUN_tTest)
# Output: 0.0096

#--------------------------------------------
# Question 1.7 Paired sign-rank test
#--------------------------------------------
# (i) Size under H_0
set.seed(3005)
reject_count_H0 = 0
for(i in 1:Nrep){
  x = rnorm(n, mean=1, sd=1)
  y = rnorm(n, mean=1, sd=1)
  test_res = wilcox.test(x, y, mu=0, paired=TRUE, exact=TRUE, alternative="two.sided")$p.value
  if(test_res < alpha) {
    reject_count_H0 = reject_count_H0 + 1
  }
  if(i %% 1000 == 0) {
    cat("Simulation process:", i, "\n")
  }
}
size_estimate = reject_count_H0 / Nrep
cat("Estimated size under H0:", size_estimate, "\n")
# Estimated size under H0: 0.0529  

# (ii) Power
theta_x_values = c(1.0, 1.5, 2, 2.5)
power_estimates = numeric(length(theta_x_values))
set.seed(3005)
for(j in seq_along(theta_x_values)){
  theta_x = theta_x_values[j]
  reject_count_H1 = 0
  for(i in 1:Nrep){
    x = rnorm(n, mean=theta_x, sd=theta_x)
    y = rnorm(n, mean=1, sd=1)
    test_res = wilcox.test(x, y, mu=0, paired=TRUE, exact=TRUE, alternative="two.sided")$p.value
    if(test_res < alpha){
      reject_count_H1 = reject_count_H1 + 1
    }
    if(i %% 1000 == 0) {
      cat("Simulation process:", i, "\n")
    }
  }
  power_estimates[j] = reject_count_H1 / Nrep
}
list(theta = theta_x_values, Power = power_estimates)

# $theta
# [1] 1.0 1.5 2.0 2.5
# 
# $Power
# [1] 0.0529 0.3031 0.6412 0.8169


# (iii) Performing Test
wilcox.test(data$x1, data$x2, mu=0, paired=TRUE, exact=TRUE, alternative="two.sided")$p.value
# Output = 0.004101777

#--------------------------------------------
# Question 1.8 Bootstrap CI
#--------------------------------------------
# Test Function
ci.bcaboot_2sample = function(x, y, FUN, alpha0=.05, B=10000){
  n1 = length(x)
  n2 = length(y)
  est = FUN(x,y)
  
  # Bootstrap replicates
  boot = numeric(B)
  for(b in 1:B){
    x_star = sample(x, n1, replace=TRUE)
    y_star = sample(y, n2, replace=TRUE)
    boot[b] = FUN(x_star, y_star)
  }
  
  jack_x = numeric(n1)
  for(i in 1:n1) jack_x[i] = FUN(x[-i], y)

  jack_y = numeric(n2)
  for(j in 1:n2) jack_y[j] = FUN(x, y[-j])
  
  jack = c(jack_x, jack_y)
  jackBar = mean(jack)
  
  A = sum((jack - jackBar)^3)/(6 * (sum((jack - jackBar)^2)^(3/2)))
  Z = qnorm(mean(boot < est))
  alpha1 = pnorm(Z + (Z + qnorm(alpha0/2)) / (1 - A*(Z + qnorm(alpha0/2))))
  alpha2 = pnorm(Z + (Z + qnorm(1 - alpha0/2)) / (1 - A*(Z + qnorm(1 - alpha0/2))))
  
  CI = quantile(boot, c(alpha1, alpha2))
  names(CI) = c("lower","upper")
  CI
}
FUN_diffvar = function(x,y) var(x)-var(y)

# (i) Size under H_0
Nrep = 500 # Lower the number of rep due to bad performance of the function
set.seed(3005)
reject_count_H0 = 0
for(i in 1:Nrep){
  x = rnorm(n, mean=1, sd=1)
  y = rnorm(n, mean=1, sd=1)
  CI = ci.bcaboot_2sample(x,y,FUN=FUN_diffvar,alpha0=alpha,B=2000) # B=2000 for speed
  if(CI[1]>0 || CI[2]<0) reject_count_H0 = reject_count_H0 + 1
  if(i %% 100 == 0) {
    cat("Simulation process:", i, "\n")
  }
}

size_estimate = reject_count_H0 / Nrep
cat("Estimated size under H0:", size_estimate, "\n")
# Estimated size under H0: 0.072

# (ii) Power
set.seed(3005)
Nrep = 500 # Lower the number of rep due to bad performance of the function
theta_x_values = c(1.0, 1.5, 2, 2.5)
power_estimates = numeric(length(theta_x_values))

for(j in seq_along(theta_x_values)){
  theta_x = theta_x_values[j]
  reject_count_H1 = 0
  for(i in 1:Nrep){
    x = rnorm(n, mean=theta_x, sd=theta_x)
    y = rnorm(n, mean=1, sd=1)
    
    CI = ci.bcaboot_2sample(x,y,FUN=FUN_diffvar,alpha0=alpha,B=2000)
    if(CI[1]>0 || CI[2]<0) reject_count_H1 = reject_count_H1 + 1
    if(i %% 100 == 0) {
      cat("Simulation process:", i, "\n")
    }
  }
  power_estimates[j] = reject_count_H1 / Nrep
}
list(theta = theta_x_values, Power = power_estimates)

# $theta
# [1] 1.0 1.5 2.0 2.5
# 
# $Power
# [1] 0.072 0.590 0.970 1.000

# (iii) Performing Test
set.seed(3005)
(CI = ci.bcaboot_2sample(data$x1,data$x2,FUN=FUN_diffvar,alpha0=alpha,B=2000))
#      lower      upper 
# -5541.1703   717.2781   


#--------------------------------------------
# Question 2 (Question Setting)
#--------------------------------------------
# mu =    c(-10, -3, -1, 0, 1, 3, 10)
# sigma = c(3, 6, 0.3, 0.3, 0.3, 6, 3)
# p =     c(1, 1, 2, 1, 2, 1, 1)/9
# t = seq(from=-20,to=20,length=301)
# f = colSums(simplify2array(lapply(t,function(t)dnorm(t, mu,sigma)))*p)
# plot(t,f,type="l")

#--------------------------------------------
# Question 2.1(b) Updated
#--------------------------------------------
ddnorm = function(x,mu,sigma,r=0){
  if(r==0){
    out = dnorm(x,mu,sigma)
  }
  if(r==1){
    out = dnorm(x,mu,sigma)*(-1/sigma^2)*(x-mu)
  }
  if(r>=2){
    out = -(r-1)/sigma^2*ddnorm(x,mu,sigma,r=r-2)-(x-mu)/sigma^2*ddnorm(x,mu,sigma,r=r-1)
  }
  out
}

# Compute lambda_K_r
compute_lambda_K_r = function(r) {
  integrate(function(t){(ddnorm(t,0,1,r=r))^2}, lower=-Inf, upper=Inf)$value
}

# Compute (f^{(r+2)}(t))^2 dt
compute_int_df_square = function(r, sigma_f) {
  integrate(function(t){(ddnorm(t,0,sigma_f,r=r+2))^2}, lower=-Inf, upper=Inf)$value
}

# Compute lambda_K_0
lambda_K_0 = compute_lambda_K_r(0)
# Output = 0.2820948

# Compute (f^{(0+2)}(t))^2 dt
int_df_2_square = compute_int_df_square(0,1)
# Output = 0.2115711

# Required coefficient
(lambda_K_0/int_df_2_square)^(1/5)
# Output = 1.059224

# Compute lambda_K_2
lambda_K_2 = compute_lambda_K_r(2)
# Output = 0.2115711

# Compute (f^{(0+2)}(t))^2 dt
int_df_4_square = compute_int_df_square(2,1)
# Output = 1.851247

# Required coefficient
(5*lambda_K_2/int_df_4_square)^(1/9)
# Output = 0.9397142

#--------------------------------------------
# Question 2.1(c)
#--------------------------------------------
kdeDerivative = function(xeval, data, r=0) {
  n = length(data)
  sigma_hat = sd(data)
  IQR_hat = IQR(data)
  s_hat = min(sigma_hat, IQR_hat/1.34)
  sigma_K_4 = 1
  
  # Compute lambda_K^{(r)}
  lambda_K_r = compute_lambda_K_r(r)
  
  # Compute integral for f^{(r+2)}
  int_df_square = compute_int_df_square(r, s_hat)
  
  # Compute ell_r*
  ell = (((2*r+1)*lambda_K_r)/(n*sigma_K_4*int_df_square))^(1/(2*r+5))*s_hat
  
  # Compute f^(r)(x) for each x in xeval
  # f^(r)(x) = (1/(n * ell^(r+1))) * sum over i of ddnorm((x - X_i)/ell,0,1,r)
  
  # Compute Z_ij = (x_j - X_i)/ell
  Z = outer(xeval, data, FUN=function(x, Xi) (x - Xi)/ell)
  
  # Compute the Summation of K^(r)(Z)
  K_r_Z = ddnorm(Z, 0, 1, r=r)
  sum_K_r = rowSums(K_r_Z)

  # Final output
  f_r_est = 1/(n*ell^(r+1))*sum_K_r
  f_r_est
}

#--------------------------------------------
# Question 2.2
#--------------------------------------------
kdeBallon = function(xeval, data) {
  n = length(data)
  lambda_K0 = compute_lambda_K_r(0)
  sigma_K_4 = 1
  
  # Compute the pilot estimates: f(x) and f''(x)
  f_est = kdeDerivative(xeval, data, r=0)
  f_2nd_est = kdeDerivative(xeval, data, r=2)
  
  # Compute ell*(x)
  ell_x = ((f_est*lambda_K0)/(n*(f_2nd_est^2)* sigma_K_4))^(1/5)
  
  # Compute the balloon estimator
  results = numeric(length(xeval))
  for (j in 1:length(xeval)) {
    x = xeval[j]
    ell_j = ell_x[j]
    Z = (x - data) / ell_j
    K_values = dnorm(Z, 0, 1)
    f_Bal_x = (1/(n * ell_j)) * sum(K_values)
    results[j] = f_Bal_x
  }
  results
}

# Plot
data = read.csv("C:\\Users\\brian\\Documents\\CUHK\\24-25 Term1\\STAT3005\\Final Project\\Q2(new).csv")
n = 3005
mu = c(-10, -3, -1, 0, 1, 3, 10)
sigma = c(3, 6, 0.3, 0.3, 0.3, 6, 3)
weights = c(1,1,2,1,2,1,1)/9
t = seq(from=-20,to=20,length=301)

# True density
f = colSums(simplify2array(lapply(t,function(t)dnorm(t, mu,sigma)))*p)

# Balloon estimate
f_ballon = kdeBallon(t, data$x)

# Rule of Thumb BW
f_rot = density(data$x, bw=bw.nrd0(data$x), from=-20, to=20, n=301)$y

# Cross Validation BW
f_cv = density(data$x, bw=bw.ucv(data$x), from=-20, to=20, n=301)$y

# Plug-in BW
f_pin = density(data$x, bw=bw.SJ(data$x), from=-20, to=20, n=301)$y

# Plot 
plot(t, f, type="l", lwd=2, col="black", ylim=c(0, max(f, f_ballon, f_rot, f_cv, f_pin)),
     main="Comparison of True Density, Balloon Estimator, and other KDEs",
     xlab="x", ylab="Density")
lines(t, f_ballon, col="red", lwd=2, lty=1)
lines(t, f_rot, col="goldenrod2", lwd=2, lty=2)
lines(t, f_cv, col="steelblue3", lwd=2, lty=3)
lines(t, f_pin, col="seagreen4", lwd=2, lty=4)
legend("topright", legend=c("True Density", "Balloon Estimator", "ROT", "CV", "Plug-in"),
       col=c("black","red","goldenrod2","steelblue3","seagreen4"), lty=c(1,2,3,4), lwd=1, cex=0.75)

#--------------------------------------------
# Question 2.3
#--------------------------------------------
kdeSample = function(xeval, data) {
  n = length(data)
  
  # Plug-in BW
  ell_0 = bw.SJ(data)
  
  # Compute the kernel input
  Z = outer(data, data, FUN=function(xi, xj) (xi - xj)/ell_0)
  phi_Z = dnorm(Z)
  
  hat_f_ell0_Xi = rowSums(phi_Z) / (n * ell_0)
  bar_f_GM = exp(mean(log(hat_f_ell0_Xi)))
  ell_i = ell_0 * sqrt(bar_f_GM / hat_f_ell0_Xi)

  f_sam = numeric(length(xeval))
  for (k in 1:length(xeval)) {
    x = xeval[k]
    Zx = (x - data)/ell_i
    phi_Zx = dnorm(Zx)
    f_sam[k] = mean(phi_Zx / ell_i)
  }
  f_sam
}

# Density for sample-point estimator
f_sam = kdeSample(t, data$x)

# Plot 
plot(t, f, type="l", lwd=2, col="black", ylim=c(0, max(f, f_sam, f_rot, f_cv, f_pin)),
     main="Comparison of True Density, Sample-point Estimator, and other KDEs",
     xlab="x", ylab="Density")
lines(t, f_sam, col="red", lwd=2, lty=1)
lines(t, f_rot, col="goldenrod2", lwd=2, lty=2)
lines(t, f_cv, col="steelblue3", lwd=2, lty=3)
lines(t, f_pin, col="seagreen4", lwd=2, lty=4)
legend("topright", legend=c("True Density", "Sample-point Estimator", "ROT", "CV", "Plug-in"),
       col=c("black","red","goldenrod2","steelblue3","seagreen4"), lty=c(1,2,3,4), lwd=1, cex=0.75)


#--------------------------------------------
# Question 3.6
#--------------------------------------------
data = read.csv("C:\\Users\\brian\\Documents\\CUHK\\24-25 Term1\\STAT3005\\Final Project\\Q3.csv")

kw.test = function(x, g, method = c("asymptotic", "exact", "permutation"), B = 1000) {
  method = match.arg(method)

  if (!is.factor(g)) {
    g = factor(g)
  }
  
  N = length(x)
  group_sizes = table(g)
  k = length(group_sizes)
  ranks = rank(x)
  
  # Compute T_obs
  Rbar_dotdot = (N+1)/2
  Rbar_j = tapply(ranks, g, mean)
  n_j = as.numeric(group_sizes)
  T_obs = (12/(N*(N+1)))*sum(n_j*(Rbar_j-Rbar_dotdot)^2)
  
  # T-stat Function
  kw_stat = function(r, grp, N) {
    Rbar_dotdot = (N+1)/2
    Rbar_j = tapply(r, grp, mean)
    n_j = table(grp)
    (12/(N*(N+1)))*sum(n_j*(Rbar_j-Rbar_dotdot)^2)
  }
  
  if (method == "asymptotic") {
    # Built-in kruskal.test() function is used for asymptotic p-value
    kt = kruskal.test(x, g)
    return(list(
      statistic = T_obs,
      p.value = kt$p.value,
      method = "Asymptotic Kruskal-Wallis Test"
    ))
  }
  
  # Permutation Test Function
  perform_permutation_test = function(ranks, group_sizes, k, B) {
    N = length(ranks)
    T_values = numeric(B)
    for (b in seq_len(B)) {
      permutated_ranks = sample(ranks, size=N, replace=FALSE)
      permutated_groups = factor(rep(seq_len(k), times=group_sizes), levels=seq_len(k))
      T_values[b] = kw_stat(permutated_ranks, permutated_groups, N)
    }
    return(T_values)
  }
  
  if (method == "exact") {
    T_values = perform_permutation_test(ranks, n_j, k, B)
    p_val = mean(T_values >= T_obs)
    return(list(
      statistic = T_obs,
      p.value = p_val,
      method = "Exact Kruskal-Wallis Test"
    ))
  }
  
  if (method == "permutation") {
    # Perform Monte Carlo permutation test
    T_values = perform_permutation_test(ranks, n_j, k, B)
    p_val = mean(T_values >= T_obs)
    return(list(
      statistic = T_obs,
      p.value = p_val,
      method = "Permutation Kruskal-Wallis Test"
    ))
  }
}

data$group = factor(data$group)
set.seed(3005)
kw.test(data$x, data$group, method="asymptotic")$p.value
kw.test(data$x, data$group, method = "exact", B=100000)$p.value
kw.test(data$x, data$group, method = "permutation", B=100000)$p.value


#--------------------------------------------
# Question 4
#--------------------------------------------
data = read.csv("C:\\Users\\brian\\Documents\\CUHK\\24-25 Term1\\STAT3005\\Final Project\\Q3.csv")

loocv = function(x,y,p=1,K=NULL,B=11,from=NULL,to=NULL,plot=TRUE,...){
  # preliminary definitions 
  #-------------------------------------------------------------
  if(is.null(K)) K = dnorm
  n = length(x)
  X = array(1,dim=c(n,p+1))
  out = rep(NA,n)
  R = diff(range(x))
  if(is.null(from)) from=R/n^(4/5)
  if(is.null(to)) to=R/10
  bw.all = seq(from=from,to=to,length=B)
  MSE = rep(NA,length(bw.all))
  L = rep(NA,n)
  # LOOCV
  #-------------------------------------------------------------
  for(i.bw in 1:length(bw.all)){
    cat(i.bw," >> ")
    bw = bw.all[i.bw]
    for(i in 1:n){
      at = x[i]
      if(p>0){
        for(j in 1:p){
          X[,j+1] = (x-at)^j
        }
      }
      W = diag(K((x-at)/bw))
      XW = t(X)%*%W
      temp = solve(XW%*%X,XW)[1,]
      out[i] = sum(temp*y)
      L[i] = temp[i]
    }
    MSE[i.bw] = mean(((y-out)/(1-mean(L)))^2)
  }
  bw.opt = bw.all[which.min(MSE)]
  # plots
  #-------------------------------------------------------------
  if(plot){ 
    plot(bw.all, MSE, type="b", ylab="MSE", xlab="Bandwidth",...)
    abline(v=bw.opt, col="red3", lty=2)
  }
  # return optimal BW
  #-------------------------------------------------------------
  list(bw=bw.all, MSE=MSE, bw.opt=bw.opt)
}
lreg = function(x,y,p=1,K=NULL,at=NULL,bw=NULL,plot=TRUE,...){
  # preliminary definitions 
  #-------------------------------------------------------------
  if(is.null(K)) K = dnorm
  n = length(x)
  X = array(1,dim=c(n,p+1))
  if(is.null(at)){ 
    at = seq(min(x),to=max(x),length=301)
  }else{
    at = sort(at)
  }
  if(is.null(bw)) bw = (max(x)-min(x))/sqrt(n)
  n.at = length(at)
  at.all = at
  out = rep(NA,n.at)
  # compute the estimate
  #-------------------------------------------------------------
  for(i.at in 1:n.at){
    at = at.all[i.at]
    if(p>0){
      for(j in 1:p){
        X[,j+1] = (x-at)^j/factorial(j)
      }
    }
    # Method 1: more transparent
    W = diag(K((x-at)/bw))
    XW = t(X)%*%W
    out[i.at] = sum(solve(XW%*%X,XW)[1,]*y)
    # Method 2: faster and easier
    # out[i.at] = lm(y~X-1, weights=K((x-at)/bw))$coef[1]
  }
  # plot the data and regression line
  #-------------------------------------------------------------
  if(plot){
    optional = list(...)
    if(is.null(optional$type)) type="l"
    if(is.null(optional$col)) col="red4"
    if(is.null(optional$lwd)) lwd=3
    if(is.null(optional$lty)) lty=1
    I = order(at.all)
    plot(x,y, pch=19, col="grey", ...)
    points(at.all[I], out[I], type=type, col=col,lwd=lwd, lty=lty)
    
  }
  # return estimates 
  #-------------------------------------------------------------
  list(at=at.all,fit=out,bw=bw)
}

x = data$x
y = data$group
K = function(t) dnorm(t)
R = diff(range(x))
n = length(x)
at = seq(from=min(x), to=max(x), length=301)

par(mfrow=c(1,2))

out.cv = loocv(x,y,p=1,K=K, from=R/n^(2/3), to=R/7, plot=TRUE, main="LOOCV")
out.cv$bw.opt

plot(x,y,pch=19,col="grey",xlab="x", ylab="Group", main="Regression")
lines(at, lreg(data$x, data$group, p=1, K=K, at=at, bw=out.cv$bw.opt, plot=FALSE)$fit, col="steelblue3", lwd=3)
lines(at, lreg(data$x, data$group, plot=FALSE)$fit, bw=0.05, col="goldenrod2", lwd=3)

linear_reg = lm(y~x)$coef
abline(a=linear_reg[1], b=linear_reg[2], col="pink", lwd=3)
legend("bottomleft", 
       c("Local Polynomial Regression Estimate (optimal)",
         "Local Polynomial Regression Estimate (ad-hoc)",
         "Linear Regression"), 
       col=c("steelblue3", "goldenrod2", "pink"),
       lty=1, lwd=3, cex=0.45,
       bg="transparent"
)





