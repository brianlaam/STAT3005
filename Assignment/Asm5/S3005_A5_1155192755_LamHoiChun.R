###-------------------------------------------------------------###
### STAT 3005  
### Assignment 5
###-------------------------------------------------------------###

# =============================================================================
# Q5.2(e) & (f)
# Using the Order Permutation Test Function for Paired Sample From Tut09
order.ptest = 
  function(x, y, FUN=NULL, B=1e4, alternative=c("two.sided","less","greater"), plot=FALSE){
    # FUN         = function of test statistic 
    # mu0         = null value of mu
    # B           = number of Monte Carlo replications
    # alternative = direction of alternative hypothesis (e.g., "less" means H1: mu(x)<mu(y))
    # plot        = logical indicating whether to plot permutation distribution
    alternative = match.arg(alternative)
    if(is.null(FUN)) FUN = function(x,y) cov(x,y)/cov(x,x)
    n = length(x)
    Tobs = FUN(x,y)
    Tb = rep(NA,B)
    for(b in 1:B){
      xb = x[sample(1:n,n)]
      Tb[b] = FUN(xb,y)
    }
    if(plot){
      hist(Tb,nclass=30,freq=FALSE,xlim=range(c(Tb,Tobs)),xlab=expression(italic(T)),
           main="Permutation distribution",
           col="pink",border=FALSE,cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
      abline(v=Tobs,col="blue4",lty=2,lwd=2)
    }
    switch(alternative, "two.sided"=mean(abs(Tb)>=abs(Tobs)),
           "less"=mean(Tb<=Tobs),
           "greater"=mean(Tb>=Tobs))
  }

# Data
data = read.csv("C:\\Users\\brian\\Documents\\CUHK\\24-25 Term1\\STAT3005\\Asm5\\clouds.csv")
A = data$A
X = data$X

# OLS estimator
FUN = function(x,y){
  n = length(x)
  sum((x[1:n]*y[1:n])/0.3 - ((1-x[1:n])*y[1:n])/0.7)
}

# Permutation p-value
set.seed(1)
order.ptest(A,X, FUN, alternative="greater", plot=TRUE) # p-value = 7e-04


# =============================================================================
# Q5.3(e) & (f)
# OLS estimator
FUN2 = function(x,y){
  n = length(x)
  sum((x[1:n]*y[1:n])/0.3 - 5 * ((1-x[1:n])*y[1:n])/0.7)
}

# Permutation p-value
set.seed(1)
order.ptest(A,X, FUN2, alternative="greater", plot=TRUE) # p-value = 7e-04

