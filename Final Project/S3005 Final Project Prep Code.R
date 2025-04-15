#-----------------------------------------------------------------------------------
# Ch5: Distribution Problem
# Simulation for Critical Values: Example 5.8 
# V1: 3 Tests: KS test, CVM test, AD test
#-----------------------------------------------------------------------------------
n.all = c(50,100,200)
test.all = c("ks","cvm","ad")
nRep = 2^16
if(1){
  out = array(NA,dim=c(nRep,length(n.all),length(test.all)))
  dimnames(out) = list(paste0("iRep=",1:nRep), paste0("n=",n.all), test.all)
  for(iRep in 1:nRep){
    set.seed(iRep)
    x.all = runif(max(n.all))
    for(i.n in 1:length(n.all)){
      n = n.all[i.n]
      x = x.all[1:n]
      out[iRep,i.n,1] = ks.test0(x, punif, exact=FALSE)$statistic
      out[iRep,i.n,2] = cvm.test0(x, punif, exact=FALSE)$statistic
      out[iRep,i.n,3] = ad.test0(x, punif, exact=FALSE)$statistic
    }
  }
}
alpha = c(0.01, 0.05, 0.1)
c.value = apply(out, 2:3, quantile, 1-alpha)
c.value

#-----------------------------------------------------------------------------------
# Ch5: Distribution Problem
# Simulation for Critical Values
# V2: 1 Test
#-----------------------------------------------------------------------------------
n.all = c(50, 100, 200)
nRep = 2^16
alpha = c(0.01, 0.05, 0.1)

out = matrix(NA, nrow = nRep, ncol = length(n.all))
colnames(out) = paste0("n=", n.all)

for (iRep in 1:nRep) {
  set.seed(iRep)
  x.all = runif(max(n.all))
  for (i.n in 1:length(n.all)) {
    n = n.all[i.n]
    x = x.all[1:n]
    test_result = ks.test(x, "punif", exact = FALSE)
    out[iRep, i.n] = test_result$statistic
  }
  if (iRep %% 10000 == 0) {
    cat("Completed replication", iRep, "out of", nRep, "\n")
  }
}

c.value = apply(out, 2, quantile, probs = 1 - alpha)
rownames(c.value) = paste0("alpha=", alpha)
c.value

#-----------------------------------------------------------------------------------
# Ch5: Lillifor's Normality Test
# H0: F = pnorm(.,mu,sigma) vs H1: not H0
#-----------------------------------------------------------------------------------
l.test0 = function(x, exact=TRUE){
  n = length(x)
  i = 1:n
  muHat = mean(x)
  sdHat = sqrt(mean(x^2)-muHat^2)
  F0 = function(x) pnorm(x,muHat,sdHat)
  u = F0(sort(x))
  Kpos = max(i/n-u)
  Kneg = max(u-(i-1)/n)
  L = sqrt(n)*max(Kpos,Kneg)
  if(exact){
    nRep = 10000
    out = rep(NA,nRep)
    for(iRep in 1:nRep) out[iRep] = l.test0(rnorm(n), exact=FALSE)$statistic
    p = mean(out>L)
   } else{
    p = NA
   }
   list(statistic=L, p.value=p)
}
nortest::lillie.test(x)
#-----------------------------------------------------------------------------------
# Ch7: Density Problem
# Deriving the Normal Reference Rule
#-----------------------------------------------------------------------------------
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
ddnorm(1.2,0,1,0)
ddnorm(1.2,0,1,1)
ddnorm(1.2,0,1,4)

p = integrate(function(x){ddnorm(x,0,1,2)^2},-10,10)$value
((1/(2*sqrt(pi)))/4/(1/4*p))^(1/5)
#-----------------------------------------------------------------------------------
# Ch9: Non-parametric Regression
# Local Polynomial Regression Estimator from scratch
# Example 9.7
# BW can be self defined
#-----------------------------------------------------------------------------------
# Function
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
# Application
set.seed(1)
n = 400
m = function(x) sin(20*x)*x    # regression function
x = (1:n)/n                    # covariate 
y = m(x)+rnorm(n)/5            # response  
data = data.frame(x,y)
head(data)
K = function(t) dnorm(t)
at = seq(from=min(x), to=max(x), length=301)
out = lreg(x, y, p=1, K=K, at=at, bw=0.05, plot=FALSE)
head(data.frame(at=out$at, fit=out$fit))
plot(x,y,col="goldenrod2",pch=19)
lines(out$at, out$fit, lwd=3,col="royalblue4")


#-----------------------------------------------------------------------------------
# Ch9: Non-parametric Regression
# Leave-one-out Cross Validation (LOOCV)
#-----------------------------------------------------------------------------------
# Function
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
# Application
set.seed(1)
n = 400
m = function(x) sin(20*x)*x    # regression function
x = (1:n)/n                    # covariate 
y = m(x)+rnorm(n)/5            # response  

# LOOCV
#---------------------------------------------------------------
par(mfrow=c(1,1))
p = 1
K = dnorm
out.cv = loocv(x,y,p=1,K=K,plot=TRUE) # BW Selection
out.cv$bw.opt
at = seq(from=min(x), to=max(x), length=301)
out = lreg(x, y, p=p, K=K, at=at, bw=out.cv$bw.opt, plot=TRUE) # Regression Plot

# BW Selection Plot
plot(out.cv$bw, out.cv$MSE, type="b", ylab="MSE", xlab="Bandwidth",
     main="Bandwidth Selection by LOOCV")
abline(v=out.cv$bw.opt, col="red3", lty=2)

# Regression Plot w/ Truth
col = c("black","royalblue4","red3","goldenrod2","pink")
pch = c(NA,1,3,0,4)
lty = c(1,1,1,1,2)
plot(x,y, pch=19,cex=0.5, xlim=range(at),col="grey",
     main="Local polynomial regression estimation",
     xlab=expression(italic(x)), ylab=expression(italic(y)))
at = seq(from=min(x), to=max(x), length=301)
points(at,m(at),col=col[1],lwd=4,type="l")
K = function(t) dnorm(t)
out = lreg(x,y,p=1,K=K,at=at,bw=out.cv$bw.opt,plot=FALSE)
points(out$at, out$fit, type="l", pch=pch[1+2],
       col=col[1+2], lwd=4, lty=lty[1+2])
out = lreg(x,y,p=1,K=K,at=at,bw=0.05,plot=FALSE)
points(out$at, out$fit, type="l", pch=pch[2+2],
       col=col[2+2], lwd=4, lty=lty[2+2])
linear_model = lm(y~x)$coef
abline(a=linear_model[1], b=linear_model[2], col="pink", lty=2, lwd=3)
legend("topleft", c("Truth","Estimate (optimal)","Estimate (ad-hoc)", "Parametric"), col=col[c(1,3,4,5)], 
       lty=lty[c(1,3,4,2)], lwd=4, cex=0.5,bg="white")


