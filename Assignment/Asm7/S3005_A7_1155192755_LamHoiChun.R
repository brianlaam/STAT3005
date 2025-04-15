#-----------------------------------------------------------------------------------
# STAT3005 Assignment 7
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# KDE (self-written)
#-----------------------------------------------------------------------------------
density0 = function(x, kernel=NULL, bw=NULL, xeval=NULL, plot=FALSE, C=c(0.9,1.06), ...){
  n = length(x) 
  # define and standardize a kernel
  if(is.null(kernel)){
    kernel = dnorm
  }else{
    kernel0 = kernel
    sigma = sqrt(integrate(function(t){kernel0(t)*t^2},-Inf,Inf)$value)
    kernel = function(t){ kernel0(t*sigma)*sigma }
  }
  # rule-of-thumb bandwidth (See Section 8.3.2)
  if(is.null(bw)){ 
    s = min(sd(x), diff(quantile(x,c(.25,.75)))/1.34)
    bw = C[1]*s*n^(-1/5)
  }
  # set xeval as the values at which KDE are evaluated, i.e., KDE(xeval) 
  if(is.null(xeval)){ 
    xeval=seq(min(x)-0.2*diff(range(x)),to=max(x)+0.2*diff(range(x)),length=301)
  }
  # compute KDE
  d = rep(NA,length(xeval))
  for(j in 1:length(xeval)){
    d[j] = mean(kernel((x-(xeval[j]))/bw))/bw
  }
  # plot the KDE 
  if(plot){ 
    plot(xeval,d,type="l",ylim=c(0,max(d)),ylab="Density",
         xlab=bquote(N==.(n)~"  Bandwidth"==.(round(bw,5))),...)
  }
  # output 
  list(xeval=xeval, density=d)
}

#-----------------------------------------------------------------------------------
# Data
#-----------------------------------------------------------------------------------
data = read.csv("C:\\Users\\brian\\Documents\\CUHK\\24-25 Term1\\STAT3005\\Asm5\\clouds.csv")
A = data$A
X = data$X

not_seeded <- data$X[data$A == 0]
seeded <- data$X[data$A == 1]

#-----------------------------------------------------------------------------------
# Q4(a)
#-----------------------------------------------------------------------------------
tilde_f0 = density0(not_seeded, C=0.9, plot=TRUE)$density
tilde_f1 = density0(seeded, C=0.9, plot=TRUE)$density

par(mfrow=c(1,2))
plot(tilde_f0, xlab = "Rainfall (acre-feet)", ylab = "Density", col = "steelblue3",
     type="l", lwd="1.75", main=expression(tilde(f[0])))
plot(tilde_f1, xlab = "Rainfall (acre-feet)", ylab = "Density", col = "steelblue3",
     type="l", lwd="1.75",main=expression(tilde(f[1])))

#-----------------------------------------------------------------------------------
# Q4(b)
#-----------------------------------------------------------------------------------
hat_f0 = density0(not_seeded, bw=bw.SJ(not_seeded), plot=TRUE)$density
hat_f1 = density0(seeded, bw=bw.SJ(seeded), plot=TRUE)$density

plot(hat_f0, xlab = "Rainfall (acre-feet)", ylab = "Density", col = "steelblue3",
     type="l", lwd="1.75", main=expression(hat(f[0])))
plot(hat_f1, xlab = "Rainfall (acre-feet)", ylab = "Density", col = "steelblue3",
     type="l", lwd="1.75",main=expression(hat(f[1])))

#-----------------------------------------------------------------------------------
# Q4(d)
#-----------------------------------------------------------------------------------
par(mfrow=c(1,1))
hat_f = density0(X, bw=bw.SJ(X), plot=TRUE)$density
plot(hat_f, xlab = "Rainfall (acre-feet)", ylab = "Density", col = "steelblue3",
     type="l", lwd="1.75", main=expression(hat(f)))

#-----------------------------------------------------------------------------------
# Q4(f)
#-----------------------------------------------------------------------------------
hat_pi <- mean(data$A)
bar_f <- hat_pi * hat_f1 + (1 - hat_pi) * hat_f0


plot(bar_f, xlab = "Rainfall (acre-feet)", ylab = "Density", col = "steelblue3", 
     type="l", lwd=1.75, main = expression(hat(f)~and~f))
lines(hat_f, col = "goldenrod2", type="l", lwd=2)
legend("topright", legend = c("f", expression(hat(f))), col = c("steelblue3","goldenrod2"), lwd=2)











