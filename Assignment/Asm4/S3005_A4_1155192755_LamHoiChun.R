###-------------------------------------------------------------###
### STAT 3005  
### Assignment 4  
###-------------------------------------------------------------###

setwd("C:\\Users\\brian\\Documents\\CUHK\\24-25 Term1\\STAT3005\\Asm4\\Images")


###-------------------------------------------------------------###
### Exercise 4.1.1
###-------------------------------------------------------------###
#install.packages("jpeg")
img = rep(list(NA),4)
img[[1]] = jpeg::readJPEG("D00.jpg")
img[[2]] = jpeg::readJPEG("D10.jpg")
img[[3]] = jpeg::readJPEG("D01.jpg")
img[[4]] = jpeg::readJPEG("D11.jpg")

# Plots of color images
par(mfrow=c(1,4),mar=c(0,0,1,0)) 
for(i in 1:4){ 
   plot(0:1,0:1,cex=0,axes=FALSE,main=c("D00","D10","D01","D11")[i])
   rasterImage(img[[i]],0,0,1,1) 
}

# Produce grayscale images 
get.gray = function(img, w=c(0.22,0.71,0.07)){
   img[,,1]*w[1] + img[,,2]*w[2] + img[,,3]*w[3]
}
img0 = rep(list(NA),4)
for(i in 1:4){
   img0[[i]] = get.gray(img[[i]])
}

# Plot of grayscale images
par(mfrow=c(1,4),mar=c(0,0,1,0)) 
for(i in 1:4){ 
   plot(0:1,0:1,cex=0,axes=FALSE,main=c("D00","D10","D01","D11")[i])
   rasterImage(img0[[i]], 0, 0, 1, 1)
}


###-------------------------------------------------------------###
### Exercise 4.1.2
###-------------------------------------------------------------###
col = c("royalblue","firebrick3","forestgreen","goldenrod2")
par(mfrow=c(1,1),mar=c(4.5,5,3,2)) 
for(j in 1:4){
   data = c(img0[[j]])
   plot(ecdf(data),
      col=col[j],add=j>1,ylab=expression(widehat(F)(x)),xlab=expression(x))
}
treatment = c("(a) Control","(b) EB 1089","(c) IR","(d) EB 1089 + IR")
legend("bottomright",treatment,col=col,lwd=1,bg="white", cex=0.65)

###-------------------------------------------------------------###
### Exercise 4.1.4
###-------------------------------------------------------------###
treatment_b = c(img0[[2]])
treatment_c = c(img0[[3]])
ks.test(treatment_b, treatment_c)$p.value # p-value < 2.2e-16

###-------------------------------------------------------------###
### Exercise 4.2
###-------------------------------------------------------------###
A = c(img0[[1]])
D = c(img0[[4]])
ks.test(A, D, alternative="greater")$p.value # p-value < 2.2e-16













###-------------------------------------------------------------###
### Remark: computation of one-sided two-sample KS test statisic
###-------------------------------------------------------------###
ks.test.statistic.2sided = function(x,y){
   m = length(x) 
   n = length(y)
   FY_at_X = ecdf(y)(sort(x))    # ecdf(y) is a function. We can compute, e.g., ecdf(y)(0.4)
   K_plus  = sqrt(n/(1+n/m))*max( (1:m)/m - FY_at_X     )
   K_minus = sqrt(n/(1+n/m))*max( FY_at_X - (0:(m-1))/m )
   K = max(K_plus, K_minus)
   c(two.sided=K, less=K_minus, greater=K_plus)
}

# simulated dataset
n = 5
m = 8
x = rnorm(n)
y = rnorm(m)

# self-written version
ks.test.statistic.2sided(A,D)

# built-in version
c(ks.test(A,D, alternative="two.sided")$statistic, 
  ks.test(A,D, alternative="less")$statistic, 
  ks.test(A,D, alternative="greater")$statistic) * sqrt(length(x)*length(y)/(length(x)+length(y)))


