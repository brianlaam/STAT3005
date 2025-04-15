###-------------------------------------------------------------###
### Exercise 3.2 (1)
###-------------------------------------------------------------###

# Step 1: Download data
stock = c("0001.hk")
from  = as.Date("2018-01-02")        
to    = as.Date("2024-10-10")
data  = tseries::get.hist.quote(stock, from, to, quote="Adjusted",quiet=TRUE)

# Step 2: Compute the rate of return 
#---------------------------------------------------------------
s = c(as.matrix(data))  
s = na.omit(s)
n = length(s)
r = head(diff(s)/lag(s), -1)

# Step 3: Check if tie exists
#---------------------------------------------------------------
sum(duplicated(r))            # <<<< Complete this line

# Step 4: Plot the prices and returns 
#---------------------------------------------------------------
par(mfrow=c(1,2))
ts.plot(s, main="Stock Price", xlab=expression(i),  ylab=expression(S[i]), col="blue")
ts.plot(r, main="Rate of Return", xlab=expression(i),  ylab=expression(r[i]), col="blue")

# Step 5: Produce autocorrelation plot 
#---------------------------------------------------------------
par(mfrow=c(1,2))
acf(r, main=expression("ACF of {"*r[i]*"}"), col="red")
acf(r^2, main=expression("ACF of {"*r[i]^2*"}"), col="red")





###-------------------------------------------------------------###
### Exercise 3.2 (2)
###-------------------------------------------------------------###

# Step 1: define x and y
#---------------------------------------------------------------
x = r[-1]
y = r[-1664]

# Step 2: define space for storing the result table 
#---------------------------------------------------------------
result = array(NA, dim=c(3,5))
rownames(result) = c("Estimate", "p-value", "time")
colnames(result) = c("P","S","K","BD","C")

# Step 3: (a) estimate of correlation 
#---------------------------------------------------------------
for(i in 1:3){  
	result[1,i] = cor(x,y,method=c("p","s","k")[i])
}
result[1,4] = cov.bd0(x,y)/(cov.bd0(x,x)*cov.bd0(y,y))         # <<<< Compute BD correlation (NOT covariance); see chp4.R
result[1,5] = cor.c1(x,y)                                      # <<<< Compute Chatterjee's correlation; see chp4.R
  
# Step 4: (b-c) p-value and time taken 
#---------------------------------------------------------------
for(i in 1:5){
	t0 = Sys.time()
	result[2,i] = cor.test0(x,y,method=c("p","s","k","b","c")[i])     # <<<<  cor.test0 is a function in chp4.R.
	t1 = Sys.time()
	result[3,i] = difftime(t1, t0, units="secs")   # <<<< Complete this line
}
result

###--------------------------------------------------------------------###
#                     P           S          K           BD           C
# Estimate 0.0092922365 0.030320416 0.02127512 0.0005463601 0.045604737
# p-value  0.7049405865 0.216526275 0.19394768 0.2210000000 0.003174619
# time     0.0008981228 0.001039028 0.02346706 8.8962469101 0.036216021
###--------------------------------------------------------------------###




###-------------------------------------------------------------###
### Exercise 3.2 (4)
###-------------------------------------------------------------###
# Step 1: (a) 
install.packages("tseries")
residuals = tseries::garch(r, order=c(1,1), trace=0)$residuals[-1] # r = rate of return 
ts.plot(residuals, main="Residuals Plot", xlab=expression(i),  ylab=expression(epsilon[i]), col="blue")

# Step 2: (b)
residuals_x = r[-1]
residuals_y = r[-length(r)]
cor.test.c1(residuals_x,residuals_y) # p-value = 0.0113211

# Rep function Version 1
rep.cor.test.c0 <- function(x,y, nRep){
  p_value <- 0
  for (iRep in 1:nRep){
    p_value <- cor.test.c0(x,y) + p_value
  }
  p_value/nRep 
}
rep.cor.test.c0(residuals_x, residuals_y, 10000) # p-value = 0.003862024

# Rep function Version 2
rep2.cor.test.c0 <- function(x,y,nRep){
  p_value <- numeric(0)
  for (iRep in 1:nRep){
    p_value <- append(p_value, cor.test.c0(x,y))
  }
  mean(p_value)
}
rep2.cor.test.c0(residuals_x, residuals_y, 10000) # p-value = 0.01403572

# cor.test.c0(residuals_x,residuals_y) # p-value = 0.002586751
# We reject H_0 since p-value < alpha = 0.05




###-------------------------------------------------------------###
### Exercise 3.3 (5)
###-------------------------------------------------------------###
n.all = c(5, 50 ,500)
test_method = c("p", "s", "k", "c")
out = array(dim = c(length(n.all), length(test_method)),
            dimnames = list(paste0("n = ", n.all), test_method))

for(i in 1:length(test_method)){
  for(i.all in 1:length(n.all)){
    A = rexp(n.all[i.all], 1)
    X = A^3
    Y = exp(A)
    out[i.all, i] = cor0(X, Y, method=test_method[i])
  }
}
out

