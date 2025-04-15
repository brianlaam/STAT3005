# Exercise 1.1 Q2
# Set parameters
nRep = 2^10
n = 30
out = array(NA, dim=c(nRep, 3, 3))
dimnames(out) = list(paste0("iRep=", 1:nRep),
                     paste0("DGP", c("(i)", "(ii)", "(iii)")),
                     paste0("estimator", 1:3))

# True values of sigma^2 for each DGP
estimand = c(5^2, 1/3, 2)

# Simulation
for(iRep in 1:nRep) {
  for(iDGP in 1:3) {
    if(iDGP==1) x = rnorm(n, 10, 5)
    if(iDGP==2) x = runif(n, -1, 1)
    if(iDGP==3) x = 100 * ((1:n)/n - 0.5)^2 + rt(n, df=4)
    
    # Estimator 1
    out[iRep, iDGP, 1] = mean(x^2)  
    
    # Estimator 2
    x_sorted = sort(x)
    q75 = x_sorted[floor(0.75 * n)]
    q25 = x_sorted[floor(0.25 * n)]
    out[iRep, iDGP, 2] = ((q75 - q25) / (2 * qnorm(0.75)))^2
    
    # Estimator 3
    differences = diff(x)  # Calculate differences between consecutive order statistics
    out[iRep, iDGP, 3] = (1 / (2 * (n - 1))) * sum(differences^2)
    
    # Log differences
    out[iRep, iDGP, ] = log(out[iRep, iDGP, ]) - log(estimand[iDGP])
  }
}

# Calculate RMSE
rmse_results = sqrt(apply(out^2, 2:3, mean))
print(rmse_results)

rm(list = ls())