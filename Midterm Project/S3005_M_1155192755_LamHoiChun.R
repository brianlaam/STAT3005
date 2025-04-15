###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###
### Stat 3005 -- Midterm (2024)
###
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------

#-------------------------------------
# Question 1
#-------------------------------------
scaleTrendTest = function(x){
	r = rank(x)
	n = length(r)
	T = sum((1:n)*abs(r/(n+1)-1/2))
	ET = (n^2-1)/8 
	varT = (n-1)*(n^2+3)/576
	t = (T-ET)/sqrt(varT) 
	t
}
# Under H0
n = 100
x = rnorm(n)
scaleTrendTest(x)

# Under H1
x = rnorm(n)*(n:1)/n 
scaleTrendTest(x)


num_samples = 100

# 100 test statistics under H0
test_statistics_H0 = numeric(num_samples)
for (i in 1:num_samples) {
  set.seed(i)
  x = rnorm(n)
  test_statistics_H0[i] = scaleTrendTest(x)
}

# 100 test statistics under H1
test_statistics_H1 = numeric(num_samples)
for (i in 1:num_samples) {
  set.seed(i)
  x = rnorm(n) * (n:1) / n
  test_statistics_H1[i] = scaleTrendTest(x)
}

hist_H0 = hist(test_statistics_H0, plot = FALSE, breaks = 20)
hist_H1 = hist(test_statistics_H1, plot = FALSE, breaks = 20)

# Plots
plot(hist_H0, col = "steelblue3", xlim = range(c(hist_H0$breaks, hist_H1$breaks)),
     ylim = range(c(hist_H0$counts, hist_H1$counts)), main = "Histogram of Test Statistics",
     xlab = "Test Statistic (t)", ylab = "Frequency", border = "black")
plot(hist_H1, col = "goldenrod2", add = TRUE, border = "black")
legend("topright", legend = c("H0", "H1"), fill = c("steelblue3", "goldenrod2"),
       border = c("black", "black"))

#-------------------------------------
# Question 2
#-------------------------------------
# download data 
id = "1bhOLQ8QzJu6P_dM7UHYK1-6zFLntaAgV" 
data = read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", id))
# I = !is.na(data[,3])

# Data Cleaning (Delete incomplete data)
data = data[8:30,]
head(data)
men   = data[,"Men"]
women = data[,"Women"]
table(men>women)
t.test(men[18:30],women[18:30])
ts.plot(data[,2:3])

# Self function for rank-sum test
append(men, women)
rank(append(men, women))
T = sum(head(rank(append(men, women)), n=length(men)))/(length(append(men, women))+1)
ET = length(men)/2
VarT = length(men)*length(women)/(12*(length(append(men, women))+1))
t = (T-ET)/VarT
p_val = 1-pnorm(t) # p-value = 0.567547

# Wilcoxon rank-sum test
(test_result = wilcox.test(men, women, alternative = "greater", exact = FALSE)) # p-value = 0.5698

# Normality of the data
par(mfrow = c(1,2))
qqnorm(men, main="QQ Plot for Men's Data")
qqline(men, col="mediumvioletred")
qqnorm(women, main="QQ Plot for Women's Data")
qqline(women, col="mediumvioletred")

#-------------------------------------
# Question 3
#-------------------------------------

# a simulated dataset in case 1 with n=50 and Delta=1.5
n = 50
Delta = 1.5
g1 = function(u,Delta=0){
	Delta*u^3
}
mu = g1((1:n)/(n+1), Delta=1.5) 
x = mu + rt(n,3)
ts.plot(x) 
points(1:n, mu, type="l", col="red")


#-------------------------------------
# Question 3.3
#-------------------------------------

# Function to compute hat{rho_C^(m)}
compute_rho_C_m = function(X, m) {
  n = length(X)
  R = rank(X)
  lambda = (n+1)*(n-m)*m
  sum_min_R = 0
  for (j in 1:m) {
    index_i = (1+m):n
    index_j = (1+m-j):(n-j)
    min_R = pmin(R[index_i], R[index_j])
    sum_min_R = sum_min_R + sum(min_R)
  }
  rho_C_m = -2 + (6/lambda)*sum_min_R
  rho_C_m
}

# Function to compute p-value for hat{rho_C^(m)}
compute_p_value = function(rho_C_m, m, n) {
  z_value = rho_C_m / sqrt(2/(5*m*n) + (8*m) / (15*n^2))
  p_value = 1 - pnorm(z_value)
  p_value
}

# Function to generate X_i for 3 cases
generate_X = function(n, Delta, case_num) {
  u = (1:n) / (n+1)
  epsilon = rt(n, df=3)
  if (case_num == 1) {
    mu = Delta*u^3
  } else if (case_num == 2) {
    mu = Delta*8*(u-0.5)^2
  } else if (case_num == 3) {
    mu = Delta*(2*u^2*cos(6*pi*u))
  }
  X = mu + epsilon
  X
}

# Function to compute p-values for all tests
compute_p_values = function(X, m_values) {
  n = length(X)
  p_values = numeric(length(m_values) + 1)
  names(p_values) = c(paste0("rho_C", m_values), "monotonic_trend")
  
  for (i in seq_along(m_values)) {
    m = m_values[i]
    rho_C_m = compute_rho_C_m(X, m)
    p_values[i] = compute_p_value(rho_C_m, m, n)
  }
  
  # Using Spearman Correlation Test to perform Monotonic Trend Test (time index against data)
  time_index = 1:n
  p_value_s = cor.test(time_index, X, method = "spearman", exact = FALSE)$p.value
  p_values["monotonic_trend"] = p_value_s
  
  p_values
}

# Simulation parameters
n = 200
nRep = 500
Delta_all = seq(0, 1.5, by = 0.1)
size = 0.05
cases = 1:3
methods = c("rho_C^(1)", "rho_C^(2)", "rho_C^(3)", "monotonic_trend")
m_values = 1:3
power_results = array(NA, dim = c(length(Delta_all), length(methods), length(cases)),
                       dimnames = list(paste0("Delta=", Delta_all), methods, paste0("Case", cases)))

# Run simulations
for (case_num in cases) {
  cat("Simulating for Case", case_num, "\n")
  for (Delta_idx in seq_along(Delta_all)) {
    Delta = Delta_all[Delta_idx]
    rejections = matrix(0, nrow = nRep, ncol = length(methods))
    colnames(rejections) = methods
    for (rep in 1:nRep) {
      set.seed(rep)
      X = generate_X(n, Delta, case_num)
      p_values = compute_p_values(X, m_values)
      rejections[rep, ] = (p_values < size)
    }
    power = colMeans(rejections)
    power_results[Delta_idx, , paste0("Case", case_num)] = power
  }
}

# Plot power curves
par(mfrow = c(1, length(cases)), mar = c(4.5, 5, 3, 2))
col = c("royalblue4", "green4", "goldenrod2", "mediumvioletred")
pch = c(1, 2, 8, 4)
lwd = rep(1.5, 4)
lty = 1

for (case_num in cases) {
  matplot(Delta_all, power_results[, , paste0("Case", case_num)], type = "b",
          lwd = lwd, lty = lty, col = col, pch = pch,
          ylim = c(0, 1.1), xlab = expression(Delta), ylab = "Power",
          main = paste("Power Curve for Case", case_num), cex.lab = 1.5, cex.main = 1.5)
  abline(h = c(0, size, 1), lty = 2, col = "gray")
  abline(v = 0, lty = 2, col = "gray")
  if (case_num == 1) {
    legend("topleft",
           lwd = 2, lty = lty, col = col, pch = pch, cex = 0.885, bty="n",
           legend = c(expression(hat(rho)[C]^{(1)}),
                      expression(hat(rho)[C]^{(2)}),
                      expression(hat(rho)[C]^{(3)}),
                      "Monotonic Trend Test"), 
           )
  }
}

#-------------------------------------
# Question 4
#-------------------------------------
par(mfrow = c(1, 1))
# Define the scaleTrendTest() function if not defined in Q1
scaleTrendTest = function(x){
  r = rank(x)
  n = length(r)
  T = sum((1:n)*abs(r/(n+1)-1/2))
  ET = (n^2-1)/8 
  varT = (n-1)*(n^2+3)/576
  t = (T-ET)/sqrt(varT) 
  t
}

# Simulation parameters
n <- 100
num_simulations <- 1000
alpha <- 0.05
delta_values <- seq(0, 3, by = 0.1)

# Storing Power
power_kimberley <- numeric(length(delta_values))
power_spearman <- numeric(length(delta_values))

# Run simulations
for (k in seq_along(delta_values)) {
  delta <- delta_values[k]
  reject_kimberley <- 0
  reject_spearman <- 0
  
  for (iRep in 1:num_simulations) {
    set.seed(iRep)
    if(iRep%%500==0) cat(">> ")
    
    sigma_i <- 1 + delta*(1-(1:n)/n)
    x <- rnorm(n, mean = 0, sd = sigma_i)
    
    # Kimberley's Test
    t_kimberley <- scaleTrendTest(x)
    p_value_kimberley <- 2 * (1 - pnorm(abs(t_kimberley)))
    if (p_value_kimberley < alpha) {
      reject_kimberley <- reject_kimberley + 1
    }
    
    # Spearman Rank Correlation Test
    abs_deviation <- abs(x - mean(x))
    spearman_test <- cor.test(1:n, abs_deviation, method = "s")
    if (spearman_test$p.value < alpha) {
      reject_spearman <- reject_spearman + 1
    }
  }
  
  # Calculate power for each test
  power_kimberley[k] <- reject_kimberley / num_simulations
  power_spearman[k] <- reject_spearman / num_simulations
}

# Plot the power curves
matplot(delta_values, cbind(power_kimberley, power_spearman), type = "b",
        lty = 1, col = c("royalblue4", "goldenrod2"), lwd = 2, pch = c(1, 8),
        xlab = expression(delta), ylab = "Power", ylim = c(0, 1.1),
        main = "Power Curves of Kimberley's Test and Spearman Correlation Test")
abline(h = c(0, alpha, 1), lty = 2, col = "gray")
abline(v = 0, lty = 2, col = "gray")
legend("topleft", legend = c("Kimberley's Test", "Spearman Test"),
       col = c("royalblue4", "goldenrod2"), lty = 1, lwd = 2, pch = c(1, 8), cex = 0.75, bty = "n")
