# Exercise 2.1
id = "1-RP9fZbyi6YIiYBi4g5zCoSGWGCA0-0L"
data = read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download",id))
head(data)
mass = data[,3]
duration = data[,4]
duration_unique_animal <- data[!duplicated(data$Animal), 4]

#=========================================================================
# Exercise 2.1 1(a)
# Method 1: Self function
n <- length(duration_unique_animal)
T_value <- (sum(sign(duration_unique_animal-21)))

Var_T <- n
t_stat <- ((T_value)-0)/sqrt(Var_T)
(p_value <- 2 * pnorm(-abs(t_stat)))
alpha_0 <- 0.05

if(p_value>alpha_0) {
  cat("Do not reject H_O since the p-value =", p_value, ">", alpha_0)
} else{
  cat("We reject H_O since the p-value =", p_value, , "<", alpha_0)
}

# Method 2: R function
binom.test((T_value+n)/2, n)$p.value

#=========================================================================
# Exercise 2.1 2(b)
# 2-sample t-test
sample1 <- subset(data, Mass.kg < 500)
sample2 <- subset(data, Mass.kg >= 500)

sample1_median <- median(sample1$Duration.s)
sample2_median <- median(sample2$Duration.s)

n1 <- nrow(sample1)
n2 <- nrow(sample2)

s1_square <- var(sample1$Duration.s)
s2_square <- var(sample2$Duration.s)
pooled_sd <- sqrt(((n1-1)*s1_square + (n2-1)*s2_square) / (n1+n2-2))

t_stat_2sample <- (sample1_median-sample2_median) / (pooled_sd*sqrt(1/n1 + 1/n2))

pnorm(t_stat_2sample)
pt(t_stat_2sample, df = 23, lower.tail = TRUE)
t.test(sample1$Duration.s, sample2$Duration.s, paired=FALSE, alternative = "less", var.equal = TRUE)$p.value

# Method 2: Wilcoxon rank sum test
wilcox.test(sample1$Duration.s, sample2$Duration.s, alternative = "less")$p.value

# QQ Plot
qqnorm(duration)
qqline(duration)

# ================================================================================= #
# Exercise 2.2.5

set.seed(123)

test.3sample = function(x1, x2, x3){ 
  x = c(x1, x2, x3)
  nj = c(length(x1), length(x2), length(x3))
  n = sum(nj)
  p = nj / sum(nj)
  R = rank(x)
  
  T = (sum(R[1:nj[1]]) - sum(R[(nj[1]+nj[2]+1):n])) / (n+1)
  muT = n/2 * (p[1] - p[3])
  varT = n^2 / 12 / (n+1)*(p[1]+p[3]-(p[1]-p[3])**2)
  t = (T - muT) / sqrt(varT)
  return(t)
} 

n1_values = c(10, 20, 30)
exact_quantiles = c(0.80, 0.90, 0.95, 0.99, 0.999)
nRep = 10000

run_simulation = function(n1) {
  n2 = 2*n1
  n3 = 3*n1
  t_values = numeric(nRep)
  
  for (i in 1:nRep) {
    x1 = 3005 + runif(n1, -2, 2)
    x2 = 3005 + runif(n2, -2, 2)
    x3 = 3005 + runif(n3, -2, 2)
    t_values[i] = test.3sample(x1, x2, x3)
  }

  return(quantile(t_values, exact_quantiles))
}

for (n1 in n1_values) {
  cat("\nQuantiles for n1 =", n1, ":\n")
  quantiles = run_simulation(n1)
  print(quantiles)
}  

set.seed(123)
asymptotic_quantiles = qnorm(exact_quantiles)
cat("Asymptotic quantiles (from N(0,1)):\n")
print(asymptotic_quantiles)