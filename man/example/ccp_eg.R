## This is an example using simulated EATS data

## Parameter values
mux = -0.30 #true mean of X
su2 = 2.69 #true variance of U
sx2 = 1.01 #true variance of X
lambda = sx2 / (sx2 + su2) #attenuation
b = 1.54 #beta_1
a = -1.33 #beta_0

## Sample size
n = 629
k = 2 # Number of replicates in external data set
## Generate data set W_ij=x_i+u_ij
set.seed(20173)
x = rnorm(n, mux, sqrt(sx2))
u = matrix(rnorm(n * k, 0, sqrt(su2)), n, k)
ww = matrix(rep(x, k), n, k, byrow = FALSE) + u # Matrix of observed W with replicates

## Generate values of the variable y
fHm <- function(x){1 / (1 + exp(-(a + b * x)))}
pr = fHm(x)
y = vector()
for(i in 1:n){y[i] = rbinom(1, 1, pr[i])}

## Apply ccp for logistic model
ccp(y = y, W_int = ww, Type = "logistic")


