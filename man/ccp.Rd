% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccp.R
\name{ccp}
\alias{ccp}
\title{Fit misspecified model with measurement error}
\usage{
ccp(y, W_int, W_ext = NULL, C = NULL, Type, print.summary = TRUE, standardize = TRUE)
}
\arguments{
\item{y}{A vector of the response variable. It can be binary for logistic regression, or continuous for linear regression.}

\item{W_int}{A \eqn{n \times R} matrix of covariate \eqn{W} (main data set), where \eqn{W_{ir} = X_i + U_{ir}}. \eqn{R} is the replicates per each observation. \eqn{X} is the true but unobserved value and \eqn{U} is iid measurement error which is independent of \eqn{X}. If \code{W_{int}} is a vector (i.e., no replicates), where \eqn{W_i = X_i + U_i}, \code{W_{ext}} is needed.}

\item{W_ext}{A \eqn{N \times K} matrix of covariates \eqn{W} from the external data set, where \eqn{W_{ik} = X_i + U_{ik}}. Note that when \code{W_{int}} is a matrix, \code{W_{ext}} will be ignored and the estimation will be only based on \code{W_{ext}}. The default setting is "NULL".}

\item{C}{A vector indicating 4 cut points, usually representing the quartiles of \eqn{X}. The default setting is "NULL", and it will be automatically calculated based on nuisance parameters.}

\item{Type}{Model type, either "logistic" or "linear".}

\item{print.summary}{Print a summary of all estimates. The default setting is \code{TRUE}.}

\item{standardize}{if standardization needs to be performed. The default is \code{TRUE}.}
}
\value{
A list of
    \item{theta5-theta1}{Estimate of \eqn{\theta_5-\theta_1}, interpreted as log relative risk in logistic regression.}
    \item{theta}{Estimates of \eqn{\Theta = (\theta_1, ..., \theta_5)}.}
    \item{nuisance}{Estimates of nuisance parameters (\eqn{\mu_x, \sigma_x^2, \sigma_u^2}) as well as parameters in the true model (\eqn{\beta_0, \beta_1}).}
    \item{se.theta}{Standard errors of \eqn{\Theta}.}
    \item{se.nuisance}{Standard errors of nuisance parameters and parameters in the true model.}
}
\description{
This function fits categorical model to a categorized continuous variable with measurement error. Though the methodology does not require specific models, we provide linear regression and logistic regression in this package. If the main data set has no replicates, an external data set with replicates is needed. Further, the nuisance parameters \eqn{(\mu_x, \sigma^2_x, \sigma^2_u)} should be transportable (see Chapter 2.2.4 and 2.2.5 of Carroll et al. (2006)).
}
\details{
{
Let \eqn{Y} be the 0/1 dependent variable, and \eqn{X} be the continuous predictor subject to measurement error. The true model is \eqn{\hbox{Pr}(Y = 1|X)=H(\beta_0+X\beta_1)}, where \eqn{H(x)=\exp(x)/\{1+\exp(x)\}} is the logistic distribution function. Let \eqn{W} be the observed variable with measurement error, \eqn{W=X+U}. \eqn{U \sim N(0,\sigma^2_u)} is the measurement error and \eqn{X \sim N(\mu_x,\sigma^2_x)} is unobserved true value. The misspecified model and the asymptotic theory for the parameters are derived on the manuscript "Categorizing a Continuous Predictor Subject to Measurement Error". If \eqn{I(X\in C_j)} means that \eqn{X} is in category \eqn{j = 1, ..., 5,} the categorical (but incorrect) model is \eqn{\hbox{Pr}(Y=1|X)=H\{\sum_{j=1}^5\theta_jI(X\in C_j)\}}.

When \eqn{Y} is a continuous variable, the true model and misspecified categorical model will be linear regression correspondingly. Other assumptions for \eqn{X, U} and \eqn{W} remain the same as described above.}
}
\section{References}{

Betsabe Blas, Tianying Wang, Victor Kipnis, Kevin Dodd and Raymond Carroll, "Categorizing a Continuous Predictor Subject to Measurement Error" (2018+).
}

\examples{
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


}
