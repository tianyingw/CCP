#' Fit misspecified model with measurement error
#'
#' @description  This function fits categorical model to a categorized continuous variable with measurement error. Though the methodology does not require specific models, we provide linear regression and logistic regression in this package. If the main data set has no replicates, an external data set with replicates is needed. Further, the nuisance parameters \eqn{(\mu_x, \sigma^2_x, \sigma^2_u)} should be transportable (see Chapter 2.2.4 and 2.2.5 of Carroll et al. (2006)).
#'
#' @usage ccp(y, W_int, W_ext = NULL, C = NULL, Type, print.summary = TRUE, standardize = TRUE)
#'
#' @param y A vector of the response variable. It can be binary for logistic regression, or continuous for linear regression.
#' @param W_int A \eqn{n \times R} matrix of covariate \eqn{W} (main data set), where \eqn{W_{ir} = X_i + U_{ir}}. \eqn{R} is the replicates per each observation. \eqn{X} is the true but unobserved value and \eqn{U} is iid measurement error which is independent of \eqn{X}. If \code{W_{int}} is a vector (i.e., no replicates), where \eqn{W_i = X_i + U_i}, \code{W_{ext}} is needed.
#' @param W_ext A \eqn{N \times K} matrix of covariates \eqn{W} from the external data set, where \eqn{W_{ik} = X_i + U_{ik}}. Note that when \code{W_{int}} is a matrix, \code{W_{ext}} will be ignored and the estimation will be only based on \code{W_{ext}}. The default setting is "NULL".
#' @param C A vector indicating 4 cut points, usually representing the quartiles of \eqn{X}. The default setting is "NULL", and it will be automatically calculated based on nuisance parameters.
#' @param Type Model type, either "logistic" or "linear".
#' @param print.summary Print a summary of all estimates. The default setting is \code{TRUE}.
#' @param standardize if standardization needs to be performed. The default is \code{TRUE}.
#'
#' @details{
#' Let \eqn{Y} be the 0/1 dependent variable, and \eqn{X} be the continuous predictor subject to measurement error. The true model is \eqn{\hbox{Pr}(Y = 1|X)=H(\beta_0+X\beta_1)}, where \eqn{H(x)=\exp(x)/\{1+\exp(x)\}} is the logistic distribution function. Let \eqn{W} be the observed variable with measurement error, \eqn{W=X+U}. \eqn{U \sim N(0,\sigma^2_u)} is the measurement error and \eqn{X \sim N(\mu_x,\sigma^2_x)} is unobserved true value. The misspecified model and the asymptotic theory for the parameters are derived on the manuscript "Categorizing a Continuous Predictor Subject to Measurement Error". If \eqn{I(X\in C_j)} means that \eqn{X} is in category \eqn{j = 1, ..., 5,} the categorical (but incorrect) model is \eqn{\hbox{Pr}(Y=1|X)=H\{\sum_{j=1}^5\theta_jI(X\in C_j)\}}.
#'
#' When \eqn{Y} is a continuous variable, the true model and misspecified categorical model will be linear regression correspondingly. Other assumptions for \eqn{X, U} and \eqn{W} remain the same as described above.}
#'
#' @return A list of
#'     \item{theta5-theta1}{Estimate of \eqn{\theta_5-\theta_1}, interpreted as log relative risk in logistic regression.}
#'     \item{theta}{Estimates of \eqn{\Theta = (\theta_1, ..., \theta_5)}.}
#'     \item{nuisance}{Estimates of nuisance parameters (\eqn{\mu_x, \sigma_x^2, \sigma_u^2}) as well as parameters in the true model (\eqn{\beta_0, \beta_1}).}
#'     \item{se.theta}{Standard errors of \eqn{\Theta}.}
#'     \item{se.nuisance}{Standard errors of nuisance parameters and parameters in the true model.}
#'
#' @section References:
#' Betsabe Blas, Tianying Wang, Victor Kipnis, Kevin Dodd and Raymond Carroll, "Categorizing a Continuous Predictor Subject to Measurement Error" (2018+).
#' @example man/example/ccp_eg.R
#'
#' @import stats MASS numDeriv
#'
#' @export
#'
ccp <- function(y, W_int, W_ext = NULL, C = NULL, Type, print.summary = TRUE,standardize = TRUE){

  #### make sure W_int is a matrix (user may provide vectors or lists)
  W_int = as.matrix(W_int)
  #### check data and compatibility
  if (Type != "logistic" & Type != "linear"){
    stop("Please choose type, either logistic or linear!")
  }
  #####check compatibility############
  if (any(is.na(y)) | any(is.na(W_int)))
    stop("Missing values found in W or Y! Please check!")
  #### check case

  #### check the dimension
  if ((dim(W_int))[1] != length(y)) {
      stop("Number of observations in Y doesn't match the number of rows in W! Please check!")
  }

  if(dim(W_int)[2] > 1){
    #### if W_int has replicates, we consider Allint case
    if (Type == "logistic"){
    out = Allint_logistic(y = y, W_int = W_int, C = C, print.summary = print.summary, standardize = standardize)}
    if (Type == "linear"){
      out = Allint_linear(y = y, W_int = W_int, C = C, print.summary = print.summary, standardize = standardize)
    }
    }else{
      #### check if external data is provided
      if(is.null(W_ext)){stop("External data need to be provided!")}

      #### make sure the external data is in matrix form
      W_ext = as.matrix(W_ext)

      #### external data must have replicates
      if(dim(W_ext)[2] == 1){stop("External data need to have replicates!")}

      if (Type == "logistic"){
      out = Extint_logistic(y = y, W_int = W_int, W_ext = W_ext, C = C, print.summary = print.summary,standardize = standardize)}
      if (Type == "linear"){
        out = Extint_linear(y = y, W_int = W_int, W_ext = W_ext, C = C, print.summary = print.summary, standardize = standardize)
      }
    }
  return(out)
}
