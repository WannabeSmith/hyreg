# Density of tobit model
dtobit <- function(y, xTbeta, sigma, left = -1, log = FALSE)
{
  d <- y > left

  if (log)
  {
    return(d * (dnorm((y - xTbeta)/sigma, log = TRUE) - log(sigma)) +
             (1 - d) * pnorm((left - xTbeta)/sigma, log.p = TRUE))
  } else
  {
    return((1/sigma * dnorm((y - xTbeta)/sigma))^d * pnorm((xTbeta - left)/sigma, lower.tail = FALSE)^(1-d))
  }
}

# log-likelihood of tobit regression model
llTobit <- function(y, xTbeta, sigma, left = -1)
{
  sum(dtobit(y, xTbeta, sigma, left = left, log = TRUE))
}

# log-likelihood of logistic regression model
llLogisticReg <- function(y, xTbeta)
{
  sum(-(1 - y) * xTbeta - log(1 + exp(-xTbeta)))
}

# log-likelihood of the hybrid tobit-logit regression model
llHybrid <- function(y.tobit, y.discrete, sigma, xTbeta.tobit, xTbeta.discrete)
{
  llTobit(y = y.tobit, xTbeta = xTbeta.tobit, sigma = sigma, left = -1) +
    llLogisticReg(y = y.discrete, xTbeta = xTbeta.discrete)
}

#' Perform hybrid tobit-logit regression
#'
#' This function gets estimates of beta, sigma, and theta in the hybrid tobit-logit
#' model. That is, when the logit betas are a scalar multiple of the tobit betas.
#'
#' @importFrom stats model.matrix dnorm pnorm optim
#' @importFrom survival survreg
#' @param formula.tobit a regression formula describing the relationship between the tobit response and the covariates
#' @param formula.discrete a regression formula describing the relationship between the bernoulli response and the covariates
#' @param data.tobit the data.frame containing the tobit responses and covariates
#' @param data.discrete the data.frame containing the bernoulli responses and covariates
#' @param start a vector of starting values for optimization. If not specified, starting values are taken from a non-hybrid tobit model, and theta is set to 1.
#' @param left a number specifying where left-censoring occurred
#' @return a list containing the following parameter estimates (from maximum likelihood):\cr
#' \item{beta}{the regression coefficients}
#' \item{sigma}{the standard deviation of the censored normal distribution}
#' \item{theta}{the multiplicative factor relating the two sets of regression coefficients}
#' @export
hyreg <- function(formula.tobit, formula.discrete, data.tobit, data.discrete, start = NULL, left = -1)
{
  response.varname.tobit <- all.vars(formula.tobit)[1]
  y.tobit <- as.vector(data.tobit[, response.varname.tobit])

  response.varname.discrete <- all.vars(formula.discrete)[1]
  y.discrete <- as.vector(data.discrete[, response.varname.discrete])

  X.tobit <- model.matrix(formula.tobit, data.tobit)
  X.discrete <- model.matrix(formula.discrete, data.discrete)

  stopifnot(ncol(X.tobit) == ncol(X.discrete))

  # Number of betas
  num.betas <- ncol(X.tobit)

  # Number of parameters is num.betas plus 2. One for sigma (std dev),
  #   one for theta (multiplicative factor between discrete and tobit betas)
  num.params <- num.betas + 2

  # Create starting values if not provided
  if(is.null(start))
  {
    model.start <- survreg(Surv(y.tobit, y.tobit > left, type = "left") ~
                             X.tobit[, -1], dist = "gaussian")
    start <- c(model.start$coefficients, model.start$scale, 1)
  }

  names(start) <- c(colnames(X.tobit), "sigma", "theta")

  # Specify lower bounds for L-BFGS-B
  # all betas can take on values from -Inf to Inf, but sigma and theta can only take on values > 0.
  #   However, setting the lower bound exactly to 0 causes numerical issues in the optimization procedure.
  lower = c(rep(-Inf, num.betas), 10e-16, 10e-16)

  # Create the objective function (this is created here so that num.betas,
  #   X.tobit, etc. are all in scope)
  objective <- function(params)
  {
    beta.tobit <- params[1:num.betas]
    sigma <- params[num.betas + 1]
    theta <- params[num.betas + 2]

    xTbeta.tobit <- X.tobit %*% beta.tobit
    xTbeta.discrete <- theta * X.discrete %*% beta.tobit

    return(-llHybrid(y.tobit = y.tobit, y.discrete = y.discrete, xTbeta.tobit = xTbeta.tobit,
                     xTbeta.discrete = xTbeta.discrete, sigma = sigma))
  }

  optimum <- optim(par = start, fn = objective,
                   lower = lower, method = "L-BFGS-B")

  beta.ests <- optimum$par[1:num.betas]
  sigma.est <- optimum$par[num.betas + 1]
  theta.est <- optimum$par[num.betas + 2]

  return(list("beta" = beta.ests,
              "sigma" = sigma.est,
              "theta" = theta.est))
}

