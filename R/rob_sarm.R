#' Estimate a Strategic Autoregressive Model with Modeled Variance
#'
#' This function is an extension of the sarm() model. It allows the user
#' to estimate the sarm() model but with model variance given as a function
#' of a grouping variable. This approach may be adopted if the user
#' suspects the data violate the assumptions of constant variance and
#' independence of observations within groups. When estimating this model,
#' it is assumed that the last variable included on the right-hand side
#' of the equation object is the grouping variable.
#'
#' @export
rob_sarm = function(eq, data = NULL){

  # Set up data for analysis

  mat = model.frame(lm(eq, data = data))
  y = mat[,1]
  R = mat[,2]
  Y = mat[,3]
  Z = as.matrix(cbind(control = 1, mat[,4:(ncol(mat) - 1)]))
  cluster_var = as.factor(mat[,ncol(mat)])

  # Write objective function
  mle = function(y, R, Y, Z, cluster_var, pars){
    y_hat = pnorm(Z%*%pars[1:ncol(Z)]) * R -
      (1 - pnorm(Z%*%pars[1:ncol(Z)])) * pars[ncol(Z)+1] * Y
    dy_hat = ave(y, cluster_var, FUN = var)
    dy_hat = rank(dy_hat)
    dy_hat = (dy_hat - mean(dy_hat))/sd(dy_hat)
    dy_mat = as.matrix(cbind(1, dy_hat))
    sigma = exp(dy_mat%*%pars[ncol(Z)+1+(1:ncol(dy_mat))])
    D = y == 0
    ll = sum(
      D*log(pnorm(-y_hat/sigma)) + (1 - D)*log(dnorm((y-y_hat)/sigma)/sigma)
    )
    return(-ll)
  }

  # Estimate model with either classic SEs or robust:
  # Optimize mle
  opt = optim(
    fn = mle,
    y = y,
    R = R,
    Y = Y,
    Z = Z,
    cluster_var = cluster_var,
    par = rep(0, len = ncol(Z)+3),
    method = "L-BFGS-B",
    control = list(maxit = 1000),
    hessian = T
  )

  # Make summary output
  term = c(colnames(Z),"strategic parameter")
  estimate = opt$par[1:(ncol(Z)+1)]
  std.error = sqrt(abs(diag(solve(opt$hessian))))[1:(ncol(Z)+1)]
  statistic = estimate/std.error
  p.value = 2*pnorm(-abs(statistic))

  # Estimate predicted value and get adjusted R-squared
  pars = opt$par[1:(ncol(Z)+1)]
  y_hat = pnorm(Z%*%pars[1:ncol(Z)]) * R -
    (1 - pnorm(Z%*%pars[1:ncol(Z)])) * pars[ncol(Z)+1] * Y
  r_sqr = cor(y[y>0],y_hat[y>0])^2 / mean((y == 0) == (y_hat <= 0))
  adj_r_sqr = 1 - ((1 - r_sqr)*(nrow(mat) - 1))/(nrow(mat) - (ncol(Z)+1) - 1)

  # Make a tidy table summarizing results
  tidy_summary = tibble::tibble(
    term = term,
    estimate = estimate,
    std.error = std.error,
    statistic = statistic,
    p.value = p.value
  )
  return(list(summary = tidy_summary,
              log_like = opt$value,
              model_frame = mat[,-ncol(mat)],
              model_eq = eq,
              r_sqr = r_sqr,
              adj_r_sqr = adj_r_sqr,
              residuals = y - y_hat,
              type = "Robust"))
}
