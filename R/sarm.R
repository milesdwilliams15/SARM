#' Estimate a Strategic Autoregressive Model
#'
#' This function takes a user supplied equation object and
#' dataframe. It assumes that the first two variables on the
#' right-hand side of the equation are the actor-specific
#' resource constraint and the spatial lag (sum of other-actor
#' expenditures on the strategic good), respectively. The
#' remaining variables are assumed to predict intrinsic
#' preference for the strategic good. NA values are allowable.
#' By default, it is assumed that estimated standard errors are not
#' robust, nor are they clustered according to some grouping variable.
#' If robust standard errors are desired, the user should set
#' robust = TRUE. If robust-clustered standard errors are desired,
#' the user should set both robust = TRUE and cluster = TRUE. If
#' clustered standard errors are desired, it is assumed that the last
#' variable included on the right-hand side of the equation is the
#' grouping variable. The grouping variable may be numeric, a factor,
#' or a character string.
#'
#' @export
sarm = function(eq, data = NULL, robust = F, cluster = F){

  # Set up data for analysis
  if(cluster == F){
    mat = model.frame(lm(eq, data = data))
    y = mat[,1]
    R = mat[,2]
    Y = mat[,3]
    Z = cbind(control = 1, mat[,4:ncol(mat)]) %>% as.matrix
  } else {
    mat = model.frame(lm(eq, data = data))
    y = mat[,1]
    R = mat[,2]
    Y = mat[,3]
    Z = as.matrix(cbind(control = 1, mat[,4:(ncol(mat) - 1)]))
    cluster_var = mat[,ncol(mat)]
  }

  # Write objective function
  mle = function(y, R, Y, Z, pars){
    y_hat = pnorm(Z%*%pars[1:ncol(Z)]) * R -
      (1 - pnorm(Z%*%pars[1:ncol(Z)])) * pars[ncol(Z)+1] * Y
    sigma = exp(pars[ncol(Z)+2])
    D = y == 0
    ll = sum(
      D*log(pnorm(-y_hat/sigma)) + (1 - D)*log(dnorm((y-y_hat)/sigma)/sigma)
    )
    return(-ll)
  }

  # Estimate model with either classic SEs or robust:
  if(robust == F){
    # Optimize mle
    opt = optim(
      fn = mle,
      y = y,
      R = R,
      Y = Y,
      Z = Z,
      par = rep(0, len = ncol(Z)+2),
      method = "L-BFGS-B",
      control = list(maxit = 1000),
      hessian = T
    )

    # Make summary output
    term = c(colnames(Z),"strategic parameter","log(scale)")
    estimate = opt$par
    std.error = sqrt(abs(diag(solve(opt$hessian))))
    statistic = estimate/std.error
    p.value = 2*pnorm(-abs(statistic))
  } else {
    # Optimize mle
    opt = optim(
      fn = mle,
      y = y,
      R = R,
      Y = Y,
      Z = Z,
      par = rep(0, len = ncol(Z)+2),
      method = "L-BFGS-B",
      control = list(maxit = 1000),
      hessian = T
    )

    # Get robust sandwich estimator
    pars = opt$par

    # Clustering?
    if(cluster == F){
      grads = list()
      for(i in 1:nrow(Z)) {
        grads[[i]] = numDeriv::grad(func = mle, y = y[i], R = R[i], Y = Y[i],
                          Z = matrix(Z[i,],nrow = 1, ncol = ncol(Z)),
                          x = opt$par)
      }
      grads = do.call(rbind, grads)
      dfa = (nrow(Z) - 1)/(nrow(Z) - length(opt$par))
    } else {
      grads = list()
      for(i in 1:length(unique(cluster_var))){
        clust = which(cluster_var == unique(cluster_var)[i])
        gradsi = numDeriv::grad(func = mle, y = y[clust], R = R[clust], Y = Y[clust],
                      Z = matrix(Z[clust,],nrow = length(clust), ncol = ncol(Z)),
                      x = opt$par)
        grads[[i]] = gradsi #apply(as.matrix(gradsi), 1, function(x) rep(x, len = length(clust)))
      }
      grads = do.call(rbind, grads)
      dfa = ((length(unique(cluster_var)))/(length(unique(cluster_var)) - 1)) *
        (nrow(Z) - 1)/(nrow(Z) - length(opt$par))
    }
    bread = solve(opt$hessian)
    meat = crossprod(grads)
    sandwich = dfa * (bread%*%meat%*%bread)

    # Make summary output
    term = c(colnames(Z),"strategic parameter","log(scale)")
    estimate = pars#opt$par[-length(opt$par)]
    std.error = sqrt(abs(diag(sandwich)))
    statistic = estimate/std.error
    p.value = 2*pnorm(-abs(statistic))
  }

  # Estimate predicted value and get adjusted R-squared
  pars = opt$par[-length(opt$par)]
  y_hat = pnorm(Z%*%pars[1:ncol(Z)]) * R -
    (1 - pnorm(Z%*%pars[1:ncol(Z)])) * pars[ncol(Z)+1] * Y
  r_sqr = cor(y[y>0],y_hat[y>0])^2 / mean((y == 0) == (y_hat <= 0))
  adj_r_sqr = 1 - ((1 - r_sqr)*(nrow(mat) - 1))/(nrow(mat) - (ncol(Z)+1) - 1)

  # Return the type of error used
  type = c("Standard","Robust","Cluster Robust")[
    c(robust == F, robust == T & cluster == F, robust == T & cluster == T)
    ]

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
              model_frame = if(cluster == F) mat else mat[,-ncol(mat)],
              model_eq = eq,
              r_sqr = r_sqr,
              adj_r_sqr = adj_r_sqr,
              residuals = y - y_hat,
              type = type))
}
