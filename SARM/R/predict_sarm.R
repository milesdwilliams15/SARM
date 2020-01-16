#' Get Predicted Values from a SARM Model
#'
#' This function reports predicted values with 95% confidence intervals
#' for a model estimated with the sarm() function. By default, it
#' generates predicted values with the data used to estimate the model.
#' The user may create or use a different dataset, and then include
#' this using the newdata command.
#'
#' @export
predict_sarm = function(model, newdata  = NULL){
  if(is.null(newdata)){
    mat = model$model_frame
    Z = as.matrix(cbind(control = 1, mat[,4:ncol(mat)]))
  } else {
    mat = newdata
    Z = as.matrix(cbind(control = 1, mat[,4:ncol(mat)]))
  }

  # Estimate forcasted values
  y = mat[,1]
  R = mat[,2]
  Y = mat[,3]
  pars = model$summary$estimate
  y_hat = pnorm(Z%*%pars[1:ncol(Z)]) * R -
    (1 - pnorm(Z%*%pars[1:ncol(Z)])) * pars[ncol(Z)+1] * Y

  # Bootstrap forcasted values
  boot_y_hat = matrix(0, nrow = nrow(mat), ncol = 999)
  for(i in 1:999) {
    boot_pars = pars +
      rnorm(n = length(pars), sd = model$summary$std.error[-length(pars)])
    boot_y_hat[,i] = pnorm(Z%*%boot_pars[1:ncol(Z)]) * R -
      (1 - pnorm(Z%*%boot_pars[1:ncol(Z)])) * boot_pars[ncol(Z)+1] * Y
  }
  se_fit = apply(boot_y_hat, 1, sd)
  preds = tibble::tibble(
    fit = y_hat,
    lwr = y_hat - 1.96*se_fit,
    upr = y_hat + 1.96*se_fit
  )
  return(preds)
}
