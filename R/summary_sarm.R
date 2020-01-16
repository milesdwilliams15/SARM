#' Report Model Summary Statistics for SARM
#'
#' This function returns summary output from a model estimated with
#' sarm(). It acts much like the summary.lm() function.
#'
#' @export
summary_sarm = function(model) {
  stars.pval = function(pval) {
    p = 0
    for(i in 1:length(pval)){
      if(pval[i] <= 0.001){
        p[i] = "***"
      } else if(pval[i] <= 0.01) {
        p[i] = "**"
      } else if(pval[i] <= 0.05) {
        p[i] = "*"
      } else if(pval[i] <= 0.1) {
        p[i] = "."
      } else {
        p[i] = " "
      }
    }
    return(p)
  }
  cat("RESULTS FOR STRATEGIC AUTOREGRESSIVE MODEL\n")
  cat("\n")
  cat("Standard Errors:",model$type,"\n")
  cat("\n")
  cat("Call:\n")
  print(model$model_eq)
  cat("\n")
  cat("Residuals:\n")
  print(quantile(model$residuals))
  cat("\n")
  cat("Coefficients:\n")
  summ = data.frame(dplyr::mutate(model$summary, x = stars.pval(p.value)))
  colnames(summ) = c("","estimate","std.error","statistic","p.value","")
  print(summ, row.names = F, right = F)
  cat("---\n")
  cat("sig. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n")
  cat("\n")
  cat("Summary Statistics:\n")
  cat(paste("RSE:", round(sd(model$residuals),3),
            "on", nrow(model$model_frame) - nrow(model$summary) + 1, "df\n"))
  cat(paste("Log-Likelihood:", round(model$log_like,3),"\n"))
  cat(paste("R-squared:", round(model$r_sqr,3),"\n"))
  cat(paste("Adj. R-squared:", round(model$adj_r_sqr,3)))
}
