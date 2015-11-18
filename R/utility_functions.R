#' @export 
plot.GFD <- function (x, title = NULL, lwd = 2, ...) {
  a <- x$input
  new_object <- GFD(formula = a$formula, data = a$data,
                    nperm = 1, alpha = a$alpha,
                    plot_CI = TRUE, Title = title, line_width = lwd)
}

#' @export
print.GFD <- function(x, ...) {
  a <- x$input
  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Wald-Type Statistic (WTS):", "\n", sep = "")
  print(x$WTS)
  cat("\n", "ANOVA-Type Statistic (ATS):", "\n", sep = "")
  print(x$ATS)
}

#' @export
summary.GFD <- function (object, ...) {
  a <- object$input
  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Descriptive:", "\n", sep = "")
  print(object$Descriptive)
  cat("\n", "Wald-Type Statistic (WTS):", "\n", sep = "")
  print(object$WTS)
  cat("\n", "ANOVA-Type Statistic (ATS):", "\n", sep = "")
  print(object$ATS)
}
