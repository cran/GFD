#' @export 
plot.GFD <- function (x, ...) {
  
  object <- x
  dots <- list(...)
  a <- object$plotting
  b <- object$Descriptive
  fac.names <- a$fac_names
  exist <- hasArg(factor) 
   
  if(length(fac.names) != 1){
    if(!exist){
      print("Please choose the factor you wish to plot (for interaction type something like group1:group2) and confirm by pressing 'Enter'")
      Faktor <- scan("", what="character")
      while(length(Faktor)==0){
        print("Please enter the name of the factor you wish to plot!")
        Faktor <- scan("", what="character")
      }
      } else {
       Faktor <- dots$factor
      } 
    x.label <- ""
  } else {
    Faktor <- fac.names
    x.label <- fac.names
  }
    
  match.arg(Faktor, fac.names)

  # default values
  args <- list(plot.object = a, descr.object = b, factor = Faktor,
               lwd =2, ylab = "Means", xlab = x.label, col = 1:length(fac.names), pch = 1:18, legendpos = "topright")
  
  args[names(dots)] <- dots
  
  do.call(plotting, args = args)
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


#' @export
print.qanova <- function(x, ...) {
  if(length(x$input$p) == 1){
    cat("\n","Variance estimation:", x$input$var_method, "\n", "Test is based on", "the", paste0(paste0(100*x$input$p, "%-",  collapse = ", "), "Quantile."), "\n")
  }
  if((length(x$input$p) > 1) && ( identical(x$lin_mat, diag(rep(1, length(x$input$p))) ) )){
    cat("\n","Variance estimation:", x$input$var_method, "\n", "Test is based on", "the", paste0(paste0(100*x$input$p, "%-",  collapse = ", "), "Quantiles"), "simultaneously.", "\n")
  }
  if( (length(x$input$p) > 1) && !( identical(x$lin_mat, diag(rep(1, length(x$input$p))) ) )){
    cat("\n","Variance estimation:",x$input$var_method, "\n", "Test is based on the linear combination of", "the","\n", paste0(paste0(100*x$input$p, "%-",  collapse = ", "), "Quantiles"), "by the matrix", "\n")
    print(x$input$lin_mat)
  }
  cat("\n","Call:", "\n")
  print(x$input$formula)
  
  cat("\n", "QANOVA: Quantile-based Analysis of Variance:","\n", sep = "")
  print(x$statistics)
}

#' @export
summary.qanova <- function (object, ...) {
  print(object)
}
