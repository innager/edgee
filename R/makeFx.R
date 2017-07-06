#' Make Edgeworth expansions
#' 
#' Create an Edgeworth series \code{F(x)} of orders 1 - 5 for various versions 
#' of t-statistic.
#' 
#' Given sample statistics or distribution parameters, produce Edgeworth 
#' expansions (EE) for t-statistic as a function of \code{x} (quantile). 
#' Includes short (for ordinary one-sample t-statistic) and general versions, 
#' which can be used for ordinary and moderated one- and two-sample t-statistics
#' and Welch t-test. Variance adjustment \code{r} is incorporated and there is 
#' not need to pass it to the resulting function.
#' 
#' @inheritParams tailDiag
#' @inheritParams makeQx
#' @param type type of Edgeworth expansions to be used: \code{"short"} for 
#'   ordinary one-sample t-statistic, \code{"one-sample"} for a more complicated
#'   version such as a moderated t-statistic, and \code{"two-sample"} for any 
#'   kind of two-sample t-statistic.
#' @param base base distribution of Edgeworth expansions. Classic expansions are
#'   based on standard normal distribution (\code{base = "normal"}), but 
#'   Student's t-distribution (\code{base = "t"}) is recommended for 
#'   approximating far tails of the sampling distribution.
#' @param df degrees of freedom for Student's t-distribution if \code{base = 
#'   "t"}. Provide a single value for the first order approximation (zero term).
#'   If not provided, \code{df} will be calculated as \code{n - 1} for 
#'   \code{type = "short" and type = "one-sample"} or \code{2n - 2} for 
#'   \code{type = "two-sample"}; for augmented degrees of freedom (\code{moder =
#'   TRUE}), the value of prior degrees of freedom \code{stats["d0"]} will be 
#'   added.
#'   
#' @return A function \code{F(x)} that takes a quantile \code{x} as an input and
#'   outputs a vector of five values for five orders of approximation. ***** 
#'   Check if \code{x} can be a vector (might be a problem with conditionals).
#'   
#' @seealso \code{\link{smpStats}} that creates \code{stats} vector from the
#'   sample, and \code{\link{makeQx}} that creates \code{q(x)} functions used in
#'   EE terms.
#'   
#' @examples 
#' 
#' @export

makeFx <- function(stats, n, r = NULL, type = "short", base = "normal", df = NULL,
                   moder = FALSE, verbose = TRUE) {
  q <- makeQx(stats, r = r, type = type, verbose = verbose)
  if (is.null(r)) {
    if (type == "short") {
      r <- sqrt((n - 1)/n)
    } else {
      r <- q$r
    }
  }
  if (base == "t") {
    if (is.null(df)) {
      if (type %in% c("one-sample", "short")) df = n - 1
      else                                    df = 2*n - 2
      if (moder) df <- df + stats["d0"]
    }
    FEdge <- function(x) {
      term0 <- pt(x, df)
      term1 <- term0 + 1/sqrt(n)*q[[1]](x)*dt(x, df + 2)
      term2 <- term1 + 1/n      *q[[2]](x)*dt(x, df + 5)
      term3 <- term2 + 1/n^(3/2)*q[[3]](x)*dt(x, df + 8)
      term4 <- term3 + 1/n^2    *q[[4]](x)*dt(x, df + 11)
      return(c(term0, term1, term2, term3, term4))
    }
  } else {
    FEdge <- function(x) {
      term0 <- pnorm(x)
      term1 <- term0 + 1/sqrt(n)*q[[1]](x)*dnorm(x)
      term2 <- term1 + 1/n      *q[[2]](x)*dnorm(x)
      term3 <- term2 + 1/n^(3/2)*q[[3]](x)*dnorm(x)
      term4 <- term3 + 1/n^2    *q[[4]](x)*dnorm(x)
      return(c(term0, term1, term2, term3, term4))
    }
  }
  return(function(x) FEdge(x/r))
}                  

