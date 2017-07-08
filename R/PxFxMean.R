#' \code{p()} functions for Edgeworth expansions terms
#' 
#' Calculate values of \code{p(x)} used in terms of Edgeworth expansions (EE) 
#' for a sample mean and a two-sample difference in means.
#' 
#' \code{p1()} is used for term 1 (second order expansion), \code{p2()} for term
#' 2 and so on.
#' 
#' @name pfuns
#' @param x a numerical value (quantile of sampling distribution) or a vector of
#'   quantiles.
#' @inheritParams qfuns
#'   
#' @return A vector of values of a mathematical function, one for each quantile 
#'   in \code{x}.
#' @seealso \code{\link{Fmean}}, \code{\link{Fdiff}} for corresponding EE
#'   values.
#'   
#' @examples 
#' 
#' @export
p1 <- function(x, lam) {
  lam3 <- lam[1]
  -1/6*He2(x)*lam3
}
#' @rdname pfuns
#' @export
p2 <- function(x, lam) {
  lam3 <- lam[1]
  lam4 <- lam[2]
  -1/72*He5(x)*lam3^2 - 1/24*He3(x)*lam4
}	
#' @rdname pfuns
#' @export
p3 <- function(x, lam) {
  lam3 <- lam[1]
  lam4 <- lam[2]
  lam5 <- lam[3]
  -1/1296*He8(x)*lam3^3 - 1/144*He6(x)*lam3*lam4 - 1/120*He4(x)*lam5
}
#' @rdname pfuns
#' @export
p4 <- function(x, lam) {
  lam3 <- lam[1]
  lam4 <- lam[2]
  lam5 <- lam[3]
  lam6 <- lam[4]
  -1/31104*He11(x)*lam3^4 - 1/1728*He9(x)*lam3^2*lam4 - 
    1/1152*He7(x)*lam4^2 - 1/720*He7(x)*lam3*lam5 - 1/720*He5(x)*lam6
}   

#' Edgeworth expansion for standardized mean
#' 
#' Calculate values of 1 - 4-term Edgeworth expansions (EE) (2nd - 5th order)
#' for the standardized mean.
#' 
#' Higher-order approximations of the cumulative distribution function of the 
#' standardized sample average.
#' 
#' @name Fmean
#' @inheritParams pfuns
#' @param n sample size.
#'   
#' @return  A vector of values of Edgeworth expansion of a corresponding order 
#'   (e.g. \code{Fm1} for a 1-term or 2nd order EE expansion and so on). The 
#'   length of the vector is the same as the length of \code{x}.
#'   
#' @seealso \code{\link{Fdiff}} for EE for two-sample standardized difference in
#'   means and \code{\link{pfuns}} for \code{p()} functions used in EE terms.
#'   
#' @examples 
#' 
#' @export
Fm1 <- function(x, n, lam) {
  pnorm(x)       + n^(-1/2)*p1(x, lam)*dnorm(x)
}
#' @rdname Fmean
#' @export
Fm2 <- function(x, n, lam) {
  Fm1(x, n, lam) + n^(-1)  *p2(x, lam)*dnorm(x)
}	
#' @rdname Fmean
#' @export
Fm3 <- function(x, n, lam) {
  Fm2(x, n, lam) + n^(-3/2)*p3(x, lam)*dnorm(x)
}	
#' @rdname Fmean
#' @export
Fm4 <- function(x, n, lam) {
  Fm3(x, n, lam) + n^(-2)  *p4(x, lam)*dnorm(x) 
}	

#' Edgeworth expansion two-sample difference in means
#' 
#' Calculate values of 1 - 4-term Edgeworth expansions (EE) (2nd - 5th order)
#' for the two-sample standardize difference in means.
#' 
#' Higher-order approximations of the cumulative distribution function of the 
#' two-sample standardized difference in means. Note that for a sample \eqn{X_1,
#' \dots, X_{n_x}, Y_1, \dots, Y_{n_y}}{X[1], ..., X[n[x]], Y[1], ..., Y[n[y]]},
#' \eqn{X} would correspond to treatment and \eqn{Y} to control.
#' 
#' @name Fdiff
#' @inheritParams pfuns
#' @param n sample size (one-sample).
#' @param nx number of observations in the first group (two-sample).
#' @param ny number of observations in the second group (two-sample).
#' @param lamx scaled cumulants of the distribution of the first group 
#'   (two-sample).
#' @param lamy scaled cumulants of the distribution of the second group 
#'   (two-sample).
#' @param varx variance of the first group (two-sample).
#' @param vary variance of the second group (two-sample).
#'   
#' @return  A vector of values of Edgeworth expansion of a corresponding order 
#'   (\code{Fm1} and \code{Fm1two} for a 1-term or 2nd order Edgeworth expansion
#'   and so on). The length of the vector is the same as the length of \code{x}.
#'   
#' @seealso \code{\link{Fmean}} for EE for one-sample standardized mean and 
#'   \code{\link{pfuns}} for \code{p()} functions used in EE terms.
#'   
#' @export
Fm1two <- function(x, nx, ny, lamx, lamy, varx, vary) {
  n <- (nx + ny)/2
  bx <- n/nx
  by <- n/ny
  j <- 3:6
  lam <- (lamx[j - 2]*bx^(j - 1)*varx^(j/2) + 
            (-1)^j*lamy[j - 2]*by^(j - 1)*vary^(j/2))/(bx*varx + by*vary)^(j/2)
  Fm1(x, n, lam)	        
}
#' @rdname Fdiff
#' @export
Fm2two <- function(x, nx, ny, lamx, lamy, varx, vary) {
  n <- (nx + ny)/2
  bx <- n/nx
  by <- n/ny
  j <- 3:6
  lam <- (lamx[j - 2]*bx^(j - 1)*varx^(j/2) + 
            (-1)^j*lamy[j - 2]*by^(j - 1)*vary^(j/2))/(bx*varx + by*vary)^(j/2)
  Fm2(x, n, lam)	        
}
#' @rdname Fdiff
#' @export
Fm3two <- function(x, nx, ny, lamx, lamy, varx, vary) {
  n <- (nx + ny)/2
  bx <- n/nx
  by <- n/ny
  j <- 3:6
  lam <- (lamx[j - 2]*bx^(j - 1)*varx^(j/2) + 
            (-1)^j*lamy[j - 2]*by^(j - 1)*vary^(j/2))/(bx*varx + by*vary)^(j/2)
  Fm3(x, n, lam)	        
}
#' @rdname Fdiff
#' @export
Fm4two <- function(x, nx, ny, lamx, lamy, varx, vary) {
  n <- (nx + ny)/2
  bx <- n/nx
  by <- n/ny
  j <- 3:6
  lam <- (lamx[j - 2]*bx^(j - 1)*varx^(j/2) + 
            (-1)^j*lamy[j - 2]*by^(j - 1)*vary^(j/2))/(bx*varx + by*vary)^(j/2)
  Fm4(x, n, lam)	        
}
