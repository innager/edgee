#' \code{p()} functions for Edgeworth expansions terms
#' 
#' Calculate values of \code{p(x)} used in terms of Edgeworth expansions (EE)
#' for a sample mean and a two-sample difference in means.
#' 
#' \code{p1()} is used for term 1 (second order expansion), \code{p2()} for term 2 and so on.
#' #' 
#' @param x a numerical value (quantile of sampling distribution) or a vector of quantiles.
#' @inheritParams q1s
#'   
#' @return A vector of values of a mathematical function, one for each quantile in \code{x}.
#' @seealso \code{\link{Fm1}}, \code{\link{Fm1two}} for corresponding EE values.
#'   
#' @examples 
#' 
#' @export
p1 <- function(x, lam) {
  lam3 <- lam[1]
  -1/6*He2(x)*lam3
}
#' @rdname p1
#' @export
p2 <- function(x, lam) {
  lam3 <- lam[1]
  lam4 <- lam[2]
  -1/72*He5(x)*lam3^2 - 1/24*He3(x)*lam4
}	
#' @rdname p1
#' @export
p3 <- function(x, lam) {
  lam3 <- lam[1]
  lam4 <- lam[2]
  lam5 <- lam[3]
  -1/1296*He8(x)*lam3^3 - 1/144*He6(x)*lam3*lam4 - 1/120*He4(x)*lam5
}
#' @rdname p1
#' @export
p4 <- function(x, lam) {
  lam3 <- lam[1]
  lam4 <- lam[2]
  lam5 <- lam[3]
  lam6 <- lam[4]
  -1/31104*He11(x)*lam3^4 - 1/1728*He9(x)*lam3^2*lam4 - 
    1/1152*He7(x)*lam4^2 - 1/720*He7(x)*lam3*lam5 - 1/720*He5(x)*lam6
}   

#' Edgeworth expansion functions for sample means
#' 
#' Calculate values of 1 - 4-term Edgeworth expansions (2 - 5 order) for a sample mean and a two-sample difference in means.
#' 
#' \code{Fm1, ..., Fm4} implement Edgeworth expansions for a sample average and \code{Fm1two, ..., Fm4two} - for a two-sample difference in means. Note that for a two-sample difference in means, the first sample (\eqn{\bar{X}}) would correspond to treatment while the second (\code{Y}) -  to control.
#' 
#' @inheritParams p1
#' @param n sample size (one-sample).
#' @param nx number of observations in the first group (two-sample).
#' @param ny number of observations in the second group (two-sample).
#' @param lamx scaled cumulants of the distribution of the first group (two-sample).
#' @param lamy scaled cumulants of the distribution of the second group (two-sample).
#' @param varx variance of the first group (two-sample). 
#' @param vary variance of the second group (two-sample). 
#'   
#' @return  A vector of values of Edgeworth expansion of a corresponding order (\code{Fm1}
#'   and \code{Fm1two} for a 1-term or 2nd order Edgeworth expansion and so on). The length of the vector is the same as the length of \code{x}.
#'   
#' @seealso \code{\link{p1}} for \code{p()} functions used in Edgeworth expansion terms.
#'   
#' @examples 
#' 
#' @export
Fm1 <- function(x, n, lam) {
  pnorm(x)       + n^(-1/2)*p1(x, lam)*dnorm(x)
}
#' @rdname Fm1
#' @export
Fm2 <- function(x, n, lam) {
  Fm1(x, n, lam) + n^(-1)  *p2(x, lam)*dnorm(x)
}	
#' @rdname Fm1
#' @export
Fm3 <- function(x, n, lam) {
  Fm2(x, n, lam) + n^(-3/2)*p3(x, lam)*dnorm(x)
}	
#' @rdname Fm1
#' @export
Fm4 <- function(x, n, lam) {
  Fm3(x, n, lam) + n^(-2)  *p4(x, lam)*dnorm(x) 
}	
#' @rdname Fm1
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
#' @rdname Fm1
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
#' @rdname Fm1
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
#' @rdname Fm1
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
