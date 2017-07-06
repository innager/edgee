#' \code{q()} functions for Edgeworth expansions terms
#' 
#' Calculate values of \code{q(x)} used in terms of Edgeworth expansions (EE)
#' for short (ordinary one-sample) and general versions of t-statistic.
#' 
#' \code{q1s, ..., q4s} implement a short version of EE that can be used for
#' ordinary one-sample t-statistic. \code{q1s()} is used for term 1 (second
#' order expansion), \code{q2s()} for term 2 and so on. They do not include
#' variance adjustment \code{r}, which is equal to 1 for a t-statistic with
#' naive biased variance estimate. For a t with unbiased variance estimate,
#' which is a standard, use \code{q1s(x/r)}, etc. \code{q1k, ..., q4k} implement
#' a general version, which can be used for any t-statistic; variance adjustment
#' is incorporated in these functions.
#' 
#' @param x numerical value, a quantile of sampling distribution. *** a vector
#'   here but what about F?
#' @param lam a vector of scaled cumulants or their estimates, starting with 3rd
#'   (skewness). Term 1 needs a third cumulant only, term 2 - third and fourth
#'   (kurtosis), and so on.
#' @param r sqare root of variance adjustment. The variance adjustment is
#'   generally equal to \code{k21}.
#' @param k12,k13,k22,k23,k31,k32,k41,k42,k51,k61 cumulant components - values
#'   calculated from sample statistics or distribution parameters.
#'   
#' @return The value of a mathematical function. (or values, one for each x?)
#' @seealso \code{\link{Ft1}}, \code{\link{Ft1gen}} for corresponding EE values.
#'   
#' @examples 
#' 
#' @export
q1s <- function(x, lam) {
  lam3 <- lam[1]
  1/6*(2*x^2 + 1)*lam3
}
#' @rdname q1s
#' @export
q2s <- function(x, lam) {
  lam3 <- lam[1]
  lam4 <- lam[2]
  -1/18*(x^5 + 2*x^3 - 3*x)*lam3^2 - 1/4*x^3 + 1/12*(x^3 - 3*x)*lam4 - 3/4*x
}	
#' @rdname q1s
#' @export
q3s <- function(x, lam) {
  lam3 <- lam[1]
  lam4 <- lam[2]
  lam5 <- lam[3]
  1/1296*(8*x^8 + 28*x^6 - 210*x^4 - 525*x^2 - 105)*lam3^3 -
    1/144*(4*x^6 - 30*x^4 - 90*x^2 - 15)*lam3*lam4 + 1/24*(2*x^6 -                                           3*x^4 - 6*x^2)*lam3 - 1/40*(2*x^4 + 8*x^2 + 1)*lam5
}    
#' @rdname q1s
#' @export
q4s <- function(x, lam) {
  lam3 <- lam[1]
  lam4 <- lam[2]
  lam5 <- lam[3]
  lam6 <- lam[4]
  -1/32*x^7 - 1/1944*(x^11 + 5*x^9 - 90*x^7 - 450*x^5 + 45*x^3 + 945*x)*lam3^4 - 
    5/96*x^5 - 1/72*(x^9 - 6*x^7 - 12*x^5 - 18*x^3 -                                                  9*x)*lam3^2 - 1/288*(x^7 - 21*x^5 + 33*x^3 + 111*x)*lam4^2 +
    1/60*(x^7 + 8*x^5 - 5*x^3 - 30*x)*lam3*lam5 - 7/96*x^3 +
    1/432*(9*x^7 - 63*x^5 + 2*(x^9 - 12*x^7 - 90*x^5 + 36*x^3 +
           261*x)*lam3^2 + 81*x^3 + 189*x)*lam4 - 1/90*(2*x^5 - 5*x^3 -                                      15*x)*lam6 - 7/32*x
}    

#' @rdname q1s
#' @export
q1k <- function(x, r, k12, k31) {
  -1/6*He2(x/r)*k31/r^3 - He0(x/r)*k12/r
}	
#' @rdname q1s
#' @export
q2k <- function(x, r, k12, k22, k31, k41) {
  -1/72*He5(x/r)*k31^2/r^6 - 1/24*(4*k12*k31 + k41)*He3(x/r)/r^4 - 
    1/2*(k12^2 + k22)*He1(x/r)/r^2
}    
#' @rdname q1s
#' @export
q3k <- function(x, r, k12, k13, k22, k31, k32, k41, k51) {
  -1/1296*He8(x/r)*k31^3/r^9 - 1/144*(2*k12*k31^2 + k31*k41)*He6(x/r)/r^7 -
    1/120*(10*k12^2*k31 + 10*k22*k31 + 5*k12*k41 + k51)*He4(x/r)/r^5 -
    1/6*(k12^3 + 3*k12*k22 + k32)*He2(x/r)/r^3 - He0(x/r)*k13/r
}  
#' @rdname q1s
#' @export
q4k <- function(x, r, k12, k13, k22, k23, k31, k32, k41, k42, k51, k61) {
  -1/31104*He11(x/r)*k31^4/r^12 - 1/5184*(4*k12*k31^3 +                                             3*k31^2*k41)*He9(x/r)/r^10 - 1/5760*(40*k12^2*k31^2 + 40*k22*k31^2 +                                         40*k12*k31*k41 + 5*k41^2 + 8*k31*k51)*He7(x/r)/r^8 - 
    1/720*(20*k12^3*k31 + 60*k12*k22*k31 + 15*k12^2*k41 + 20*k31*k32 + 
             15*k22*k41 + 6*k12*k51 + k61)*He5(x/r)/r^6 -
    1/24*(k12^4 + 6*k12^2*k22 + 3*k22^2 + 4*k13*k31 + 4*k12*k32 +
            k42)*He3(x/r)/r^4 - 1/2*(2*k12*k13 + k23)*He1(x/r)/r^2
}  

#' Standalone Edgeworth expansion functions for a t-statistic
#' 
#' Calculate values of 1 - 4-term Edgeworth expansions (2 - 5 order) for short
#' (ordinary one-sample) and general versions of t-statistic.
#' 
#' \code{Ft1, ..., Ft4} implement a short version of EE that can be used for 
#' ordinary one-sample t-statistic. Note that in this case variance adjustment
#' \eqn{r^2} is equal to 1 for a t-statistic with naive biased variance estimate
#' and \eqn{(n - 1)/n} for a standard t with unbiased variance estimate.
#' \code{Ft1gen, ..., Ft4gen} implement a general version, which can be used for
#' any t-statistic.
#' 
#' @inheritParams q1s
#' @param k12,k13,k22,k23,k31,k32,k41,k42,k51,k61 cumulant components - values 
#'   calculated from sample statistics or distribution parameters.
#' @param norm if \code{TRUE}, normal distribution is used as a base for
#'   Edgeworth expansions, if \code{FALSE} - Student's t-distribution. Defaults
#'   to \code{FALSE}.
#' @param df degrees of freedom for Student's t-distribution if \code{norm =
#'   FALSE}. For a general version, provide a single value for the first order
#'   approximation (zero term).
#'   
#' @return The value of Edgeworth expansion of a corresponding order (\code{Ft1}
#'   and \code{Ft1gen} for a 1-term or 2nd order EE and so on). (or values, one
#'   for each x?)
#'   
#' @seealso \code{\link{q1s}}, \code{\link{q1k}} for \code{q()} functions used in
#'   EE terms.
#'   
#' @examples 
#' 
#' @export
Ft1 <- function(x, n, lam, r = 1, norm = FALSE, df = n - 1) {
  if (norm) pnorm(x/r)  + n^(-1/2)*q1s(x/r, lam)*dnorm(x/r)
  else      pt(x/r, df) + n^(-1/2)*q1s(x/r, lam)*dt(x/r, df + 2)                  
}
#' @rdname Ft1
#' @export
Ft2 <- function(x, n, lam, r = 1, norm = FALSE, df = n - 1) {
  Ft1(x, n, r, lam, norm, df) + 
    n^(-1)  *q2s(x/r, lam)*(norm*dnorm(x/r) + (1 - norm)*dt(x/r, df + 5))
}
#' @rdname Ft1
#' @export
Ft3 <- function(x, n, lam, r = 1, norm = FALSE, df = n - 1) {
  Ft2(x, n, lam, r, norm, df) + 
    n^(-3/2)*q3s(x/r, lam)*(norm*dnorm(x/r) + (1 - norm)*dt(x/r, df + 8))
}	
#' @rdname Ft1
#' @export
Ft4 <- function(x, n, lam, r = 1, norm = FALSE, df = n - 1) {
  Ft3(x, n, lam, norm, df) + 
    n^(-2)  *q4s(x/r, lam)*(norm*dnorm(x/r) + (1 - norm)*dt(x/r, df + 11))
}	

Ft1gen <- function(x, n, r, k12, k31, norm = FALSE, df = NULL) {
  if (norm) pnorm(x/r)  + n^(-1/2)*q1k(x, r, k12, k31)*dnorm(x/r)
  else      pt(x/r, df) + n^(-1/2)*q1k(x, r, k12, k31)*dt(x/r, df + 2)                  
}
#' @rdname Ft1
#' @export
Ft2gen <- function(x, n, r, k12, k22, k31, k41, norm = FALSE, df = NULL) {
  Ft1gen(x, n, r, k12, k31, norm, df) + n^(-1)*q2k(x, r, k12, k22, k31, k41)*(norm*dnorm(x/r) + (1 - norm)*dt(x/r, df + 5))
}	
#' @rdname Ft1
#' @export
Ft3gen <- function(x, n, r, k12, k13, k22, k31, k32, k41, k51, norm = FALSE, df = NULL) {
  Ft2gen(x, n, r, k12, k22, k31, k41, norm, df) + n^(-3/2)*q3k(x, r, k12, k13, k22, k31, k32, k41, k51)*(norm*dnorm(x/r) + (1 - norm)*dt(x/r, df + 8))
}	
#' @rdname Ft1
#' @export
Ft4gen <- function(x, n, r, k12, k13, k22, k23, k31, k32, k41, k42, k51, k61, norm = FALSE, df = NULL) {
  Ft3gen(x, n, r, k12, k13, k22, k31, k32, k41, k51, norm, df) + n^(-2)*q4k(x, r, k12, k13, k22, k23, k31, k32, k41, k42, k51, k61)*(norm*dnorm(x/r) + (1 - norm)*dt(x/r, df + 11))
}	

