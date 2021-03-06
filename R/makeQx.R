#' Make terms for Edgeworth expansions
#' 
#' Create functions \code{q(x)} used in terms of Edgeworth expansions (EE) for 
#' various versions of t-statistic.
#' 
#' Given sample statistics or distribution parameters, produce the defining 
#' portion of EE terms as a function of \code{x} (quantile). Includes short (for
#' ordinary one-sample t-statistic) and general versions, which can be used for 
#' ordinary and moderated one- and two-sample t-statistics and Welch t-test. 
#' When used in Edgeworth expansions, the argument to a \code{q()} function 
#' should incorporate variance adjustment (see examples); this adjustment is
#' equal to 1 if naive biased variance estimates are used in a one-sample
#' t-statistic or a Welch test.
#' 
#' @inheritParams tailDiag
#' @param type type of Edgeworth expansions to be used: \code{"short"} for 
#'   ordinary one-sample t-statistic, \code{"one-sample"} for a more complicated
#'   version such as a moderated t-statistic, and \code{"two-sample"} for any 
#'   kind of two-sample t-statistic.
#' @param r square root of variance adjustment. Provide if wish to use a value 
#'   that is different from the one that will be calculated from \code{stats} 
#'   for a general version.
#'   
#' @return A list of functions named \code{"q1"}, \code{"q2"}, \code{"q3"}, and 
#'   \code{"q4"} to be used in respective terms of an Edgeworth expansion 
#'   (orders 2 - 5). These functions take a vector of quantiles \code{x} as an 
#'   input and output a vector of the same length as \code{x}. For a general 
#'   version (\code{type = "one-sample"} and \code{type = "two-sample"}), the 
#'   list also includes a calculated value for variance adjustment \code{r}.
#'   
#' @seealso \code{\link{smpStats}} that creates \code{stats} vector from the 
#'   sample and \code{\link{makeFx}} that creates a function returning the 
#'   values of Edgeworth series of orders 1 - 5 as a function of \code{x}.
#'   
#' @examples 
#' n <- 10                 # sample size
#' # Gamma distribution with shape parameter \code{shp}
#' shp <- 3
#' ord <- 3:6              # orders of scaled cumulants
#' lambdas <- factorial(ord - 1)/shp^((ord - 2)/2)
#' qt <- makeQx(lambdas)
#' quant <- -(3:5)
#' qt[[1]](quant) 
#' 
#' # from sample
#' smp <- rgamma(n, shape = shp) - shp
#' stats <- smpStats(smp)
#' t <- sqrt(n)*mean(smp)/sd(smp)
#' r <- sqrt((n - 1)/n)
#' qt <- makeQx(stats)
#' sapply(qt, function(f) f(t/r))
#' # 2nd order (1-term) Edgeworth expansion, standard normal base
#' pnorm(t/r) + n^(-1/2)*qt[[1]](t/r)*dnorm(t/r)
#' 
#' @export

makeQx <- function(stats, type = "short", r = NULL, verbose = TRUE) {  	
  # if stats vector is unnamed, assume it contains lambda's (standardized
  # cumulants)
  if (is.null(names(stats))) {
    if (type != "short") {
      stop("no names provided for stats")
    } else {
      if (verbose) {
        message("no names provided for stats; assumed to be scaled cumulants")
      }	
      lam3 <- stats[1]
      lam4 <- stats[2]
      lam5 <- stats[3]
      lam6 <- stats[4]
    }
  } else {
    # remove everything after dot in names(stats)
    names(stats) <- gsub("\\..*", "", names(stats))  
    if ("d0" %in% names(stats) && is.infinite(stats['d0'])) {
      stop("prior degrees of freedom not finite")
    }
  }
  
  # one-sample, biased or unbiased (for unbiased, use x1 = sqrt(C)*x)
  if (type == "short") {
    lam3 <- stats['lam3']  # not required to be in order if named
    lam4 <- stats['lam4']
    lam5 <- stats['lam5']
    lam6 <- stats['lam6']
    names(lam3) <- names(lam4) <- names(lam5) <- names(lam6) <- NULL
    q1 <- function(x) 1/6*(2*x^2 + 1)*lam3
    q2 <- function(x) -1/18*(x^5 + 2*x^3 - 3*x)*lam3^2 - 1/4*x^3 + 
      1/12*(x^3 - 3*x)*lam4 - 3/4*x
    q3 <- function(x) 1/1296*(8*x^8 + 28*x^6 - 210*x^4 - 525*x^2 - 105)*lam3^3 -
      1/144*(4*x^6 - 30*x^4 - 90*x^2 - 15)*lam3*lam4 + 1/24*(2*x^6 -
      3*x^4 - 6*x^2)*lam3 - 1/40*(2*x^4 + 8*x^2 + 1)*lam5
    q4 <- function(x) -1/32*x^7 - 1/1944*(x^11 + 5*x^9 - 90*x^7 - 450*x^5 + 45*x^3 +
      945*x)*lam3^4 - 5/96*x^5 - 1/72*(x^9 - 6*x^7 - 12*x^5 - 18*x^3 -
      9*x)*lam3^2 - 1/288*(x^7 - 21*x^5 + 33*x^3 + 111*x)*lam4^2 +
      1/60*(x^7 + 8*x^5 - 5*x^3 - 30*x)*lam3*lam5 - 7/96*x^3 +
      1/432*(9*x^7 - 63*x^5 + 2*(x^9 - 12*x^7 - 90*x^5 + 36*x^3 +
      261*x)*lam3^2 + 81*x^3 + 189*x)*lam4 - 1/90*(2*x^5 - 5*x^3 -
      15*x)*lam6 - 7/32*x
    return(list(q1 = q1, q2 = q2, q3 = q3, q4 = q4))                   
  }
  
  if        (type == "one-sample") {
    k <- calculateK1smp(stats)
  } else if (type == "two-sample") {
    k <- calculateK2smp(stats)
  }
  names(k) <- NULL  # to remove names from output
  k12 <- k[1]
  k13 <- k[2]
  if (is.null(r)) r <- sqrt(k[3])  # k21 = k[3]
  k22 <- k[4]
  k23 <- k[5]
  k31 <- k[6]
  k32 <- k[7]
  k41 <- k[8]
  k42 <- k[9]
  k51 <- k[10]
  k61 <- k[11]

  q1 <- function(x) -1/6*He2(x)*k31/r^3 - He0(x)*k12/r
  q2 <- function(x) -1/72*He5(x)*k31^2/r^6 - 1/24*(4*k12*k31 + k41)*He3(x)/r^4 - 
    1/2*(k12^2 + k22)*He1(x)/r^2
  q3 <- function(x) -1/1296*He8(x)*k31^3/r^9 - 1/144*(2*k12*k31^2 + k31*k41)*He6(x)/r^7 -
    1/120*(10*k12^2*k31 + 10*k22*k31 + 5*k12*k41 + k51)*He4(x)/r^5 -
    1/6*(k12^3 + 3*k12*k22 + k32)*He2(x)/r^3 - He0(x)*k13/r
  q4 <- function(x) -1/31104*He11(x)*k31^4/r^12 - 1/5184*(4*k12*k31^3 + 
    3*k31^2*k41)*He9(x)/r^10 - 1/5760*(40*k12^2*k31^2 + 40*k22*k31^2 +                                40*k12*k31*k41 + 5*k41^2 + 8*k31*k51)*He7(x)/r^8 - 
    1/720*(20*k12^3*k31 + 60*k12*k22*k31 + 15*k12^2*k41 + 20*k31*k32 + 
             15*k22*k41 + 6*k12*k51 + k61)*He5(x)/r^6 -
    1/24*(k12^4 + 6*k12^2*k22 + 3*k22^2 + 4*k13*k31 + 4*k12*k32 +
            k42)*He3(x)/r^4 - 1/2*(2*k12*k13 + k23)*He1(x)/r^2
  
  return(list(q1 = q1, q2 = q2, q3 = q3, q4 = q4, r = r))
}
