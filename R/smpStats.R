#' Sample statistics for Edgeworth expansions
#' 
#' Calculate sample statistics needed for Edgeworth expansions.
#' 
#' @param smp sample.
#' @param a vector of the same length as \code{smp} specifying categories of 
#'   observations (should contain two unique values). Treatment code is assumed
#'   to have a higher numeric value than control (relevant for \code{type =
#'   "Welch"}).
#' @param unbiased.mom \code{logical} value indicating if unbiased estimators 
#'   for third through sixth central moments should be used.
#' @inheritParams empEdge
#' @param moder \code{logical} value indicating if Edgeworth expansions for a
#'   moderated t-statistic will be used. If \code{TRUE}, prior information
#'   (\code{d0} and \code{s20}) and posterior variance should be provided.
#' @param d0 prior degrees of freedom (needed if \code{moder = TRUE}).
#' @param s20 prior value for variance (needed if \code{moder = TRUE}).
#' @param varpost posterior variance (needed if \code{moder = TRUE}).
#'   
#' @return A named vector of sample statistics to be used in Edgeworth
#'   exansions. The calculated statistics and corresponding names are: \itemize{
#'   \item for ordinary one-sample t-statistic: scaled cumulants named
#'   \code{"lam3", "lam4", "lam5", "lam6"}; \item for moderated one-sample
#'   t-statistic: central moment estimates named \code{"mu2", "mu3", "mu4",
#'   "mu5", "mu6"}, \code{A}, \code{B}, and prior degrees of freedom named
#'   \code{"d0"}; \item for ordinary two-sample t-statistic: central moment
#'   estimates repeated twice since the same distribution is assumed for two
#'   groups, named \code{"mu_x2", "mu_x3", "mu_x4", "mu_x5", "mu_x6"} and
#'   \code{"mu_y2", "mu_y3", "mu_y4", "mu_y5", "mu_y6"}, \code{A, B_x, B_y,
#'   b_x}, and \code{b_y}; \item  for moderated two-sample t-statistic:
#'   estimates of the same quantities as for ordinary t (with different
#'   estimators); additionally, prior degrees of freedom named \code{"d0"} is
#'   included; \item for Welch t-test: estimates of the same quantities as for
#'   ordinary t-statistic (with different estimators). In this case, central
#'   moment estimates for treatment and control groups are different.}
#'   
#' @seealso \code{\link{tailDiag}}, \code{\link{makeFx}}, and
#'   \code{\link{makeQx}} for functions that use \code{stats} argument
#'   corresponding to the output of \code{smpStats()}.
#' @export
smpStats <- function(smp, a = NULL, unbiased.mom = TRUE, type = NULL, 
                     moder = FALSE, d0 = NULL, s20 = NULL, varpost = NULL) {
  n <- length(smp)
  if (is.null(a)) {
    type <- "one-sample"
  } else {
    if (length(a) != n) stop("design does not match data")
    if (length(unique(a)) == 1) {
      type <- "one-sample"
    } else if (length(unique(a)) == 2) {
      nx <- sum(a == max(a))
      ny <- n - nx
      n <- n/2
      if (is.null(type)) {
        type <- "two-sample"
      }
    } else {
      stop("more than two categories in the sample")
    }
  }
  
  if (type == "one-sample") {
    if (unbiased.mom) {
      mu <- getMomEdgeUnb(smp)
    } else {
      mu <- getMomEdgeBias(smp)
    }
    if (!moder) {
      return(getLam(mu))
    } else {
      if (is.null(d0))      stop("Please provide d0")
      if (is.null(s20))     stop("Please provide s20")
      if (is.null(varpost)) stop("Please provide posterior variance")
      mu[1] <- varpost
      A <- (d0*s20 + n*mu[1])/(d0 + n - 1)
      names(A) <- NULL  # remove carryover from mu
      B <- n/(d0 + n - 1)
      other.stats <- c(A, B, d0)
      return(c(mu, A = A, B = B, d0 = d0))
    }
  }
  
  if (type == "two-sample") {
    bx <- n/nx
    by <- n/ny
    Cxy <- (nx + ny)/(nx + ny - 2)
    Bx <- Cxy*by
    By <- Cxy*bx
    if (unbiased.mom) {
      mu <- getMomEdgeUnb2(smp, a)
    } else {
      mu <- getMomEdgeBias2(smp, a)
    }
    if (!moder) {
      A <- Cxy*(bx + by)*mu[1]
      names(A) <- NULL  # remove carryover from mu
      stats <- rep(mu, 2)
      names(stats) <- c(paste("mu_x", 2:6, sep = ""), 
                        paste("mu_y", 2:6, sep = ""))
      return(c(stats, A = A, B_x = Bx, B_y = By, b_x = bx, b_y = by))
    } else {
      if (is.null(d0))      stop("Please provide d0")
      if (is.null(s20))     stop("Please provide s20")
      if (is.null(varpost)) stop("Please provide posterior variance")
      mu[1] <- varpost
      dg <- nx + ny - 2
      A <- (bx + by)*(d0*s20 + Cxy*dg*mu[1])/(d0 + dg)
      names(A) <- NULL  # remove carryover from mu
      Bx <- Bx*dg/(d0 + dg)
      By <- By*dg/(d0 + dg)
      stats <- rep(mu, 2)
      names(stats) <- c(paste("mu_x", 2:6, sep = ""), 
                        paste("mu_y", 2:6, sep = ""))
      return(c(stats, A = A, B_x = Bx, B_y = By, b_x = bx, b_y = by, d0 = d0))
    }
  }
  
  if (type %in% c("welch", "Welch")) {
    if (is.null(a) | length(unique(a)) != 2) stop("a does not match test type")
    Cx <- nx/(nx - 1)
    Cy <- ny/(ny - 1)
    Bx <- Cx*bx
    By <- Cy*by
    treat <- a == max(a)
    smpx <- smp[treat]
    smpy <- smp[!treat]
    if (unbiased.mom) {
      mux <- getMomEdgeUnb(smpx)
      muy <- getMomEdgeUnb(smpy)
    } else {
      mux <- getMomEdgeBias(smpx)
      muy <- getMomEdgeBias(smpy)
    }
    A <- Cx*bx*mux[1] + Cy*by*muy[1]
    names(A) <- NULL  # remove carryover from mu
    stats <- c(mux, muy)
    names(stats) <- c(paste("mu_x", 2:6, sep = ""), 
                      paste("mu_y", 2:6, sep = ""))
    return(c(stats, A = A, B_x = Bx, B_y = By, b_x = bx, b_y = by))
  }
}