#' Unbiased central moment estimates
#' 
#' Calculate unbiased estimates of central moments and their powers and products
#' up to specified order.
#' 
#' Unbiased estimates up to the 6th order can be calculated. Second and third 
#' orders contain estimates of the variance and third central moment, fourth 
#' order includes estimates of fourth moment and squared variance 
#' (\eqn{\mu_2^2}{\mu[2]^2}), fifth order - of fifth moment and a product of 
#' second and third moments (\eqn{\mu_2 \mu_3}{\mu[2] \mu[3]}), sixth order - of
#' sixth moment, a product of second and fourth moments (\eqn{\mu_2
#' \mu_4}{\mu[2] \mu[4]}), squared third moment (\eqn{\mu_3^2}{\mu[3]^2}), and
#' cubed variance (\eqn{\mu_2^3}{\mu[2]^3}).
#' 
#' @param smp sample.
#' @param order highest order of the estimates to calclulate. Estimates of lower
#'   orders will be included.
#'   
#' @return A named vector of estimates of central moments and their powers and
#'   products up to \code{order}. The highest order available is 6th. The names
#'   of the elements are \code{"M2", "M3", "M4", "M5", "M6"} for corresponding
#'   central moments, \code{"M2M3", "M2M4"} for products of the moments (second
#'   and third, second and fourth), and \code{"M2pow2", "M2pow3", "M3pow2"} for
#'   powers of the moments - corresponding to estimates of squared variance,
#'   cubed variance, and squared third moment.
#' @seealso \code{\link{unbMom2smp}} for two-sample pooled estimates.
#' @examples 
#' smp <- rgamma(10, shape = 3)
#' unbMom(smp, 6)
#' @export
unbMom <- function(smp, order) {
  n <- length(smp)
  order <- trunc(order)
  if (order < 0) {
    return()
  }
  if (order == 0) {
    return(1)
  }
  if (order == 1) {
    return(0)
  }
  res <- c(M2 = var(smp))
  if (order == 2) {
    return(res)
  }	
  m1 <- mean(smp)
  m3 <- mean((smp - m1)^3)
  M3 <- M3one(m3, n)
  res <- c(res, M3 = M3)
  if (order == 3) {
    return(res)
  }	
  m2 <- mean((smp - m1)^2)
  m4 <- mean((smp - m1)^4)
  M2pow2 <- M2pow2one(m2, m4, n)
  M4     <- M4one(    m2, m4, n)
  res <- c(res, M2pow2 = M2pow2, M4 = M4)
  if (order == 4) {
    return(res)
  }
  m5 <- mean((smp - m1)^5)
  M2M3 <- M2M3one(m2, m3, m5, n)
  M5   <- M5one(  m2, m3, m5, n)
  res <- c(res, M2M3 = M2M3, M5 = M5)
  if (order == 5) {
    return(res)
  }
  m6 <- mean((smp - m1)^6)
  M2pow3 <- M2pow3one(m2, m3, m4, m6, n)
  M3pow2 <- M3pow2one(m2, m3, m4, m6, n)
  M2M4   <- M2M4one(  m2, m3, m4, m6, n) 
  M6     <- M6one(    m2, m3, m4, m6, n)
  res <- c(res, M2pow3 = M2pow3, M3pow2 = M3pow2, M2M4 = M2M4, M6 = M6)
  if (order == 6) {
    return(res)
  }
  if (order > 6) {
    warning("orders higher than 6th are not available")
    return(res)
  }
}

#' Pooled central moment estimates - two-sample
#' 
#' Calculate unbiased pooled estimates of central moments and their powers and 
#' products up to specified order.
#' 
#' Pooled estimates up to the 6th order can be calculated. Second and third 
#' orders contain estimates of the variance and third central moment, fourth 
#' order includes estimates of fourth moment and squared variance 
#' (\eqn{\mu_2^2}{\mu[2]^2}), fifth order - of fifth moment and a product of 
#' second and third moments (\eqn{\mu_2 \mu_3}{\mu[2] \mu[3]}), sixth order - of
#' sixth moment, a product of second and fourth moments (\eqn{\mu_2 
#' \mu_4}{\mu[2] \mu[4]}), squared third moment (\eqn{\mu_3^2}{\mu[3]^2}), and 
#' cubed variance (\eqn{\mu_2^3}{\mu[2]^3}).
#' 
#' @inherit unbMom params return
#' @param a vector of the same length as \code{smp} specifying categories of 
#'   observations (should contain two unique values).
#' @seealso \code{\link{unbMom}} for one-sample unbiased estimates.
#' @examples 
#' nsmp <- 23
#' smp2 <- rgamma(nsmp, shape = 3)
#' treatment <- sample(0:1, size = nsmp, replace = TRUE)  
#' unbMom2smp(smp2, treatment, 6)
#' @export
unbMom2smp <- function(smp, a, order) {
  if (length(unique(a)) != 2 | length(a) != length(smp)) {
    stop("design does not match sample")
  }
  smpx <- smp[a == max(a)]
  smpy <- smp[a == min(a)]
  nx <- length(smpx)
  ny <- length(smpy)
  order <- trunc(order)
  if (order < 0) {
    return()
  }
  if (order == 0) {
    return(1)
  }
  if (order == 1) {
    return(0)
  }
  mx1 <- mean(smpx)
  my1 <- mean(smpy)
  m2 <- mean(c((smpx - mx1)^2, (smpy - my1)^2))
  M2 <- M2two(m2, nx, ny)
  res <- c(M2 = M2)
  if (order == 2) {
    return(res)
  }	
  m3 <- mean(c((smpx - mx1)^3, (smpy - my1)^3))
  M3 <- M3two(m3, nx, ny)
  res <- c(res, M3 = M3)
  if (order == 3) {
    return(res)
  }	
  m4 <-  mean(c((smpx - mx1)^4, (smpy - my1)^4))
  M2pow2 <- M2pow2two(m2, m4, nx, ny)
  M4     <- M4two(    m2, m4, nx, ny)
  res <- c(res, M2pow2 = M2pow2, M4 = M4)
  if (order == 4) {
    return(res)
  }
  m5 <- mean(c((smpx - mx1)^5, (smpy - my1)^5))
  M2M3 <- M2M3two(m2, m3, m5, nx, ny)
  M5   <- M5two(  m2, m3, m5, nx, ny)
  res <- c(res, M2M3 = M2M3, M5 = M5)
  if (order == 5) {
    return(res)
  }
  m6 <- mean(c((smpx - mx1)^6, (smpy - my1)^6))
  M2pow3 <- M2pow3two(m2, m3, m4, m6, nx, ny)
  M3pow2 <- M3pow2two(m2, m3, m4, m6, nx, ny)
  M2M4   <- M2M4two(  m2, m3, m4, m6, nx, ny) 
  M6     <- M6two(    m2, m3, m4, m6, nx, ny)
  res <- c(res, M2pow3 = M2pow3, M3pow2 = M3pow2, M2M4 = M2M4, M6 = M6)
  if (order == 6) {
    return(res)
  }
  if (order > 6) {
    warning("orders higher than 6th are not available")
    return(res)
  }
}

#' Central moment estimates for Edgeworth expansions
#' 
#' Calculate central moment estimates for use in Edgeworth expansions for one- 
#' and two-sample t-tests.
#' 
#' For one-sample estimates: \code{getMomEdgeBias()} calculates regular unbiased
#' sample variance and naive biased estimates for 3rd to 6th central moments; 
#' \code{getMomEdgeUnb()} calculates unbiased estimates of 2nd to 6th central 
#' moments. For two-sample estimates, where the two populations are assumed to 
#' have the same variance and higher central moments: \code{getMomEdgeBias2()} 
#' calculates unbiased pooled variance and naive biased 3rd to 6th central 
#' moments; \code{getMomEdgeUnb2()} provides unbiased pooled estimates of 2nd to
#' 6th moments.
#' 
#' @name momEdge
#' @inheritParams unbMom2smp
#' @return A named vector of length \code{5}. The names of the elements 
#'   correspond to the estimands and are \code{"mu2", "mu3", "mu4", "mu5", 
#'   "mu6"}.
#' @seealso \code{\link{getLam}} for calculating scaled cumulants from moments.
#'   
#' @examples 
#' n     <- 10 
#' n2smp <- 23
#' smp <- rgamma(n, shape = 3)        
#' getMomEdgeBias(smp)  # var unbiased, moments 3 - 6 naive biased
#' getMomEdgeUnb( smp)
#' 
#' smp2 <- rgamma(n2smp, shape = 3)
#' treatment <- sample(0:1, size = n2smp, replace = TRUE)  
#' getMomEdgeBias2(smp2, treatment)  # pooled var, moments 3 - 6 naive biased
#' getMomEdgeUnb2( smp2, treatment)
#'  
#' @export
getMomEdgeBias <- function(smp) {
  n <- length(smp)
  mu <- numeric(5)
  mu[1] <- var(smp)
  m1 <- mean(smp)
  for (j in 3:6) {
    mu[j - 1] <-  mean((smp - m1)^j)
  }
  names(mu) <- paste("mu", 2:6, sep = "")
  return(mu)
}
#' @rdname momEdge
#' @export
getMomEdgeUnb <- function(smp) {
  n <- length(smp)
  m1 <- mean(smp)
  for (j in 2:6) {
    assign(paste("m", j, sep = ""), mean((smp - m1)^j))
  }
  M2 <- M2one(m2, n)
  M3 <- M3one(m3, n)
  M4 <- M4one(m2, m4, n)
  M5 <- M5one(m2, m3, m5, n)
  M6 <- M6one(m2, m3, m4, m6, n)
  return(c(mu2 = M2, mu3 = M3, mu4 = M4, mu5 = M5, mu6 = M6))
}
#' @rdname momEdge
#' @export
getMomEdgeBias2 <- function(smp, a) {
  smpx <- smp[a == max(a)]
  smpy <- smp[a == min(a)]
  nx <- length(smpx)
  ny <- length(smpy)
  mx1 <- mean(smpx)
  my1 <- mean(smpy)
  mu <- numeric(5)  
  for (j in 2:6) {
    mu[j - 1] <- mean(c((smpx - mx1)^j, (smpy - my1)^j))
  }
  mu[1] <- M2two(mu[1], nx, ny)
  names(mu) <- paste("mu", 2:6, sep = "")
  return(mu)
}
#' @rdname momEdge
#' @export
getMomEdgeUnb2 <- function(smp, a) {
  smpx <- smp[a == max(a)]
  smpy <- smp[a == min(a)]
  nx <- length(smpx)
  ny <- length(smpy)
  mx1 <- mean(smpx)
  my1 <- mean(smpy)
  for (j in 2:6) {
    assign(paste("m", j, sep = ""), mean(c((smpx - mx1)^j, (smpy - my1)^j)))
  }
  M2 <- M2two(m2, nx, ny)
  M3 <- M3two(m3, nx, ny)
  M4 <- M4two(m2, m4, nx, ny)
  M5 <- M5two(m2, m3, m5, nx, ny)
  M6 <- M6two(m2, m3, m4, m6, nx, ny)
  return(c(mu2 = M2, mu3 = M3, mu4 = M4, mu5 = M5, mu6 = M6))
}

#' Scaled cumulants
#' 
#' Calculate skewness, kurtosis, and 5th and 6th scaled cumulants from central
#' moments or their estimates.
#' 
#' @param mu vector of 2nd - 6th central moments.
#' @return A named vector of 3rd to 6th scaled cumulants. The names of the
#'   elements are \code{"lam3", "lam4", "lam5", "lam6"}.
#' @examples 
#' n <- 10
#' smp <- rgamma(n, shape = 3)
#' getLam(getMomEdgeUnb(smp))   
#' @export
getLam <- function(mu) {
  mu <- c(0, mu)  # to make index = order
  lam <- numeric(4)
  lam[1] <- mu[3]/mu[2]^(3/2)
  lam[2] <- mu[4]/mu[2]^2 - 3
  lam[3] <- mu[5]/mu[2]^(5/2) - 10*lam[1]
  lam[4] <- mu[6]/mu[2]^3 - 15*mu[4]/mu[2]^2 - 10*mu[3]^2/mu[2]^3 + 30
  names(lam) <- paste("lam", 3:6, sep = "")
  return(lam)
}

