#' Empirical Edgeworth expansions for high-dimensional data analysis
#' 
#' Higher order inference for one- and two-sample t-tests in high-dimensional 
#' data. Includes ordinary and moderated t-statistic and Welch t-test.
#' 
#' Unadjusted p-values are calculated for five orders of approximation for 
#' ordinary and moderated (empirical Bayes method) t-statistics; prior 
#' information and moderated t-statistics are calculated with \code{limma} 
#' package. If prior degrees of freedom is \code{Inf}, higher orders are 
#' provided for ordinary t-statistic only. In a two-sample test, when the 
#' variances (and distributions) are not assumed to be equal and Welch t-test is
#' performed, only results for ordinary t-statistic are provided. Variance 
#' adjustment is used for all the orders (see the paper) and therefore even 
#' first-order results might differ slightly from the regular Student's 
#' t-distribution approximation. When a first-order p-value (for moderated 
#' t-statistic if relevant) is greater than provided significance level 
#' \code{alpha}, no higher order inference is calculated.
#' 
#' Tail diagnostic investigating Edgeworth expansion (EE) tail behavior is 
#' performed for each relevant feature (row of data); if EE of a particular 
#' order is not determined to be helpful, p-value of a previous order is 
#' provided in its place.
#' 
#' For better performance of a second order, using \code{unbiased.mom = TRUE} is
#' recommended. For variance estimate, posterior variance is used for moderated
#' t-statistic and unbiased/pooled variance for ordinary t.
#' 
#' @param dat data matrix with rows corresponding to features. The number of 
#'   columns is a sample size and number of rows is a number of tests. If the 
#'   number of tests is \code{1}, \code{dat} can be a vector.
#' @param a treatment vector; the length has to correspond to the number of 
#'   columns in \code{dat}. Treatment code is assumed to have a higher numeric 
#'   value than control.
#' @param side the test can be one-sided or two-sided. For a one-sided test, the
#'   values are \code{"left"} or \code{"right"}.
#' @param type type of the test with possible values \code{"one-sample"}, 
#'   \code{"two-sample"}, and \code{"Welch"}. For regular one- and two-sample 
#'   tests the value is inferred from \code{a} but for Welch t-test it needs to 
#'   be specified.
#' @param unbiased.mom \code{logical} value indicating if unbiased estimators 
#'   for third through sixth central moments should be used. Defaults to 
#'   \code{TRUE}.
#' @param alpha significance level.
#' @param ncheck number of intervals for tail diagnostic.
#' @param lim tail region for tail diagnostic. Provide the endpoints for the 
#'   right tail (positive values).
#'   
#' @return A matrix with the same number of rows as \code{dat}, each row 
#'   providing p-values for five orders of Edgeworth expansions (0 - 4-term 
#'   expansions) for a corresponding feature (row of data). Where applicable, 
#'   p-values will be provided for both ordinary and moderated t-statistics (10 
#'   columns, five orders each); for Welch t-test the matrix will have five 
#'   columns, and if prior degrees of freedom is \code{Inf}, only first order 
#'   p-values are returned for moderated t-statistic (six columns); note that
#'   variance adjustment \eqn{r^2} is 1 in that case.
#'   
#' @seealso \code{\link{tailDiag}} for tail daignostic, \code{\link{makeFx}}, 
#'   \code{\link{Ftshort}}, and \code{\link{Ftgen}} for calculating Edgeworth 
#'   expansions of orders 1 to 5, and \code{\link{smpStats}} for extracting 
#'   statistics needed to calculate EE from a sample.
#'   
#' @examples
#' # simulate a data set
#' nx <- 10           # sample size
#' m  <- 1e4          # number of tests
#' ns <- 0.05*m       # number of significant features
#' dat <- matrix(rgamma(m*nx, shape = 3) - 3, nrow = m)
#' shifts <- runif(ns, 1, 5)
#' dat[1:ns, ] <- dat[1:ns, ] - shifts
#' # run
#' res <- empEdge(dat)
#' head(res, 3)
#' 
#' # one test (data not high-dimensional)
#' empEdge(dat[1, ], side = "left", unbiased.mom = FALSE, alpha = 0.1)
#' 
#' # Welch test
#' ny <- 12
#' dat2 <- cbind(matrix(rnorm(m*ny), nrow = m), dat)
#' treat <- rep(0:1, c(ny, nx))
#' res <- empEdge(dat2, treat, type = "Welch", ncheck = 50, lim = c(1, 10))
#' head(res, 3)
#' 
#' # prior degrees of freedom not finite
#' if (require(limma)) {
#'   d0 <- 0
#'   while (is.finite(d0)) {
#'     dat <- matrix(rnorm(m*nx), nrow = m)
#'     dat[1:ns, ] <- dat[1:ns, ] + shifts
#'     fit <- lmFit(dat, rep(1, nx)) 
#'     d0 <- ebayes(fit)$df.prior
#'   }
#' } 
#' res <- empEdge(dat, side = "right")
#' head(res, 3)  
#' 
#' @export
#' @useDynLib edgee

empEdge <- function(dat, a = NULL, side = "two-sided", type = NULL, 
                    unbiased.mom = TRUE, alpha = 0.05, ncheck = 30, 
                    lim = c(1, 7)) { 
  
  if (is.null(dim(dat)) || any(dim(dat) == 1)) {
    dat <- matrix(dat, nrow = 1)
  } 
  n <- ncol(dat)
  m <- nrow(dat)

  if (is.null(a)) {
  	type <- "one-sample"
  	design.mat <- rep(1, n)
  	n <- c(n, 0)
  } else {
		if (length(a) != n) stop("design does not match data")
		if (length(unique(a)) == 1) {
			type <- "one-sample"
			design.mat <- rep(1, n)
		} else if (length(unique(a)) == 2) {
			# assume treatment code has higher numeric value than control
		  treat <- a == max(a)
		  n <- c(sum(treat), sum(!treat))	
		  if (m == 1) {
		    dat <- matrix(c(dat[, treat], dat[, !treat]), nrow = 1)
		  } else {
		    dat <- cbind(dat[, treat], dat[, !treat])
		  }
		  design.mat <- model.matrix(~ rep(1:0, n))
			if (is.null(type)) {
				type <- "two-sample"
			}
		} else {
			stop("more than two categories in the sample")
		}
  }
	
#	if (is.null(unbiased.mom)) {
#	  unbiased.mom <- ifelse(type %in% c("one-sample", "Welch", "welch"), 
#	                         TRUE, FALSE)
#	}

	if (type %in% c("welch", "Welch")) {
	  if (is.null(a) | length(unique(a)) != 2) stop("a does not match test type")
	  co <- .C("empEdgeWelch", dat = as.double(dat), 
	           nc = as.integer(n), nr = as.integer(m), alpha = as.double(alpha), 
	           side = as.character(side), ncheck = as.integer(ncheck), 
	           lim = as.double(lim), unbmom = as.integer(unbiased.mom), 
	           pval = as.double(rep(0, m*5)))
	  if (side == "two-sided") {
	    co$pval <- 2*co$pval
	  }
	  co$pval  <- matrix(co$pval, nrow = m, ncol = 5)
	  colnames(co$pval)  <- c("t-dist",  paste("term",  1:4, sep = ""))
	  return(co$pval)
	}
  
  if (!require(limma)) stop("Please install package 'limma' from Bioconductor")

	if (type %in% c("one-sample", "two-sample")) {
		fit <- lmFit(dat, design.mat, weights = NULL) 
		nf <- sum(1/n[n != 0])             # one-smp: 1/n, two-smp: 1/nx + 1/ny
		t.ord <- fit$coefficients[, ncol(fit$coefficients)]/(sqrt(nf)*fit$sigma)
		if (m == 1) {
		  co <- .C("empEdgeOrd", dat = as.double(dat),
		           nc = as.integer(n), nr = as.integer(m), alpha = as.double(alpha),
		           onesmp = as.integer(type == "one-sample"), 
		           side = as.character(side), ncheck = as.integer(ncheck), 
		           lim = as.double(lim), unbmom = as.integer(unbiased.mom), 
		           t.ord = as.double(t.ord), pval = as.double(rep(0, m*5)))
		  if (side == "two-sided") {
		    co$pval  <- 2*co$pval
		  }
		  names(co$pval) <- c("t-dist",  paste("term",  1:4, sep = ""))
		  return(co$pval)
		}
		
		fbay <- ebayes(fit)                           
		s20 <- fbay$s2.prior
		d0  <- fbay$df.prior
		if (is.infinite(d0)) {
		  co <- .C("empEdgeOrd", dat = as.double(dat),
		           nc = as.integer(n), nr = as.integer(m), alpha = as.double(alpha),
		           onesmp = as.integer(type == "one-sample"), 
		           side = as.character(side), ncheck = as.integer(ncheck), 
		           lim = as.double(lim), unbmom = as.integer(unbiased.mom), 
		           t.ord = as.double(t.ord), pval = as.double(rep(0, m*5)))
		  tM <- fbay$t[, ncol(fbay$t)]
		  pM <- fbay$p.value[, ncol(fbay$p.value)]
		  if (side == "two-sided") {
		    co$pval <- 2*co$pval
		  } else {
		    pM <- pM/2
		    pM <- mapply(function(t, p, side) {
		      if((t < 0 & side == "left") | (t > 0 & side == "right")) {
		        return(p)
		      } else {
		        return(1 - p)
		      }}, tM, pM, side)
		  }
		  co$pval <- matrix(co$pval, nrow = m, ncol = 5)
		  colnames(co$pval) <- c("t-dist",  paste("term",  1:4, sep = ""))
		  return(cbind(co$pval, "tM-dist" = pM))
		}
		
		co <- .C("empEdge", dat = as.double(dat),
		         nc = as.integer(n), nr = as.integer(m),
		         d0 = as.double(d0), s20 = as.double(s20), alpha = as.double(alpha),
		         onesmp = as.integer(type == "one-sample"), 
		         side = as.character(side), ncheck = as.integer(ncheck), 
		         lim = as.double(lim), unbmom = as.integer(unbiased.mom), 
		         varpost = as.double(fbay$s2.post), t.ord = as.double(t.ord), 
		         t.mod = as.double(fbay$t[, ncol(fbay$t)]),
		         pval = as.double(rep(0, m*5)), pvalM = as.double(rep(0, m*5)))
		
		if (side == "two-sided") {
		  co$pval  <- 2*co$pval
		  co$pvalM <- 2*co$pvalM
		}
		co$pval  <- matrix(co$pval, nrow = m, ncol = 5)
		co$pvalM <- matrix(co$pvalM, nrow = m, ncol = 5)
		colnames(co$pval)  <- c("t-dist",  paste("term",  1:4, sep = ""))
		colnames(co$pvalM) <- c("tM-dist", paste("termM", 1:4, sep = ""))
		return(cbind(co$pval, co$pvalM)) 
	} else {
	  stop("type is not recognized")
	}
}

