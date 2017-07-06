#' Tail diagnostic for Edgeworth expansions
#' 
#' Evaluate tail behavior and usability for 2nd - 5th order (terms 1-4)
#' Edgeworth expansions for a t-statistic.
#' 
#' Determine if the tails of the Edgeworth expansion behave similarly to a
#' cumulative distribution function and therefore can be used for approximating
#' sampling distribution. The results also provide information on the thickness
#' of the tails of that distribution (thicker tail will behave nicely). The
#' function is evaluated according to three criteria: monotonicity, boundedness
#' by 0 and 1, and being no less conservative than the first order approximation
#' (Student's t-distribution).
#' 
#' Both left and right tails are evaluated. t-distribution-based expansions are
#' used.
#' 
#' @param stats named vector of distribution parameters or estimates needed for
#'   Edgeworth expansion. If no names for the vector are provided and the type
#'   is \code{"short"}, elements are assumed to be scaled cumulants
#'   \eqn{\lambda} (3 - 6). If the elements are named, the order of the elements
#'   is arbitrary. Required names: \itemize{ \item for ordinary one-sample
#'   t-statistic only (\code{type = "short"}): \code{"lam3", "lam4", "lam5",
#'   "lam6"} (scaled cumulants); \item for more complicated one-sample versions,
#'   such as a moderated t-statistic (\code{type = "one-sample"}): \code{"mu2",
#'   "mu3", "mu4", "mu5", "mu6"} (central moments), \code{"A"}, and \code{"B"}; 
#'   \item for a two-sample t-statistic (\code{type = "two-sample"}):
#'   \code{"mu_x2", "mu_x3", "mu_x4", "mu_x5", "mu_x6"} (central moments for a
#'   treatment group), \code{"mu_y2", "mu_y3", "mu_y4", "mu_y5", "mu_y6"}
#'   (central moments for a control group), \code{"A", "B_x", "B_y", "b_x"}, and
#'   \code{"b_y"}. Note that
#'   if the same distribution is assumed for the two groups, the values for
#'   central moments for these groups should be the same (e.g. pooled variance
#'   for a second moment). \item  optionally for a moderated t-statistic: \code{"d0"}
#'   (prior degrees of freedom). Ignored if \code{df} is provided. } 
#'   
#' @param n a single value for a sample size summary to be used in Edgeworth
#'   expansion. \strong{Important:} an average (not sum!) of two group sizes for a
#'   two-sample test.
#' @param type which Edgeworth expansions are to be evaluated: \code{"short"}
#'   for ordinary one-sample t-statistic, \code{"one-sample"} for a more
#'   complicated version such as a moderated t-statistic, and
#'   \code{"two-sample"} for any kind of two-sample t-statistic.
#' @param df degrees of freedom for a first order approximation, a parameter of
#'   Student's t-distribution. If not provided, the value will be calculated
#'   based on arguments \code{type} and \code{moder}.
#' @param moder \code{logical} value for calculating augmented degrees of
#'   freedom for a moderated t-statistic; if \code{TRUE}, the value for prior
#'   degrees of freedom shoud be included in \code{stats}. Ignored if \code{df}
#'   is provided.
#' @param ncheck number of intervals for tail diagnostic.
#' @param lim Tail region for tail diagnostic. Provide the endpoints for the
#'   right tail (positive values).
#' @param verbose if \code{TRUE}, the warning message will be printed when the
#'   elements of \code{stats} are not named.
#'   
#' @return A matrix of logical values: four columnts for orders 2 through 5, two
#'   rows for both tails.
#'   
#' @seealso \code{\link{smpStats}} that creates \code{stats} vector from the sample.
#'   
#' @examples
#' # Gamma distribution with shape parameter \code{shp}, ordinary t-statistic:
#' shp <- 3
#' n <- 10
#' lambda <- factorial(2:5)/shp^((1:4)/2)
#' tailDiag(lambda, n, "short")
#' 
#' # Sample statistics:
#' smp <- two-smp?
#' df.prior <- 5  # e.g. obtained from limma::eBayes()
#' smp.stats <- calculateStats(smp, a)
#' tailDiag(c(smp.stats, d0 = df.prior), mean(a), type = "two-sample", moder = TRUE)
#' 
#' @export

# r for q() not needed (reparameterization of quantile x: x = x_1/r)
tailDiag <- function(stats, n, type = "short", df = NULL, moder = FALSE,
                     ncheck = 50, lim = c(1, 10), verbose = TRUE) {
	if (is.null(names(stats))) {
		if (type != "short") {
			stop("no names provided for stats")
		} else {
			if (verbose) {
				warning("no names provided for stats; assumed to be scaled cumulants")
			}
		}
	} else {
		# extract all the needed variables from the named vector
	  vars <- gsub("\\..*", "", names(stats))  # remove everything after dot
		for (i in 1:length(stats)) {
			assign(vars[i], stats[i])
		}
	}
	if (is.null(df)) {
		if (type %in% c("one-sample", "short")) {
			df = n - 1
		} else {
			df = 2*n - 2
		}
		if (moder) {
			df <- df + d0  # d0 should be unpacked from stats
		}
	}

	if (type == "short") {
		if (is.null(names(stats))) {
			qargs <- stats
		} else {
			qargs <- c(lam3, lam4, lam5, lam6)
		}
	} else if (type == "one-sample") {
		qargs <- c(K12one(A, B, mu2, mu3, mu4, mu5, mu6),
		       K13one(A, B, mu2, mu3, mu4, mu5, mu6),
		       K21one(A, B, mu2, mu3, mu4, mu5, mu6),
		       K22one(A, B, mu2, mu3, mu4, mu5, mu6),
		       K23one(A, B, mu2, mu3, mu4, mu5, mu6),
		       K31one(A, B, mu2, mu3, mu4, mu5, mu6),
		       K32one(A, B, mu2, mu3, mu4, mu5, mu6),
		       K41one(A, B, mu2, mu3, mu4, mu5, mu6),
		       K42one(A, B, mu2, mu3, mu4, mu5, mu6),
		       K51one(A, B, mu2, mu3, mu4, mu5, mu6),
		       K61one(A, B, mu2, mu3, mu4, mu5, mu6))
	} else if (type == "two-sample") {
		k <- c(K12two(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                      mu_y2, mu_y3, mu_y4, mu_y5, mu_y6 ),
		       K13two(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                      mu_y2, mu_y3, mu_y4, mu_y5, mu_y6 ),
                K21two(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                       mu_y2, mu_y3, mu_y4, mu_y5, mu_y6 ),
               K22two(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                      mu_y2, mu_y3, mu_y4, mu_y5, mu_y6 ),
               K23two(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                      mu_y2, mu_y3, mu_y4, mu_y5, mu_y6 ),
               K31two(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                      mu_y2, mu_y3, mu_y4, mu_y5, mu_y6 ),
		       K32two(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                      mu_y2, mu_y3, mu_y4, mu_y5, mu_y6 ),
		       K41two(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                      mu_y2, mu_y3, mu_y4, mu_y5, mu_y6 ),
		       K42two(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                      mu_y2, mu_y3, mu_y4, mu_y5, mu_y6 ),
		       K51two(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                      mu_y2, mu_y3, mu_y4, mu_y5, mu_y6 ),
		       K61two(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                      mu_y2, mu_y3, mu_y4, mu_y5, mu_y6 ))
	}

	xright <- seq(lim[1], lim[2], length.out = ncheck + 1)
	xleft  <- -xright
	df <- cumsum(c(df, 2, 3, 3, 3))

	co <- .C("tailDiagR", type = as.character(type),
	                      qargs = as.double(qargs),
	         xleft = as.double(xleft), xright = as.double(xright),
	         n = as.double(n), nch = as.integer(ncheck + 1), df = as.double(df),
	         lthick = as.integer(rep(0, 4)), rthick = as.integer(rep(0, 4)))

	thick <- rbind(as.logical(co$lthick), as.logical(co$rthick))
	colnames(thick) <- paste("term", 1:4)
	rownames(thick) <- paste(c("left  tail", "right tail"), "nice")
	return(thick)
}
