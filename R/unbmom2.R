#' @family unbiased estimates
#' @inherit M6two title description params
#' @return Pooled variance estimate.
#' @examples 
#' n1 <- 10
#' n2 <- 8
#' shp <- 3
#' smp1 <- rgamma(n1, shape = shp) - shp
#' smp2 <- rgamma(n2, shape = shp)
#' for (j in 2:6) {
#'   assign(paste("m", j, sep = ""), 
#'          mean(c((smp1 - mean(smp1))^j, (smp2 - mean(smp2))^j)))
#' }          
#' M2two(m2, n1, n2)
#' @export
M2two <- function(m2, n_x, n_y) m2*(n_x + n_y)/(n_x + n_y - 2)
#' @family unbiased estimates
#' @inherit M6two title description params
#' @return Pooled estimate of a third central moment. 
#' @examples 
#' n1 <- 10
#' n2 <- 8
#' shp <- 3
#' smp1 <- rgamma(n1, shape = shp) - shp
#' smp2 <- rgamma(n2, shape = shp)
#' for (j in 2:6) {
#'   assign(paste("m", j, sep = ""), 
#'          mean(c((smp1 - mean(smp1))^j, (smp2 - mean(smp2))^j)))
#' }          
#' M3two(m3, n1, n2)
#' @export
M3two <- function(m3, n_x, n_y) m3*n_x*n_y*(n_x + n_y)/(n_x^2*n_y + n_x*n_y^2 - 6*n_x*n_y + 2*n_x + 2*n_y)
#' @family unbiased estimates
#' @inherit M6two title description params
#' @return Pooled estimate of a fourth central moment. 
#' @examples 
#' n1 <- 10
#' n2 <- 8
#' shp <- 3
#' smp1 <- rgamma(n1, shape = shp) - shp
#' smp2 <- rgamma(n2, shape = shp)
#' for (j in 2:6) {
#'   assign(paste("m", j, sep = ""), 
#'          mean(c((smp1 - mean(smp1))^j, (smp2 - mean(smp2))^j)))
#' }          
#' M4two(m2, m4, n1, n2)
#' @export
M4two <- function(m2, m4, n_x, n_y) n_x*n_y*(3*m2^2*(n_x^2 + 2*n_x*n_y + n_y^2)*(4*n_x^2*n_y^2 - 5*n_x^2*n_y + 3*n_x^2 - 5*n_x*n_y^2 + 3*n_y^2) - m4*n_x*n_y*(n_x + n_y)*(n_x^3*n_y + 2*n_x^2*n_y^2 - 5*n_x^2*n_y + n_x*n_y^3 - 5*n_x*n_y^2 + 12*n_x*n_y - 3*n_x - 3*n_y))/(3*(n_x^2*n_y + n_x*n_y^2 - 4*n_x*n_y + n_x + n_y)*(4*n_x^2*n_y^2 - 5*n_x^2*n_y + 3*n_x^2 - 5*n_x*n_y^2 + 3*n_y^2) - (n_x^3*n_y^2 + n_x^2*n_y^3 - 8*n_x^2*n_y^2 + 6*n_x^2*n_y - 3*n_x^2 + 6*n_x*n_y^2 - 3*n_y^2)*(n_x^3*n_y + 2*n_x^2*n_y^2 - 5*n_x^2*n_y + n_x*n_y^3 - 5*n_x*n_y^2 + 12*n_x*n_y - 3*n_x - 3*n_y))
#' @family unbiased estimates
#' @inherit M6two title description params
#' @return Pooled estimate of squared variance \eqn{\mu_2^2}{\mu[2]^2}, where
#'   \eqn{\mu_2}{\mu[2]} is a variance.
#' @examples 
#' n1 <- 10
#' n2 <- 8
#' shp <- 3
#' smp1 <- rgamma(n1, shape = shp) - shp
#' smp2 <- rgamma(n2, shape = shp)
#' for (j in 2:6) {
#'   assign(paste("m", j, sep = ""), 
#'          mean(c((smp1 - mean(smp1))^j, (smp2 - mean(smp2))^j)))
#' }          
#' M2pow2two(m2, m4, n1, n2)
#' @export
M2pow2two <- function(m2, m4, n_x, n_y) n_x*n_y*(-m2^2*(n_x^2 + 2*n_x*n_y + n_y^2)*(n_x^3*n_y^2 + n_x^2*n_y^3 - 8*n_x^2*n_y^2 + 6*n_x^2*n_y - 3*n_x^2 + 6*n_x*n_y^2 - 3*n_y^2) + m4*n_x*n_y*(n_x + n_y)*(n_x^2*n_y + n_x*n_y^2 - 4*n_x*n_y + n_x + n_y))/(3*(n_x^2*n_y + n_x*n_y^2 - 4*n_x*n_y + n_x + n_y)*(4*n_x^2*n_y^2 - 5*n_x^2*n_y + 3*n_x^2 - 5*n_x*n_y^2 + 3*n_y^2) - (n_x^3*n_y^2 + n_x^2*n_y^3 - 8*n_x^2*n_y^2 + 6*n_x^2*n_y - 3*n_x^2 + 6*n_x*n_y^2 - 3*n_y^2)*(n_x^3*n_y + 2*n_x^2*n_y^2 - 5*n_x^2*n_y + n_x*n_y^3 - 5*n_x*n_y^2 + 12*n_x*n_y - 3*n_x - 3*n_y))
#' @family unbiased estimates
#' @inherit M6two title description params
#' @inheritParams M2M3two
#' @return Pooled estimate of a fifth central moment. 
#' @examples 
#' n1 <- 10
#' n2 <- 8
#' shp <- 3
#' smp1 <- rgamma(n1, shape = shp) - shp
#' smp2 <- rgamma(n2, shape = shp)
#' for (j in 2:6) {
#'   assign(paste("m", j, sep = ""), 
#'          mean(c((smp1 - mean(smp1))^j, (smp2 - mean(smp2))^j)))
#' }          
#' M5two(m2, m3, m5, n1, n2)
#' @export
M5two <- function(m2, m3, m5, n_x, n_y) n_x^2*n_y^2*(10*m2*m3*(n_x^2 + 2*n_x*n_y + n_y^2)*(-2*n_x^3*n_y^3 + 5*n_x^3*n_y^2 - 8*n_x^3*n_y + 4*n_x^3 + 5*n_x^2*n_y^3 - 8*n_x*n_y^3 + 4*n_y^3) + m5*n_x*n_y*(n_x + n_y)*(n_x^4*n_y^2 + 2*n_x^3*n_y^3 - 12*n_x^3*n_y^2 + 2*n_x^3*n_y + n_x^2*n_y^4 - 12*n_x^2*n_y^3 + 60*n_x^2*n_y^2 - 42*n_x^2*n_y + 20*n_x^2 + 2*n_x*n_y^3 - 42*n_x*n_y^2 + 20*n_y^2))/(10*(n_x^3*n_y^2 + n_x^2*n_y^3 - 8*n_x^2*n_y^2 + 5*n_x^2*n_y - 2*n_x^2 + 5*n_x*n_y^2 - 2*n_y^2)*(-2*n_x^3*n_y^3 + 5*n_x^3*n_y^2 - 8*n_x^3*n_y + 4*n_x^3 + 5*n_x^2*n_y^3 - 8*n_x*n_y^3 + 4*n_y^3) + (n_x^4*n_y^3 + n_x^3*n_y^4 - 10*n_x^3*n_y^3 + 10*n_x^3*n_y^2 - 10*n_x^3*n_y + 4*n_x^3 + 10*n_x^2*n_y^3 - 10*n_x*n_y^3 + 4*n_y^3)*(n_x^4*n_y^2 + 2*n_x^3*n_y^3 - 12*n_x^3*n_y^2 + 2*n_x^3*n_y + n_x^2*n_y^4 - 12*n_x^2*n_y^3 + 60*n_x^2*n_y^2 - 42*n_x^2*n_y + 20*n_x^2 + 2*n_x*n_y^3 - 42*n_x*n_y^2 + 20*n_y^2))
#' @family unbiased estimates
#' @inherit M6two title description params
#' @param m5 naive biased fifth central moment estimate \eqn{m_5 = 1/(n_x + n_y)
#'   \sum_{i = 1}^{n_x} ((X_i - \bar{X})^5 + \sum_{i = 1}^{n_y} ((Y_i - 
#'   \bar{Y})^5}{m[5] = mean(c((X - X-bar)^5, (Y - Y-bar)^5))} for vectors 
#'   \code{X} and \code{Y}.
#' @return Pooled estimate of a product of second and third central moments 
#'   \eqn{\mu_2 \mu_3}{\mu[2] \mu[3]}, where \eqn{\mu_2}{\mu[2]} and 
#'   \eqn{\mu_3}{\mu[3]} are second and third central moments respectively.
#' @examples 
#' n1 <- 10
#' n2 <- 8
#' shp <- 3
#' smp1 <- rgamma(n1, shape = shp) - shp
#' smp2 <- rgamma(n2, shape = shp)
#' for (j in 2:6) {
#'   assign(paste("m", j, sep = ""), 
#'          mean(c((smp1 - mean(smp1))^j, (smp2 - mean(smp2))^j)))
#' }          
#' M2M3two(m2, m3, m5, n1, n2)
#' @export
M2M3two <- function(m2, m3, m5, n_x, n_y) n_x^2*n_y^2*(m2*m3*(n_x^2 + 2*n_x*n_y + n_y^2)*(n_x^4*n_y^3 + n_x^3*n_y^4 - 10*n_x^3*n_y^3 + 10*n_x^3*n_y^2 - 10*n_x^3*n_y + 4*n_x^3 + 10*n_x^2*n_y^3 - 10*n_x*n_y^3 + 4*n_y^3) - m5*n_x*n_y*(n_x + n_y)*(n_x^3*n_y^2 + n_x^2*n_y^3 - 8*n_x^2*n_y^2 + 5*n_x^2*n_y - 2*n_x^2 + 5*n_x*n_y^2 - 2*n_y^2))/(10*(n_x^3*n_y^2 + n_x^2*n_y^3 - 8*n_x^2*n_y^2 + 5*n_x^2*n_y - 2*n_x^2 + 5*n_x*n_y^2 - 2*n_y^2)*(-2*n_x^3*n_y^3 + 5*n_x^3*n_y^2 - 8*n_x^3*n_y + 4*n_x^3 + 5*n_x^2*n_y^3 - 8*n_x*n_y^3 + 4*n_y^3) + (n_x^4*n_y^3 + n_x^3*n_y^4 - 10*n_x^3*n_y^3 + 10*n_x^3*n_y^2 - 10*n_x^3*n_y + 4*n_x^3 + 10*n_x^2*n_y^3 - 10*n_x*n_y^3 + 4*n_y^3)*(n_x^4*n_y^2 + 2*n_x^3*n_y^3 - 12*n_x^3*n_y^2 + 2*n_x^3*n_y + n_x^2*n_y^4 - 12*n_x^2*n_y^3 + 60*n_x^2*n_y^2 - 42*n_x^2*n_y + 20*n_x^2 + 2*n_x*n_y^3 - 42*n_x*n_y^2 + 20*n_y^2))
  

