% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unbmom1.R
\name{M2M3one}
\alias{M2M3one}
\title{Unbiased central moment estimates}
\usage{
M2M3one(m2, m3, m5, n)
}
\arguments{
\item{m2}{naive biased variance estimate \eqn{m_2 = 1/n \sum_{i = 1}^n ((X_i 
- \bar{X})^2}{m[2] = mean((X - X-bar)^2)} for a vector \code{X}.}

\item{m3}{naive biased third central moment estimate \eqn{m_3 = 1/n \sum_{i =
1}^n ((X_i - \bar{X})^3}{m[3] = mean((X - X-bar)^3)} for a vector \code{X}.}

\item{m5}{naive biased fifth central moment estimate \eqn{m_5 = \sum_{i =
1}^n ((X_i - \bar{X})^5}{m[5] = mean((X - X-bar)^5)} for a vector \code{X}.}

\item{n}{sample size.}
}
\value{
Unbiased estimate of a product of second and third central moments
  \eqn{\mu_2 \mu_3}{\mu[2] \mu[3]}, where \eqn{\mu_2}{\mu[2]} and
  \eqn{\mu_3}{\mu[3]} are second and third central moments respectively.
}
\description{
Calculate unbiased estimates of central moments and their powers and 
products.
}
\examples{
n <- 10
smp <- rgamma(n, shape = 3)
for (j in 2:6) {
  assign(paste("m", j, sep = ""), mean((smp - mean(smp))^j))
}
M2M3one(m2, m3, m5, n) 
}
\seealso{
Other unbiased estimates: 
\code{\link{M2M3two}()},
\code{\link{M2M4one}()},
\code{\link{M2M4two}()},
\code{\link{M2one}()},
\code{\link{M2pow2one}()},
\code{\link{M2pow2two}()},
\code{\link{M2pow3one}()},
\code{\link{M2pow3two}()},
\code{\link{M2two}()},
\code{\link{M3one}()},
\code{\link{M3pow2one}()},
\code{\link{M3pow2two}()},
\code{\link{M3two}()},
\code{\link{M4one}()},
\code{\link{M4two}()},
\code{\link{M5one}()},
\code{\link{M5two}()},
\code{\link{M6one}()},
\code{\link{M6two}()}
}
\concept{unbiased estimates}
