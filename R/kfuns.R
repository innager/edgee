#' Probabilists' Hermite polynomials
#' 
#' @name Hermite
#' @param x numeric vector.
#' 
#' @return A vector of the same length as \code{x}.
#'
He0  <- function(x) 1
#' @rdname Hermite
He1  <- function(x) x
#' @rdname Hermite
He2  <- function(x) x^2 - 1
#' @rdname Hermite
He3  <- function(x) x^3 - 3*x
#' @rdname Hermite
He4  <- function(x) x^4 - 6*x^2 + 3
#' @rdname Hermite
He5  <- function(x) x^5 - 10*x^3 + 15*x
#' @rdname Hermite
He6  <- function(x) x^6 - 15*x^4 + 45*x^2 - 15
#' @rdname Hermite
He7  <- function(x) x^7 - 21*x^5 + 105*x^3 - 105*x
#' @rdname Hermite
He8  <- function(x) x^8 - 28*x^6 + 210*x^4 - 420*x^2 + 105
#' @rdname Hermite
He9  <- function(x) x^9 - 36*x^7 + 378*x^5 - 1260*x^3 + 945*x
#' @rdname Hermite
He10 <- function(x) x^10 - 45*x^8 + 630*x^6 - 3150*x^4 + 4725*x^2 - 945
#' @rdname Hermite
He11 <- function(x) x^11 - 55*x^9 + 990*x^7 - 6930*x^5 + 17325*x^3 - 10395*x

#' \code{k()} functions for Edgeworth expansions - one-sample
#' 
#' Calculate \code{k}'s (cumulant components) for a general version of Edgeworth
#' expansions (EE) for one-sample t-statistic.
#' 
#' Variance adjustment \eqn{r^2} is equal to the output of \code{K21one()},
#' unless different variance estimates are used for \code{A}, numerator of
#' \code{k}, and \code{r}.
#' 
#' @name kfuns1
#' @family \code{k()} functions
#'   
#' @param A value of \code{A} (depends on the type of the test).
#' @param B value of \code{B} (depends on the type of the test).
#' @param mu2,mu3,mu4,mu5,mu6 central moments (2 - 6) or their estimates.
#'   
#' @return A calculated value for a respective component.
#' @rdname kfuns1
K12one <- function(A, B, mu2, mu3, mu4, mu5, mu6) {
  -1/2*B*mu3/A^(3/2)
}
#' @rdname kfuns1
K13one <- function(A, B, mu2, mu3, mu4, mu5, mu6) {
  -1/16*(6*(8*mu2*mu3 - mu5)*A*B^2 - 15*(mu2^2*mu3 - mu3*mu4)*B^3 -
           8*A^2*B*mu3)/A^(7/2)
}
#' @rdname kfuns1
K21one <- function(A, B, mu2, mu3, mu4, mu5, mu6) {
  mu2/A
}	
#' @rdname kfuns1
K22one <- function(A, B, mu2, mu3, mu4, mu5, mu6) {
  1/4*(4*(4*mu2^2 - mu4)*A*B - (4*mu2^3 - 7*mu3^2 - 4*mu2*mu4)*B^2)/A^3
}
#' @rdname kfuns1
K23one <- function(A, B, mu2, mu3, mu4, mu5, mu6) {
  -1/16*(16*(3*mu2^2 - mu4)*A^3*B - 8*(58*mu2^3 - 19*mu3^2 - 30*mu2*mu4 +
                                         2*mu6)*A^2*B^2 + 2*(112*mu2^4 - 360*mu2*mu3^2 - 144*mu2^2*mu4 + 24*mu4^2 + 
                                                               45*mu3*mu5 + 8*mu2*mu6)*A*B^3 - 3*(16*mu2^5 - 59*mu2^2*mu3^2 +
                                                                                                    16*mu2*mu4^2 - (32*mu2^3 - 59*mu3^2)*mu4)*B^4)/A^5
}
#' @rdname kfuns1
K31one <- function(A, B, mu2, mu3, mu4, mu5, mu6) {
  -(3*B*mu2 - A)*mu3/A^(5/2)
}	
#' @rdname kfuns1
K32one <- function(A, B, mu2, mu3, mu4, mu5, mu6) {
  1/8*(12*(13*mu2*mu3 - mu5)*A^2*B - 3*(175*mu2^2*mu3 - 31*mu3*mu4 -
                                          12*mu2*mu5)*A*B^2 + (123*mu2^3*mu3 - 83*mu3^3 -
                                                                 123*mu2*mu3*mu4)*B^3)/A^(9/2)
}	
#' @rdname kfuns1
K41one <- function(A, B, mu2, mu3, mu4, mu5, mu6) {
  -((3*mu2^2 - mu4)*A^2 - 6*(3*mu2^3 - mu3^2 - mu2*mu4)*A*B + 3*(mu2^4 -
                                                                   6*mu2*mu3^2 - mu2^2*mu4)*B^2)/A^4
}
#' @rdname kfuns1
K42one <- function(A, B, mu2, mu3, mu4, mu5, mu6) {
  -1/8*(16*(42*mu2^3 - 13*mu3^2 - 19*mu2*mu4 + mu6)*A^3*B - 12*(284*mu2^4 - 
                                                                  266*mu2*mu3^2 - 172*mu2^2*mu4 + 12*mu4^2 + 21*mu3*mu5 +
                                                                  8*mu2*mu6)*A^2*B^2 + 12*(112*mu2^5 - 702*mu2^2*mu3^2 + 32*mu2*mu4^2 +
                                                                                             63*mu2*mu3*mu5 + 4*mu2^2*mu6 - 2*(74*mu2^3 - 51*mu3^2)*mu4)*A*B^3 -
          3*(64*mu2^6 - 624*mu2^3*mu3^2 + 233*mu3^4 + 64*mu2^2*mu4^2 - 16*(8*mu2^4 -
                                                                             39*mu2*mu3^2)*mu4)*B^4)/A^6
}	
#' @rdname kfuns1
K51one <- function(A, B, mu2, mu3, mu4, mu5, mu6) {
  -1/2*(2*(10*mu2*mu3 - mu5)*A^3 - 10*(35*mu2^2*mu3 - 5*mu3*mu4 -
                                         2*mu2*mu5)*A^2*B + 15*(56*mu2^3*mu3 - 7*mu3^3 - 16*mu2*mu3*mu4 -
                                                                  2*mu2^2*mu5)*A*B^2 - 15*(10*mu2^4*mu3 - 21*mu2*mu3^3 -
                                                                                             10*mu2^2*mu3*mu4)*B^3)/A^(11/2)
}	
#' @rdname kfuns1
K61one <- function(A, B, mu2, mu3, mu4, mu5, mu6) {
  1/2*(2*(30*mu2^3 - 10*mu3^2 - 15*mu2*mu4 + mu6)*A^4 - 30*(48*mu2^4 -
                                                              40*mu2*mu3^2 - 27*mu2^2*mu4 + 2*mu4^2 + 3*mu3*mu5 + mu2*mu6)*A^3*B +
         15*(336*mu2^5 - 623*mu2^2*mu3^2 + 28*mu2*mu4^2 + 44*mu2*mu3*mu5 +
               6*mu2^2*mu6 - (226*mu2^3 - 67*mu3^2)*mu4)*A^2*B^2 - 30*(56*mu2^6 -
                                                                         597*mu2^3*mu3^2 + 40*mu3^4 + 18*mu2^2*mu4^2 + 33*mu2^2*mu3*mu5 +
                                                                         mu2^3*mu6 - 3*(25*mu2^4 - 51*mu2*mu3^2)*mu4)*A*B^3 + 45*(4*mu2^7 -
                                                                                                                                    73*mu2^4*mu3^2 + 80*mu2*mu3^4 + 4*mu2^3*mu4^2 - (8*mu2^5 -
                                                                                                                                                                                       73*mu2^2*mu3^2)*mu4)*B^4)/A^7
}	

#' \code{k()} functions for Edgeworth expansions - two-sample 
#' 
#' Calculate \code{k}'s (cumulant components) for a general version of Edgeworth expansions (EE) for two-sample t-statistic.
#' 
#' Note that the test statistic for this Edgeworth expansion is defined as \eqn{\sqrt{n}(\bar{X} - \bar{Y})/s}{sqrt(n)(X-bar - Y-bar)/s} and therefore \code{X} would normally represent a treatment group and \code{Y} - control group. Variance adjustment \eqn{r^2} is equal to the output of \code{K21two()}, unless different variance estimates are used for \code{A}, numerator of \code{k}, and \code{r}.
#' 
#' @name kfuns2
#' @family \code{k()} functions
#' 
#' @param A value of \code{A}.
#' @param B_x value of \eqn{B_x}{B[x]} (depends on the type of the test).
#' @param B_y value of \eqn{B_y}{B[y]} (depends on the type of the test). 
#' @param b_x value of \eqn{b_x}{b[x]} - equal to \eqn{n/n_x}{n/n[x]}, where \eqn{n = (n_x + n_y)/2}{n = (n[x] + n[y])/2}.
#' @param b_y value of \eqn{b_y}{b[y]} - equal to \eqn{n/n_y}{n/n[y]}, where \eqn{n = (n_x + n_y)/2}{n = (n[x] + n[y])/2}.
#' @param mu_x2,mu_x3,mu_x4,mu_x5,mu_x6 central moments (2 - 6) for a treatment group or their estimates.
#' @param mu_y2,mu_y3,mu_y4,mu_y5,mu_y6 central moments (2 - 6) for a control group or their estimates. 
#' 
#' @return A calculated value for a respective component.
#' @rdname kfuns2
K12two <- function(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                   mu_y2, mu_y3, mu_y4, mu_y5, mu_y6) {
  -1/2*(B_x*b_x*mu_x3 - B_y*b_y*mu_y3)/A^(3/2)
}
#' @rdname kfuns2
K13two <- function(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                   mu_y2, mu_y3, mu_y4, mu_y5, mu_y6) {
  1/16*(15*B_x^3*b_x^2*mu_x2^2*mu_x3 + 15*B_x*B_y^2*b_x*b_y*mu_x3*mu_y2^2
        - 15*B_x^3*b_x^2*mu_x3*mu_x4 + 8*(B_x*b_x^2*mu_x3 - B_y*b_y^2*mu_y3)*A^2
        - 6*(8*B_x^2*b_x^2*mu_x2*mu_x3 + 2*B_x*B_y*b_x*b_y*mu_x3*mu_y2 -
               B_x^2*b_x^2*mu_x5 + B_y^2*b_y^2*mu_y5 - 2*(B_x*B_y*b_x*b_y*mu_x2 +
                                                            4*B_y^2*b_y^2*mu_y2)*mu_y3)*A - 15*(B_x^2*B_y*b_x*b_y*mu_x2^2 +
                                                                                                  B_y^3*b_y^2*mu_y2^2 - B_x^2*B_y*b_x*b_y*mu_x4)*mu_y3 -
          15*(B_x*B_y^2*b_x*b_y*mu_x3 - B_y^3*b_y^2*mu_y3)*mu_y4)/A^(7/2)
}
#' @rdname kfuns2
K21two <- function(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                   mu_y2, mu_y3, mu_y4, mu_y5, mu_y6) {
  (b_x*mu_x2 + b_y*mu_y2)/A
}
#' @rdname kfuns2
K22two <- function(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                   mu_y2, mu_y3, mu_y4, mu_y5, mu_y6) {
  -1/4*(4*B_x^2*b_x^2*mu_x2^3 + 4*B_y^2*b_x*b_y*mu_x2*mu_y2^2 +
          4*B_y^2*b_y^2*mu_y2^3 - 7*B_x^2*b_x^2*mu_x3^2 -
          4*B_x^2*b_x^2*mu_x2*mu_x4 + 14*B_x*B_y*b_x*b_y*mu_x3*mu_y3 -
          7*B_y^2*b_y^2*mu_y3^2 - 4*(4*B_x*b_x^2*mu_x2^2 + 4*B_y*b_y^2*mu_y2^2 -
                                       B_x*b_x^2*mu_x4 - B_y*b_y^2*mu_y4 + (B_x*b_x*b_y +
                                                                              B_y*b_x*b_y)*mu_x2*mu_y2)*A + 4*(B_x^2*b_x*b_y*mu_x2^2 -
                                                                                                                 B_x^2*b_x*b_y*mu_x4)*mu_y2 - 4*(B_y^2*b_x*b_y*mu_x2 +
                                                                                                                                                   B_y^2*b_y^2*mu_y2)*mu_y4)/A^3
}
#' @rdname kfuns2
K23two <- function(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                   mu_y2, mu_y3, mu_y4, mu_y5, mu_y6) {
  1/16*(48*B_x^4*b_x^3*mu_x2^5 + 48*B_y^4*b_x*b_y^2*mu_x2*mu_y2^4 +
          48*B_y^4*b_y^3*mu_y2^5 - 177*B_x^4*b_x^3*mu_x2^2*mu_x3^2 +
          48*B_x^4*b_x^3*mu_x2*mu_x4^2 - 16*(3*B_x*b_x^3*mu_x2^2 +
                                               3*B_y*b_y^3*mu_y2^2 - B_x*b_x^3*mu_x4 - B_y*b_y^3*mu_y4)*A^3 +
          96*(B_x^2*B_y^2*b_x*b_y^2*mu_x2^2 - B_x^2*B_y^2*b_x*b_y^2*mu_x4)*mu_y2^3
        + 8*(58*B_x^2*b_x^3*mu_x2^3 + 58*B_y^2*b_y^3*mu_y2^3 -
               19*B_x^2*b_x^3*mu_x3^2 - 30*B_x^2*b_x^3*mu_x2*mu_x4 -
               19*B_y^2*b_y^3*mu_y3^2 + 2*B_x^2*b_x^3*mu_x6 + 2*B_y^2*b_y^3*mu_y6 +
               7*(b_x^2*b_y + b_x*b_y^2)*B_x*B_y*mu_x3*mu_y3 + 2*(8*B_x*B_y*b_x*b_y^2 +
                                                                    5*B_y^2*b_x*b_y^2)*mu_x2*mu_y2^2 + 2*((5*B_x^2*b_x^2*b_y +
                                                                                                             8*B_x*B_y*b_x^2*b_y)*mu_x2^2 - 2*(B_x^2*b_x^2*b_y +
                                                                                                                                                 B_x*B_y*b_x^2*b_y)*mu_x4)*mu_y2 - 2*(15*B_y^2*b_y^3*mu_y2 +
                                                                                                                                                                                        2*(B_x*B_y*b_x*b_y^2 + B_y^2*b_x*b_y^2)*mu_x2)*mu_y4)*A^2 +
          3*(32*B_x^2*B_y^2*b_x^2*b_y*mu_x2^3 - 59*B_x^2*B_y^2*b_x^2*b_y*mu_x3^2 -
               32*B_x^2*B_y^2*b_x^2*b_y*mu_x2*mu_x4)*mu_y2^2 -
          177*(B_x^2*B_y^2*b_x*b_y^2*mu_x2^2 + B_y^4*b_y^3*mu_y2^2 -
                 B_x^2*B_y^2*b_x*b_y^2*mu_x4)*mu_y3^2 + 48*(B_y^4*b_x*b_y^2*mu_x2 +
                                                              B_y^4*b_y^3*mu_y2)*mu_y4^2 - 2*(112*B_x^3*b_x^3*mu_x2^4 +
                                                                                                112*B_y^3*b_y^3*mu_y2^4 - 360*B_x^3*b_x^3*mu_x2*mu_x3^2 -
                                                                                                144*B_x^3*b_x^3*mu_x2^2*mu_x4 + 24*B_x^3*b_x^3*mu_x4^2 +
                                                                                                45*B_x^3*b_x^3*mu_x3*mu_x5 + 8*B_x^3*b_x^3*mu_x2*mu_x6 +
                                                                                                24*B_y^3*b_y^3*mu_y4^2 + 8*(3*B_x*B_y^2*b_x*b_y^2 +
                                                                                                                              5*B_y^3*b_x*b_y^2)*mu_x2*mu_y2^3 + 24*(4*(B_x*B_y^2*b_x^2*b_y +
                                                                                                                                                                          B_x^2*B_y*b_x*b_y^2)*mu_x2^2 - (B_x*B_y^2*b_x^2*b_y +
                                                                                                                                                                                                            4*B_x^2*B_y*b_x*b_y^2)*mu_x4)*mu_y2^2 - 6*(60*B_y^3*b_y^3*mu_y2 +
                                                                                                                                                                                                                                                         (7*B_x*B_y^2*b_x*b_y^2 + 8*B_y^3*b_x*b_y^2)*mu_x2)*mu_y3^2 +
                                                                                                2*(4*B_x^3*b_x^2*b_y*mu_x6 + 4*(5*B_x^3*b_x^2*b_y +
                                                                                                                                  3*B_x^2*B_y*b_x^2*b_y)*mu_x2^3 - 3*(8*B_x^3*b_x^2*b_y +
                                                                                                                                                                        7*B_x^2*B_y*b_x^2*b_y)*mu_x3^2 - 12*(2*B_x^3*b_x^2*b_y +
                                                                                                                                                                                                               B_x^2*B_y*b_x^2*b_y)*mu_x2*mu_x4)*mu_y2 +
                                                                                                3*(118*B_x^2*B_y*b_x^2*b_y*mu_x2*mu_x3 +
                                                                                                     118*B_x*B_y^2*b_x*b_y^2*mu_x3*mu_y2 -
                                                                                                     15*B_x^2*B_y*b_x^2*b_y*mu_x5)*mu_y3 - 24*(6*B_y^3*b_y^3*mu_y2^2 +
                                                                                                                                                 (4*B_x*B_y^2*b_x^2*b_y + B_x^2*B_y*b_x*b_y^2)*mu_x2^2 +
                                                                                                                                                 (B_x*B_y^2*b_x*b_y^2 + 2*B_y^3*b_x*b_y^2)*mu_x2*mu_y2 -
                                                                                                                                                 (B_x*B_y^2*b_x^2*b_y + B_x^2*B_y*b_x*b_y^2)*mu_x4)*mu_y4 -
                                                                                                45*(B_x*B_y^2*b_x*b_y^2*mu_x3 - B_y^3*b_y^3*mu_y3)*mu_y5 +
                                                                                                8*(B_y^3*b_x*b_y^2*mu_x2 + B_y^3*b_y^3*mu_y2)*mu_y6)*A -
          3*(32*B_x^4*b_x^3*mu_x2^3 - 59*B_x^4*b_x^3*mu_x3^2)*mu_x4 +
          48*(B_x^4*b_x^2*b_y*mu_x2^4 - 2*B_x^4*b_x^2*b_y*mu_x2^2*mu_x4 +
                B_x^4*b_x^2*b_y*mu_x4^2)*mu_y2 + 354*(B_x^3*B_y*b_x^2*b_y*mu_x2^2*mu_x3
                                                      + B_x*B_y^3*b_x*b_y^2*mu_x3*mu_y2^2 -
                                                        B_x^3*B_y*b_x^2*b_y*mu_x3*mu_x4)*mu_y3 -
          3*(32*B_x^2*B_y^2*b_x^2*b_y*mu_x2^3 + 32*B_y^4*b_x*b_y^2*mu_x2*mu_y2^2 +
               32*B_y^4*b_y^3*mu_y2^3 - 59*B_x^2*B_y^2*b_x^2*b_y*mu_x3^2 -
               32*B_x^2*B_y^2*b_x^2*b_y*mu_x2*mu_x4 +
               118*B_x*B_y^3*b_x*b_y^2*mu_x3*mu_y3 - 59*B_y^4*b_y^3*mu_y3^2 +
               32*(B_x^2*B_y^2*b_x*b_y^2*mu_x2^2 -
                     B_x^2*B_y^2*b_x*b_y^2*mu_x4)*mu_y2)*mu_y4)/A^5
}
#' @rdname kfuns2
K31two <- function(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                   mu_y2, mu_y3, mu_y4, mu_y5, mu_y6) {
  -(3*B_x*b_x^2*mu_x2*mu_x3 + 3*B_x*b_x*b_y*mu_x3*mu_y2 - (b_x^2*mu_x3 -
                                                             b_y^2*mu_y3)*A - 3*(B_y*b_x*b_y*mu_x2 + B_y*b_y^2*mu_y2)*mu_y3)/A^(5/2)
}
K32two <- function(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                   mu_y2, mu_y3, mu_y4, mu_y5, mu_y6) {
  1/8*(123*B_x^3*b_x^3*mu_x2^3*mu_x3 +
         123*B_x*B_y^2*b_x^2*b_y*mu_x2*mu_x3*mu_y2^2 +
         123*B_x*B_y^2*b_x*b_y^2*mu_x3*mu_y2^3 - 83*B_x^3*b_x^3*mu_x3^3 -
         123*B_x^3*b_x^3*mu_x2*mu_x3*mu_x4 -
         249*B_x*B_y^2*b_x*b_y^2*mu_x3*mu_y3^2 + 83*B_y^3*b_y^3*mu_y3^3 +
         12*(13*B_x*b_x^3*mu_x2*mu_x3 - B_x*b_x^3*mu_x5 + B_y*b_y^3*mu_y5 +
               (2*B_x*b_x^2*b_y + B_y*b_x^2*b_y)*mu_x3*mu_y2 - (13*B_y*b_y^3*mu_y2 +
                                                                  (B_x*b_x*b_y^2 + 2*B_y*b_x*b_y^2)*mu_x2)*mu_y3)*A^2 -
         3*(175*B_x^2*b_x^3*mu_x2^2*mu_x3 - 31*B_x^2*b_x^3*mu_x3*mu_x4 -
              12*B_x^2*b_x^3*mu_x2*mu_x5 + (5*B_y^2*b_x^2*b_y +
                                              98*B_x*B_y*b_x*b_y^2)*mu_x3*mu_y2^2 - 4*(3*B_x^2*b_x^2*b_y*mu_x5 -
                                                                                         (23*B_x^2*b_x^2*b_y + 5*B_x*B_y*b_x^2*b_y)*mu_x2*mu_x3)*mu_y2 -
              (175*B_y^2*b_y^3*mu_y2^2 + (98*B_x*B_y*b_x^2*b_y +
                                            5*B_x^2*b_x*b_y^2)*mu_x2^2 + 4*(5*B_x*B_y*b_x*b_y^2 +
                                                                              23*B_y^2*b_x*b_y^2)*mu_x2*mu_y2 - (26*B_x*B_y*b_x^2*b_y +
                                                                                                                   5*B_x^2*b_x*b_y^2)*mu_x4)*mu_y3 + (31*B_y^2*b_y^3*mu_y3 -
                                                                                                                                                        (5*B_y^2*b_x^2*b_y + 26*B_x*B_y*b_x*b_y^2)*mu_x3)*mu_y4 +
              12*(B_y^2*b_x*b_y^2*mu_x2 + B_y^2*b_y^3*mu_y2)*mu_y5)*A +
         123*(B_x^3*b_x^2*b_y*mu_x2^2*mu_x3 - B_x^3*b_x^2*b_y*mu_x3*mu_x4)*mu_y2
       - 3*(41*B_x^2*B_y*b_x^2*b_y*mu_x2^3 + 41*B_y^3*b_x*b_y^2*mu_x2*mu_y2^2 +
              41*B_y^3*b_y^3*mu_y2^3 - 83*B_x^2*B_y*b_x^2*b_y*mu_x3^2 -
              41*B_x^2*B_y*b_x^2*b_y*mu_x2*mu_x4 + 41*(B_x^2*B_y*b_x*b_y^2*mu_x2^2 -
                                                         B_x^2*B_y*b_x*b_y^2*mu_x4)*mu_y2)*mu_y3 -
         123*(B_x*B_y^2*b_x^2*b_y*mu_x2*mu_x3 + B_x*B_y^2*b_x*b_y^2*mu_x3*mu_y2 -
                (B_y^3*b_x*b_y^2*mu_x2 + B_y^3*b_y^3*mu_y2)*mu_y3)*mu_y4)/A^(9/2)
}
#' @rdname kfuns2
K41two <- function(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                   mu_y2, mu_y3, mu_y4, mu_y5, mu_y6) {
  -(3*B_x^2*b_x^3*mu_x2^4 + 6*B_y^2*b_x*b_y^2*mu_x2*mu_y2^3 +
      3*B_y^2*b_y^3*mu_y2^4 - 18*B_x^2*b_x^3*mu_x2*mu_x3^2 -
      3*B_x^2*b_x^3*mu_x2^2*mu_x4 + (3*b_x^3*mu_x2^2 + 3*b_y^3*mu_y2^2 -
                                       b_x^3*mu_x4 - b_y^3*mu_y4)*A^2 - 3*(B_x^2*b_x*b_y^2*mu_x4 -
                                                                             (B_y^2*b_x^2*b_y + B_x^2*b_x*b_y^2)*mu_x2^2)*mu_y2^2 -
      18*(B_y^2*b_x*b_y^2*mu_x2 + B_y^2*b_y^3*mu_y2)*mu_y3^2 -
      6*(3*B_x*b_x^3*mu_x2^3 + 3*B_y*b_x*b_y^2*mu_x2*mu_y2^2 +
           3*B_y*b_y^3*mu_y2^3 - B_x*b_x^3*mu_x3^2 - B_x*b_x^3*mu_x2*mu_x4 -
           B_y*b_y^3*mu_y3^2 + (B_y*b_x^2*b_y + B_x*b_x*b_y^2)*mu_x3*mu_y3 +
           (3*B_x*b_x^2*b_y*mu_x2^2 - B_x*b_x^2*b_y*mu_x4)*mu_y2 -
           (B_y*b_x*b_y^2*mu_x2 + B_y*b_y^3*mu_y2)*mu_y4)*A +
      6*(B_x^2*b_x^2*b_y*mu_x2^3 - 3*B_x^2*b_x^2*b_y*mu_x3^2 -
           B_x^2*b_x^2*b_y*mu_x2*mu_x4)*mu_y2 + 36*(B_x*B_y*b_x^2*b_y*mu_x2*mu_x3 +
                                                      B_x*B_y*b_x*b_y^2*mu_x3*mu_y2)*mu_y3 - 3*(B_y^2*b_x^2*b_y*mu_x2^2 +
                                                                                                  2*B_y^2*b_x*b_y^2*mu_x2*mu_y2 + B_y^2*b_y^3*mu_y2^2)*mu_y4)/A^4
}
#' @rdname kfuns2
K42two <- function(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                   mu_y2, mu_y3, mu_y4, mu_y5, mu_y6) {
  1/8*(192*B_x^4*b_x^4*mu_x2^6 + 384*B_y^4*b_x*b_y^3*mu_x2*mu_y2^5 +
         192*B_y^4*b_y^4*mu_y2^6 - 1872*B_x^4*b_x^4*mu_x2^3*mu_x3^2 +
         699*B_x^4*b_x^4*mu_x3^4 + 192*B_x^4*b_x^4*mu_x2^2*mu_x4^2 -
         2796*B_x*B_y^3*b_x*b_y^3*mu_x3*mu_y3^3 + 699*B_y^4*b_y^4*mu_y3^4 -
         192*(2*B_x^2*B_y^2*b_x*b_y^3*mu_x4 - (B_y^4*b_x^2*b_y^2 +
                                                 2*B_x^2*B_y^2*b_x*b_y^3)*mu_x2^2)*mu_y2^4 - 16*(42*B_x*b_x^4*mu_x2^3 +
                                                                                                   42*B_y*b_y^4*mu_y2^3 - 13*B_x*b_x^4*mu_x3^2 - 19*B_x*b_x^4*mu_x2*mu_x4 -
                                                                                                   13*B_y*b_y^4*mu_y3^2 + B_x*b_x^4*mu_x6 + B_y*b_y^4*mu_y6 +
                                                                                                   3*(B_x*b_x*b_y^3 + 3*B_y*b_x*b_y^3)*mu_x2*mu_y2^2 + 3*(B_x*b_x^2*b_y^2 +
                                                                                                                                                            B_y*b_x^2*b_y^2)*mu_x3*mu_y3 + (3*(3*B_x*b_x^3*b_y +
                                                                                                                                                                                                 B_y*b_x^3*b_y)*mu_x2^2 - (3*B_x*b_x^3*b_y + B_y*b_x^3*b_y)*mu_x4)*mu_y2
                                                                                                 - (19*B_y*b_y^4*mu_y2 + (B_x*b_x*b_y^3 +
                                                                                                                            3*B_y*b_x*b_y^3)*mu_x2)*mu_y4)*A^3 +
         48*(16*B_x^2*B_y^2*b_x^2*b_y^2*mu_x2^3 -
               39*B_x^2*B_y^2*b_x^2*b_y^2*mu_x3^2 -
               16*B_x^2*B_y^2*b_x^2*b_y^2*mu_x2*mu_x4)*mu_y2^3 +
         12*(284*B_x^2*b_x^4*mu_x2^4 + 284*B_y^2*b_y^4*mu_y2^4 -
               266*B_x^2*b_x^4*mu_x2*mu_x3^2 - 172*B_x^2*b_x^4*mu_x2^2*mu_x4 +
               12*B_x^2*b_x^4*mu_x4^2 + 21*B_x^2*b_x^4*mu_x3*mu_x5 +
               8*B_x^2*b_x^4*mu_x2*mu_x6 + 12*B_y^2*b_y^4*mu_y4^2 +
               4*(9*B_x*B_y*b_x*b_y^3 + 49*B_y^2*b_x*b_y^3)*mu_x2*mu_y2^3 +
               2*((126*B_x*B_y*b_x^2*b_y^2 + (4*b_x^2*b_y^2 + 3*b_x*b_y^3)*B_x^2 +
                     (3*b_x^3*b_y + 4*b_x^2*b_y^2)*B_y^2)*mu_x2^2 - (B_y^2*b_x^3*b_y +
                                                                       36*B_x*B_y*b_x^2*b_y^2 + (2*b_x^2*b_y^2 +
                                                                                                   3*b_x*b_y^3)*B_x^2)*mu_x4)*mu_y2^2 - 2*(133*B_y^2*b_y^4*mu_y2 +
                                                                                                                                             6*(B_x*B_y*b_x*b_y^3 + 6*B_y^2*b_x*b_y^3)*mu_x2)*mu_y3^2 +
               4*(2*B_x^2*b_x^3*b_y*mu_x6 + (49*B_x^2*b_x^3*b_y +
                                               9*B_x*B_y*b_x^3*b_y)*mu_x2^3 - 3*(6*B_x^2*b_x^3*b_y +
                                                                                   B_x*B_y*b_x^3*b_y)*mu_x3^2 - 3*(9*B_x^2*b_x^3*b_y +
                                                                                                                     B_x*B_y*b_x^3*b_y)*mu_x2*mu_x4)*mu_y2 + (2*(27*B_x^2*b_x^2*b_y^2 +
                                                                                                                                                                   4*(22*b_x^3*b_y + 3*b_x^2*b_y^2)*B_x*B_y)*mu_x2*mu_x3 +
                                                                                                                                                                2*(27*B_y^2*b_x^2*b_y^2 + 4*(3*b_x^2*b_y^2 +
                                                                                                                                                                                               22*b_x*b_y^3)*B_x*B_y)*mu_x3*mu_y2 - 7*(2*B_x*B_y*b_x^3*b_y +
                                                                                                                                                                                                                                         B_x^2*b_x^2*b_y^2)*mu_x5)*mu_y3 - 2*(86*B_y^2*b_y^4*mu_y2^2 +
                                                                                                                                                                                                                                                                                (36*B_x*B_y*b_x^2*b_y^2 + B_x^2*b_x*b_y^3 + (3*b_x^3*b_y +
                                                                                                                                                                                                                                                                                                                               2*b_x^2*b_y^2)*B_y^2)*mu_x2^2 + 6*(B_x*B_y*b_x*b_y^3 +
                                                                                                                                                                                                                                                                                                                                                                    9*B_y^2*b_x*b_y^3)*mu_x2*mu_y2 - (B_y^2*b_x^3*b_y +
                                                                                                                                                                                                                                                                                                                                                                                                        10*B_x*B_y*b_x^2*b_y^2 + B_x^2*b_x*b_y^3)*mu_x4)*mu_y4 +
               7*(3*B_y^2*b_y^4*mu_y3 - (B_y^2*b_x^2*b_y^2 +
                                           2*B_x*B_y*b_x*b_y^3)*mu_x3)*mu_y5 + 8*(B_y^2*b_x*b_y^3*mu_x2 +
                                                                                    B_y^2*b_y^4*mu_y2)*mu_y6)*A^2 -
         48*(39*B_x^2*B_y^2*b_x^3*b_y*mu_x2*mu_x3^2 - 4*B_x^4*b_x^2*b_y^2*mu_x4^2
             - 4*(2*B_x^2*B_y^2*b_x^3*b_y + B_x^4*b_x^2*b_y^2)*mu_x2^4 +
               8*(B_x^2*B_y^2*b_x^3*b_y + B_x^4*b_x^2*b_y^2)*mu_x2^2*mu_x4)*mu_y2^2 -
         18*(104*B_x^2*B_y^2*b_x^2*b_y^2*mu_x2^3 +
               104*B_y^4*b_x*b_y^3*mu_x2*mu_y2^2 + 104*B_y^4*b_y^4*mu_y2^3 -
               233*B_x^2*B_y^2*b_x^2*b_y^2*mu_x3^2 -
               104*B_x^2*B_y^2*b_x^2*b_y^2*mu_x2*mu_x4 +
               104*(B_x^2*B_y^2*b_x*b_y^3*mu_x2^2 -
                      B_x^2*B_y^2*b_x*b_y^3*mu_x4)*mu_y2)*mu_y3^2 +
         192*(B_y^4*b_x^2*b_y^2*mu_x2^2 + 2*B_y^4*b_x*b_y^3*mu_x2*mu_y2 +
                B_y^4*b_y^4*mu_y2^2)*mu_y4^2 - 12*(112*B_x^3*b_x^4*mu_x2^5 +
                                                     112*B_y^3*b_y^4*mu_y2^5 - 702*B_x^3*b_x^4*mu_x2^2*mu_x3^2 +
                                                     32*B_x^3*b_x^4*mu_x2*mu_x4^2 + 63*B_x^3*b_x^4*mu_x2*mu_x3*mu_x5 +
                                                     4*B_x^3*b_x^4*mu_x2^2*mu_x6 + 8*(B_x*B_y^2*b_x*b_y^3 +
                                                                                        16*B_y^3*b_x*b_y^3)*mu_x2*mu_y2^4 + 8*((14*B_x*B_y^2*b_x^2*b_y^2 +
                                                                                                                                  2*B_y^3*b_x^2*b_y^2 + 13*B_x^2*B_y*b_x*b_y^3)*mu_x2^2 -
                                                                                                                                 (4*B_x*B_y^2*b_x^2*b_y^2 + 13*B_x^2*B_y*b_x*b_y^3)*mu_x4)*mu_y2^3 +
                                                     (4*B_x^3*b_x^2*b_y^2*mu_x6 + 8*(13*B_x*B_y^2*b_x^3*b_y +
                                                                                       2*B_x^3*b_x^2*b_y^2 + 14*B_x^2*B_y*b_x^2*b_y^2)*mu_x2^3 -
                                                        3*(9*B_x*B_y^2*b_x^3*b_y + 8*B_x^3*b_x^2*b_y^2 +
                                                             91*B_x^2*B_y*b_x^2*b_y^2)*mu_x3^2 - 4*(8*B_x*B_y^2*b_x^3*b_y +
                                                                                                      5*B_x^3*b_x^2*b_y^2 + 28*B_x^2*B_y*b_x^2*b_y^2)*mu_x2*mu_x4)*mu_y2^2 -
                                                     3*(234*B_y^3*b_y^4*mu_y2^2 + (91*B_x*B_y^2*b_x^2*b_y^2 +
                                                                                     8*B_y^3*b_x^2*b_y^2 + 9*B_x^2*B_y*b_x*b_y^3)*mu_x2^2 +
                                                          2*(8*B_x*B_y^2*b_x*b_y^3 + 79*B_y^3*b_x*b_y^3)*mu_x2*mu_y2 -
                                                          (25*B_x*B_y^2*b_x^2*b_y^2 + 9*B_x^2*B_y*b_x*b_y^3)*mu_x4)*mu_y3^2 +
                                                     32*(B_y^3*b_x*b_y^3*mu_x2 + B_y^3*b_y^4*mu_y2)*mu_y4^2 -
                                                     2*(74*B_x^3*b_x^4*mu_x2^3 - 51*B_x^3*b_x^4*mu_x3^2)*mu_x4 +
                                                     (32*B_x^3*b_x^3*b_y*mu_x4^2 + 63*B_x^3*b_x^3*b_y*mu_x3*mu_x5 +
                                                        8*B_x^3*b_x^3*b_y*mu_x2*mu_x6 + 8*(16*B_x^3*b_x^3*b_y +
                                                                                             B_x^2*B_y*b_x^3*b_y)*mu_x2^4 - 6*(79*B_x^3*b_x^3*b_y +
                                                                                                                                 8*B_x^2*B_y*b_x^3*b_y)*mu_x2*mu_x3^2 - 8*(21*B_x^3*b_x^3*b_y +
                                                                                                                                                                             B_x^2*B_y*b_x^3*b_y)*mu_x2^2*mu_x4)*mu_y2 -
                                                     3*(21*B_x^2*B_y*b_x^3*b_y*mu_x2*mu_x5 - (317*B_x^2*B_y*b_x^3*b_y +
                                                                                                9*B_x^3*b_x^2*b_y^2)*mu_x2^2*mu_x3 - (9*B_y^3*b_x^2*b_y^2 +
                                                                                                                                        317*B_x*B_y^2*b_x*b_y^3)*mu_x3*mu_y2^2 + (59*B_x^2*B_y*b_x^3*b_y +
                                                                                                                                                                                    9*B_x^3*b_x^2*b_y^2)*mu_x3*mu_x4 + (21*B_x^2*B_y*b_x^2*b_y^2*mu_x5 -
                                                                                                                                                                                                                          158*(B_x^2*B_y*b_x^2*b_y^2 +
                                                                                                                                                                                                                                 B_x*B_y^2*b_x^2*b_y^2)*mu_x2*mu_x3)*mu_y2)*mu_y3 -
                                                     (148*B_y^3*b_y^4*mu_y2^3 - 102*B_y^3*b_y^4*mu_y3^2 +
                                                        8*(13*B_x*B_y^2*b_x^3*b_y + 4*B_x^2*B_y*b_x^2*b_y^2)*mu_x2^3 +
                                                        8*(B_x*B_y^2*b_x*b_y^3 + 21*B_y^3*b_x*b_y^3)*mu_x2*mu_y2^2 -
                                                        3*(9*B_x*B_y^2*b_x^3*b_y + 25*B_x^2*B_y*b_x^2*b_y^2)*mu_x3^2 -
                                                        32*(B_x*B_y^2*b_x^3*b_y + B_x^2*B_y*b_x^2*b_y^2)*mu_x2*mu_x4 +
                                                        3*(9*B_y^3*b_x^2*b_y^2 + 59*B_x*B_y^2*b_x*b_y^3)*mu_x3*mu_y3 +
                                                        4*((28*B_x*B_y^2*b_x^2*b_y^2 + 5*B_y^3*b_x^2*b_y^2 +
                                                              8*B_x^2*B_y*b_x*b_y^3)*mu_x2^2 - 8*(B_x*B_y^2*b_x^2*b_y^2 +
                                                                                                    B_x^2*B_y*b_x*b_y^3)*mu_x4)*mu_y2)*mu_y4 -
                                                     63*(B_x*B_y^2*b_x^2*b_y^2*mu_x2*mu_x3 + B_x*B_y^2*b_x*b_y^3*mu_x3*mu_y2
                                                         - (B_y^3*b_x*b_y^3*mu_x2 + B_y^3*b_y^4*mu_y2)*mu_y3)*mu_y5 +
                                                     4*(B_y^3*b_x^2*b_y^2*mu_x2^2 + 2*B_y^3*b_x*b_y^3*mu_x2*mu_y2 +
                                                          B_y^3*b_y^4*mu_y2^2)*mu_y6)*A - 48*(8*B_x^4*b_x^4*mu_x2^4 -
                                                                                                39*B_x^4*b_x^4*mu_x2*mu_x3^2)*mu_x4 + 48*(8*B_x^4*b_x^3*b_y*mu_x2^5 -
                                                                                                                                            39*B_x^4*b_x^3*b_y*mu_x2^2*mu_x3^2 + 8*B_x^4*b_x^3*b_y*mu_x2*mu_x4^2 -
                                                                                                                                            (16*B_x^4*b_x^3*b_y*mu_x2^3 - 39*B_x^4*b_x^3*b_y*mu_x3^2)*mu_x4)*mu_y2 +
         12*(312*B_x^3*B_y*b_x^3*b_y*mu_x2^3*mu_x3 +
               312*B_x*B_y^3*b_x^2*b_y^2*mu_x2*mu_x3*mu_y2^2 +
               312*B_x*B_y^3*b_x*b_y^3*mu_x3*mu_y2^3 - 233*B_x^3*B_y*b_x^3*b_y*mu_x3^3
             - 312*B_x^3*B_y*b_x^3*b_y*mu_x2*mu_x3*mu_x4 +
               312*(B_x^3*B_y*b_x^2*b_y^2*mu_x2^2*mu_x3 -
                      B_x^3*B_y*b_x^2*b_y^2*mu_x3*mu_x4)*mu_y2)*mu_y3 -
         48*(8*B_x^2*B_y^2*b_x^3*b_y*mu_x2^4 + 16*B_y^4*b_x*b_y^3*mu_x2*mu_y2^3 +
               8*B_y^4*b_y^4*mu_y2^4 - 39*B_x^2*B_y^2*b_x^3*b_y*mu_x2*mu_x3^2 -
               8*B_x^2*B_y^2*b_x^3*b_y*mu_x2^2*mu_x4 - 8*(B_x^2*B_y^2*b_x*b_y^3*mu_x4 -
                                                            (B_y^4*b_x^2*b_y^2 + B_x^2*B_y^2*b_x*b_y^3)*mu_x2^2)*mu_y2^2 -
               39*(B_y^4*b_x*b_y^3*mu_x2 + B_y^4*b_y^4*mu_y2)*mu_y3^2 +
               (16*B_x^2*B_y^2*b_x^2*b_y^2*mu_x2^3 - 39*B_x^2*B_y^2*b_x^2*b_y^2*mu_x3^2
                - 16*B_x^2*B_y^2*b_x^2*b_y^2*mu_x2*mu_x4)*mu_y2 +
               78*(B_x*B_y^3*b_x^2*b_y^2*mu_x2*mu_x3 +
                     B_x*B_y^3*b_x*b_y^3*mu_x3*mu_y2)*mu_y3)*mu_y4)/A^6
}
#' @rdname kfuns2
K51two <- function(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                   mu_y2, mu_y3, mu_y4, mu_y5, mu_y6) {
  1/2*(150*B_x^3*b_x^4*mu_x2^4*mu_x3 +
         300*B_x*B_y^2*b_x^2*b_y^2*mu_x2*mu_x3*mu_y2^3 +
         150*B_x*B_y^2*b_x*b_y^3*mu_x3*mu_y2^4 - 315*B_x^3*b_x^4*mu_x2*mu_x3^3 -
         150*B_x^3*b_x^4*mu_x2^2*mu_x3*mu_x4 - 2*(10*b_x^4*mu_x2*mu_x3 -
                                                    10*b_y^4*mu_y2*mu_y3 - b_x^4*mu_x5 + b_y^4*mu_y5)*A^3 +
         315*(B_y^3*b_x*b_y^3*mu_x2 + B_y^3*b_y^4*mu_y2)*mu_y3^3 +
         10*(35*B_x*b_x^4*mu_x2^2*mu_x3 - 5*B_x*b_x^4*mu_x3*mu_x4 -
               2*B_x*b_x^4*mu_x2*mu_x5 + 3*(3*B_y*b_x^2*b_y^2 +
                                              2*B_x*b_x*b_y^3)*mu_x3*mu_y2^2 + 2*(10*B_x*b_x^3*b_y*mu_x2*mu_x3 -
                                                                                    B_x*b_x^3*b_y*mu_x5)*mu_y2 - (20*B_y*b_x*b_y^3*mu_x2*mu_y2 +
                                                                                                                    35*B_y*b_y^4*mu_y2^2 + 3*(2*B_y*b_x^3*b_y + 3*B_x*b_x^2*b_y^2)*mu_x2^2 -
                                                                                                                    (2*B_y*b_x^3*b_y + 3*B_x*b_x^2*b_y^2)*mu_x4)*mu_y3 + (5*B_y*b_y^4*mu_y3
                                                                                                                                                                          - (3*B_y*b_x^2*b_y^2 + 2*B_x*b_x*b_y^3)*mu_x3)*mu_y4 +
               2*(B_y*b_x*b_y^3*mu_x2 + B_y*b_y^4*mu_y2)*mu_y5)*A^2 -
         150*(B_x^3*b_x^2*b_y^2*mu_x3*mu_x4 - (B_x*B_y^2*b_x^3*b_y +
                                                 B_x^3*b_x^2*b_y^2)*mu_x2^2*mu_x3)*mu_y2^2 -
         945*(B_x*B_y^2*b_x^2*b_y^2*mu_x2*mu_x3 +
                B_x*B_y^2*b_x*b_y^3*mu_x3*mu_y2)*mu_y3^2 -
         15*(56*B_x^2*b_x^4*mu_x2^3*mu_x3 - 7*B_x^2*b_x^4*mu_x3^3 -
               16*B_x^2*b_x^4*mu_x2*mu_x3*mu_x4 - 2*B_x^2*b_x^4*mu_x2^2*mu_x5 +
               7*B_y^2*b_y^4*mu_y3^3 + 2*(B_y^2*b_x^2*b_y^2 +
                                            21*B_x*B_y*b_x*b_y^3)*mu_x3*mu_y2^3 - 7*(B_y^2*b_x^2*b_y^2 +
                                                                                       2*B_x*B_y*b_x*b_y^3)*mu_x3*mu_y3^2 - 2*(B_x^2*b_x^2*b_y^2*mu_x5 -
                                                                                                                                 (B_y^2*b_x^3*b_y + 6*B_x^2*b_x^2*b_y^2 +
                                                                                                                                    21*B_x*B_y*b_x^2*b_y^2)*mu_x2*mu_x3)*mu_y2^2 +
               4*(17*B_x^2*b_x^3*b_y*mu_x2^2*mu_x3 - 4*B_x^2*b_x^3*b_y*mu_x3*mu_x4 -
                    B_x^2*b_x^3*b_y*mu_x2*mu_x5)*mu_y2 - (68*B_y^2*b_x*b_y^3*mu_x2*mu_y2^2 +
                                                            56*B_y^2*b_y^4*mu_y2^3 + 2*(21*B_x*B_y*b_x^3*b_y +
                                                                                          B_x^2*b_x^2*b_y^2)*mu_x2^3 - 7*(2*B_x*B_y*b_x^3*b_y +
                                                                                                                            B_x^2*b_x^2*b_y^2)*mu_x3^2 - 2*(7*B_x*B_y*b_x^3*b_y +
                                                                                                                                                              B_x^2*b_x^2*b_y^2)*mu_x2*mu_x4 + 2*((21*B_x*B_y*b_x^2*b_y^2 +
                                                                                                                                                                                                     6*B_y^2*b_x^2*b_y^2 + B_x^2*b_x*b_y^3)*mu_x2^2 - (7*B_x*B_y*b_x^2*b_y^2
                                                                                                                                                                                                                                                       + B_x^2*b_x*b_y^3)*mu_x4)*mu_y2)*mu_y3 - 2*((B_y^2*b_x^3*b_y +
                                                                                                                                                                                                                                                                                                      7*B_x*B_y*b_x^2*b_y^2)*mu_x2*mu_x3 + (B_y^2*b_x^2*b_y^2 +
                                                                                                                                                                                                                                                                                                                                              7*B_x*B_y*b_x*b_y^3)*mu_x3*mu_y2 - 8*(B_y^2*b_x*b_y^3*mu_x2 +
                                                                                                                                                                                                                                                                                                                                                                                      B_y^2*b_y^4*mu_y2)*mu_y3)*mu_y4 + 2*(B_y^2*b_x^2*b_y^2*mu_x2^2 +
                                                                                                                                                                                                                                                                                                                                                                                                                             2*B_y^2*b_x*b_y^3*mu_x2*mu_y2 + B_y^2*b_y^4*mu_y2^2)*mu_y5)*A +
         15*(20*B_x^3*b_x^3*b_y*mu_x2^3*mu_x3 - 21*B_x^3*b_x^3*b_y*mu_x3^3 -
               20*B_x^3*b_x^3*b_y*mu_x2*mu_x3*mu_x4)*mu_y2 -
         15*(10*B_x^2*B_y*b_x^3*b_y*mu_x2^4 + 20*B_y^3*b_x*b_y^3*mu_x2*mu_y2^3 +
               10*B_y^3*b_y^4*mu_y2^4 - 63*B_x^2*B_y*b_x^3*b_y*mu_x2*mu_x3^2 -
               10*B_x^2*B_y*b_x^3*b_y*mu_x2^2*mu_x4 - 10*(B_x^2*B_y*b_x*b_y^3*mu_x4 -
                                                            (B_y^3*b_x^2*b_y^2 + B_x^2*B_y*b_x*b_y^3)*mu_x2^2)*mu_y2^2 +
               (20*B_x^2*B_y*b_x^2*b_y^2*mu_x2^3 - 63*B_x^2*B_y*b_x^2*b_y^2*mu_x3^2 -
                  20*B_x^2*B_y*b_x^2*b_y^2*mu_x2*mu_x4)*mu_y2)*mu_y3 -
         150*(B_x*B_y^2*b_x^3*b_y*mu_x2^2*mu_x3 +
                2*B_x*B_y^2*b_x^2*b_y^2*mu_x2*mu_x3*mu_y2 +
                B_x*B_y^2*b_x*b_y^3*mu_x3*mu_y2^2 - (B_y^3*b_x^2*b_y^2*mu_x2^2 +
                                                       2*B_y^3*b_x*b_y^3*mu_x2*mu_y2 +
                                                       B_y^3*b_y^4*mu_y2^2)*mu_y3)*mu_y4)/A^(11/2)
}
#' @rdname kfuns2
K61two <- function(A, B_x, B_y, b_x, b_y, mu_x2, mu_x3, mu_x4, mu_x5, mu_x6,
                   mu_y2, mu_y3, mu_y4, mu_y5, mu_y6) {
  1/2*(180*B_x^4*b_x^5*mu_x2^7 + 540*B_y^4*b_x*b_y^4*mu_x2*mu_y2^6 +
         180*B_y^4*b_y^5*mu_y2^7 - 3285*B_x^4*b_x^5*mu_x2^4*mu_x3^2 +
         3600*B_x^4*b_x^5*mu_x2*mu_x3^4 + 180*B_x^4*b_x^5*mu_x2^3*mu_x4^2 -
         180*(2*B_x^2*B_y^2*b_x*b_y^4*mu_x4 - (3*B_y^4*b_x^2*b_y^3 +
                                                 2*B_x^2*B_y^2*b_x*b_y^4)*mu_x2^2)*mu_y2^5 + 2*(30*b_x^5*mu_x2^3 +
                                                                                                  30*b_y^5*mu_y2^3 - 10*b_x^5*mu_x3^2 - 15*b_x^5*mu_x2*mu_x4 -
                                                                                                  10*b_y^5*mu_y3^2 - 15*b_y^5*mu_y2*mu_y4 + b_x^5*mu_x6 + b_y^5*mu_y6)*A^4
       - 45*(73*B_x^2*B_y^2*b_x^2*b_y^3*mu_x3^2 +
               24*B_x^2*B_y^2*b_x^2*b_y^3*mu_x2*mu_x4 - 4*(B_y^4*b_x^3*b_y^2 +
                                                             6*B_x^2*B_y^2*b_x^2*b_y^3)*mu_x2^3)*mu_y2^4 +
         3600*(B_y^4*b_x*b_y^4*mu_x2 + B_y^4*b_y^5*mu_y2)*mu_y3^4 -
         30*(48*B_x*b_x^5*mu_x2^4 + 30*B_y*b_x*b_y^4*mu_x2*mu_y2^3 +
               48*B_y*b_y^5*mu_y2^4 - 40*B_x*b_x^5*mu_x2*mu_x3^2 -
               27*B_x*b_x^5*mu_x2^2*mu_x4 + 2*B_x*b_x^5*mu_x4^2 +
               3*B_x*b_x^5*mu_x3*mu_x5 + B_x*b_x^5*mu_x2*mu_x6 + 2*B_y*b_y^5*mu_y4^2 +
               6*(3*(B_y*b_x^3*b_y^2 + B_x*b_x^2*b_y^3)*mu_x2^2 - (B_y*b_x^3*b_y^2 +
                                                                     B_x*b_x^2*b_y^3)*mu_x4)*mu_y2^2 - 10*(B_y*b_x*b_y^4*mu_x2 +
                                                                                                             4*B_y*b_y^5*mu_y2)*mu_y3^2 + (30*B_x*b_x^4*b_y*mu_x2^3 -
                                                                                                                                             10*B_x*b_x^4*b_y*mu_x3^2 - 15*B_x*b_x^4*b_y*mu_x2*mu_x4 +
                                                                                                                                             B_x*b_x^4*b_y*mu_x6)*mu_y2 + (10*(B_y*b_x^4*b_y +
                                                                                                                                                                                 2*B_x*b_x^3*b_y^2)*mu_x2*mu_x3 + 10*(2*B_y*b_x^2*b_y^3 +
                                                                                                                                                                                                                        B_x*b_x*b_y^4)*mu_x3*mu_y2 - (B_y*b_x^4*b_y +
                                                                                                                                                                                                                                                        2*B_x*b_x^3*b_y^2)*mu_x5)*mu_y3 - (15*B_y*b_x*b_y^4*mu_x2*mu_y2 +
                                                                                                                                                                                                                                                                                             27*B_y*b_y^5*mu_y2^2 + 6*(B_y*b_x^3*b_y^2 + B_x*b_x^2*b_y^3)*mu_x2^2 -
                                                                                                                                                                                                                                                                                             2*(B_y*b_x^3*b_y^2 + B_x*b_x^2*b_y^3)*mu_x4)*mu_y4 + (3*B_y*b_y^5*mu_y3
                                                                                                                                                                                                                                                                                                                                                   - (2*B_y*b_x^2*b_y^3 + B_x*b_x*b_y^4)*mu_x3)*mu_y5 +
               (B_y*b_x*b_y^4*mu_x2 + B_y*b_y^5*mu_y2)*mu_y6)*A^3 -
         90*(73*B_x^2*B_y^2*b_x^3*b_y^2*mu_x2*mu_x3^2 -
               2*B_x^4*b_x^2*b_y^3*mu_x4^2 - 2*(6*B_x^2*B_y^2*b_x^3*b_y^2 +
                                                  B_x^4*b_x^2*b_y^3)*mu_x2^4 + 4*(3*B_x^2*B_y^2*b_x^3*b_y^2 +
                                                                                    B_x^4*b_x^2*b_y^3)*mu_x2^2*mu_x4)*mu_y2^3 -
         14400*(B_x*B_y^3*b_x^2*b_y^3*mu_x2*mu_x3 +
                  B_x*B_y^3*b_x*b_y^4*mu_x3*mu_y2)*mu_y3^3 + 15*(336*B_x^2*b_x^5*mu_x2^5 +
                                                                   444*B_y^2*b_x*b_y^4*mu_x2*mu_y2^4 + 336*B_y^2*b_y^5*mu_y2^5 -
                                                                   623*B_x^2*b_x^5*mu_x2^2*mu_x3^2 + 28*B_x^2*b_x^5*mu_x2*mu_x4^2 +
                                                                   44*B_x^2*b_x^5*mu_x2*mu_x3*mu_x5 + 6*B_x^2*b_x^5*mu_x2^2*mu_x6 +
                                                                   4*(3*(36*B_x*B_y*b_x^2*b_y^3 + B_x^2*b_x*b_y^4 + (b_x^3*b_y^2 +
                                                                                                                       9*b_x^2*b_y^3)*B_y^2)*mu_x2^2 - (B_y^2*b_x^3*b_y^2 +
                                                                                                                                                          36*B_x*B_y*b_x^2*b_y^3 + 3*B_x^2*b_x*b_y^4)*mu_x4)*mu_y2^3 +
                                                                   (6*B_x^2*b_x^3*b_y^2*mu_x6 + 12*(B_y^2*b_x^4*b_y +
                                                                                                      36*B_x*B_y*b_x^3*b_y^2 + (9*b_x^3*b_y^2 + b_x^2*b_y^3)*B_x^2)*mu_x2^3 -
                                                                      3*(B_y^2*b_x^4*b_y + 48*B_x*B_y*b_x^3*b_y^2 + 4*(3*b_x^3*b_y^2 +
                                                                                                                         4*b_x^2*b_y^3)*B_x^2)*mu_x3^2 - 2*(2*B_y^2*b_x^4*b_y +
                                                                                                                                                              72*B_x*B_y*b_x^3*b_y^2 + 3*(11*b_x^3*b_y^2 +
                                                                                                                                                                                            2*b_x^2*b_y^3)*B_x^2)*mu_x2*mu_x4)*mu_y2^2 -
                                                                   (464*B_y^2*b_x*b_y^4*mu_x2*mu_y2 + 623*B_y^2*b_y^5*mu_y2^2 +
                                                                      3*(48*B_x*B_y*b_x^2*b_y^3 + B_x^2*b_x*b_y^4 + 4*(4*b_x^3*b_y^2 +
                                                                                                                         3*b_x^2*b_y^3)*B_y^2)*mu_x2^2 - (16*B_y^2*b_x^3*b_y^2 +
                                                                                                                                                            48*B_x*B_y*b_x^2*b_y^3 + 3*B_x^2*b_x*b_y^4)*mu_x4)*mu_y3^2 +
                                                                   28*(B_y^2*b_x*b_y^4*mu_x2 + B_y^2*b_y^5*mu_y2)*mu_y4^2 -
                                                                   (226*B_x^2*b_x^5*mu_x2^3 - 67*B_x^2*b_x^5*mu_x3^2)*mu_x4 +
                                                                   4*(111*B_x^2*b_x^4*b_y*mu_x2^4 - 116*B_x^2*b_x^4*b_y*mu_x2*mu_x3^2 -
                                                                        73*B_x^2*b_x^4*b_y*mu_x2^2*mu_x4 + 7*B_x^2*b_x^4*b_y*mu_x4^2 +
                                                                        11*B_x^2*b_x^4*b_y*mu_x3*mu_x5 + 3*B_x^2*b_x^4*b_y*mu_x2*mu_x6)*mu_y2 +
                                                                   2*((280*B_x*B_y*b_x^4*b_y + 111*B_x^2*b_x^3*b_y^2)*mu_x2^2*mu_x3 +
                                                                        (111*B_y^2*b_x^2*b_y^3 + 280*B_x*B_y*b_x*b_y^4)*mu_x3*mu_y2^2 -
                                                                        (40*B_x*B_y*b_x^4*b_y + 27*B_x^2*b_x^3*b_y^2)*mu_x3*mu_x4 -
                                                                        2*(8*B_x*B_y*b_x^4*b_y + 3*B_x^2*b_x^3*b_y^2)*mu_x2*mu_x5 +
                                                                        2*(2*(9*B_y^2*b_x^3*b_y^2 + 9*B_x^2*b_x^2*b_y^3 + 40*(b_x^3*b_y^2 +
                                                                                                                                b_x^2*b_y^3)*B_x*B_y)*mu_x2*mu_x3 - (8*B_x*B_y*b_x^3*b_y^2 +
                                                                                                                                                                       3*B_x^2*b_x^2*b_y^3)*mu_x5)*mu_y2)*mu_y3 -
                                                                   (292*B_y^2*b_x*b_y^4*mu_x2*mu_y2^2 + 226*B_y^2*b_y^5*mu_y2^3 -
                                                                      67*B_y^2*b_y^5*mu_y3^2 + 4*(3*B_y^2*b_x^4*b_y + 36*B_x*B_y*b_x^3*b_y^2 +
                                                                                                    B_x^2*b_x^2*b_y^3)*mu_x2^3 - (3*B_y^2*b_x^4*b_y + 48*B_x*B_y*b_x^3*b_y^2
                                                                                                                                  + 16*B_x^2*b_x^2*b_y^3)*mu_x3^2 - 4*(B_y^2*b_x^4*b_y +
                                                                                                                                                                         12*B_x*B_y*b_x^3*b_y^2 + B_x^2*b_x^2*b_y^3)*mu_x2*mu_x4 +
                                                                      2*(27*B_y^2*b_x^2*b_y^3 + 40*B_x*B_y*b_x*b_y^4)*mu_x3*mu_y3 +
                                                                      2*((72*B_x*B_y*b_x^2*b_y^3 + 2*B_x^2*b_x*b_y^4 + 3*(2*b_x^3*b_y^2 +
                                                                                                                            11*b_x^2*b_y^3)*B_y^2)*mu_x2^2 - 2*(B_y^2*b_x^3*b_y^2 +
                                                                                                                                                                  12*B_x*B_y*b_x^2*b_y^3 + B_x^2*b_x*b_y^4)*mu_x4)*mu_y2)*mu_y4 -
                                                                   4*((3*B_y^2*b_x^3*b_y^2 + 8*B_x*B_y*b_x^2*b_y^3)*mu_x2*mu_x3 +
                                                                        (3*B_y^2*b_x^2*b_y^3 + 8*B_x*B_y*b_x*b_y^4)*mu_x3*mu_y2 -
                                                                        11*(B_y^2*b_x*b_y^4*mu_x2 + B_y^2*b_y^5*mu_y2)*mu_y3)*mu_y5 +
                                                                   6*(B_y^2*b_x^2*b_y^3*mu_x2^2 + 2*B_y^2*b_x*b_y^4*mu_x2*mu_y2 +
                                                                        B_y^2*b_y^5*mu_y2^2)*mu_y6)*A^2 + 45*(12*B_x^4*b_x^3*b_y^2*mu_x2*mu_x4^2
                                                                                                              + 4*(2*B_x^2*B_y^2*b_x^4*b_y + 3*B_x^4*b_x^3*b_y^2)*mu_x2^5 -
                                                                                                                73*(B_x^2*B_y^2*b_x^4*b_y + B_x^4*b_x^3*b_y^2)*mu_x2^2*mu_x3^2 +
                                                                                                                (73*B_x^4*b_x^3*b_y^2*mu_x3^2 - 8*(B_x^2*B_y^2*b_x^4*b_y +
                                                                                                                                                     3*B_x^4*b_x^3*b_y^2)*mu_x2^3)*mu_x4)*mu_y2^2 -
         45*(73*B_x^2*B_y^2*b_x^3*b_y^2*mu_x2^4 +
               146*B_y^4*b_x*b_y^4*mu_x2*mu_y2^3 + 73*B_y^4*b_y^5*mu_y2^4 -
               480*B_x^2*B_y^2*b_x^3*b_y^2*mu_x2*mu_x3^2 -
               73*B_x^2*B_y^2*b_x^3*b_y^2*mu_x2^2*mu_x4 -
               73*(B_x^2*B_y^2*b_x*b_y^4*mu_x4 - (B_y^4*b_x^2*b_y^3 +
                                                    B_x^2*B_y^2*b_x*b_y^4)*mu_x2^2)*mu_y2^2 +
               2*(73*B_x^2*B_y^2*b_x^2*b_y^3*mu_x2^3 -
                    240*B_x^2*B_y^2*b_x^2*b_y^3*mu_x3^2 -
                    73*B_x^2*B_y^2*b_x^2*b_y^3*mu_x2*mu_x4)*mu_y2)*mu_y3^2 +
         180*(B_y^4*b_x^3*b_y^2*mu_x2^3 + 3*B_y^4*b_x^2*b_y^3*mu_x2^2*mu_y2 +
                3*B_y^4*b_x*b_y^4*mu_x2*mu_y2^2 + B_y^4*b_y^5*mu_y2^3)*mu_y4^2 -
         30*(56*B_x^3*b_x^5*mu_x2^6 + 114*B_y^3*b_x*b_y^4*mu_x2*mu_y2^5 +
               56*B_y^3*b_y^5*mu_y2^6 - 597*B_x^3*b_x^5*mu_x2^3*mu_x3^2 +
               40*B_x^3*b_x^5*mu_x3^4 + 18*B_x^3*b_x^5*mu_x2^2*mu_x4^2 +
               33*B_x^3*b_x^5*mu_x2^2*mu_x3*mu_x5 + B_x^3*b_x^5*mu_x2^3*mu_x6 +
               40*B_y^3*b_y^5*mu_y3^4 + 6*((9*B_x*B_y^2*b_x^2*b_y^3 +
                                              10*B_y^3*b_x^2*b_y^3 + 9*B_x^2*B_y*b_x*b_y^4)*mu_x2^2 -
                                             3*(B_x*B_y^2*b_x^2*b_y^3 + 3*B_x^2*B_y*b_x*b_y^4)*mu_x4)*mu_y2^4 -
               40*(B_y^3*b_x^2*b_y^3 + 3*B_x*B_y^2*b_x*b_y^4)*mu_x3*mu_y3^3 +
               (B_x^3*b_x^2*b_y^3*mu_x6 + 2*(54*B_x*B_y^2*b_x^3*b_y^2 +
                                               B_y^3*b_x^3*b_y^2 + B_x^3*b_x^2*b_y^3 +
                                               54*B_x^2*B_y*b_x^2*b_y^3)*mu_x2^3 - 3*(11*B_x*B_y^2*b_x^3*b_y^2 +
                                                                                        2*B_x^3*b_x^2*b_y^3 + 120*B_x^2*B_y*b_x^2*b_y^3)*mu_x3^2 -
                  3*(12*B_x*B_y^2*b_x^3*b_y^2 + B_x^3*b_x^2*b_y^3 +
                       36*B_x^2*B_y*b_x^2*b_y^3)*mu_x2*mu_x4)*mu_y2^3 +
               3*(6*B_x^3*b_x^3*b_y^2*mu_x4^2 + 11*B_x^3*b_x^3*b_y^2*mu_x3*mu_x5 +
                    B_x^3*b_x^3*b_y^2*mu_x2*mu_x6 + 2*(9*B_x*B_y^2*b_x^4*b_y +
                                                         10*B_x^3*b_x^3*b_y^2 + 9*B_x^2*B_y*b_x^3*b_y^2)*mu_x2^4 -
                    (11*B_x*B_y^2*b_x^4*b_y + 72*B_x^3*b_x^3*b_y^2 +
                       120*B_x^2*B_y*b_x^3*b_y^2)*mu_x2*mu_x3^2 - 3*(2*B_x*B_y^2*b_x^4*b_y +
                                                                       9*B_x^3*b_x^3*b_y^2 + 6*B_x^2*B_y*b_x^3*b_y^2)*mu_x2^2*mu_x4)*mu_y2^2 -
               3*(269*B_y^3*b_x*b_y^4*mu_x2*mu_y2^2 + 199*B_y^3*b_y^5*mu_y2^3 +
                    (120*B_x*B_y^2*b_x^3*b_y^2 + 2*B_y^3*b_x^3*b_y^2 +
                       11*B_x^2*B_y*b_x^2*b_y^3)*mu_x2^3 - 40*(B_x*B_y^2*b_x^3*b_y^2 +
                                                                 B_x^2*B_y*b_x^2*b_y^3)*mu_x3^2 - (40*B_x*B_y^2*b_x^3*b_y^2 +
                                                                                                     11*B_x^2*B_y*b_x^2*b_y^3)*mu_x2*mu_x4 + ((120*B_x*B_y^2*b_x^2*b_y^3 +
                                                                                                                                                 72*B_y^3*b_x^2*b_y^3 + 11*B_x^2*B_y*b_x*b_y^4)*mu_x2^2 -
                                                                                                                                                (40*B_x*B_y^2*b_x^2*b_y^3 +
                                                                                                                                                   11*B_x^2*B_y*b_x*b_y^4)*mu_x4)*mu_y2)*mu_y3^2 +
               18*(B_y^3*b_x^2*b_y^3*mu_x2^2 + 2*B_y^3*b_x*b_y^4*mu_x2*mu_y2 +
                     B_y^3*b_y^5*mu_y2^2)*mu_y4^2 - 3*(25*B_x^3*b_x^5*mu_x2^4 -
                                                         51*B_x^3*b_x^5*mu_x2*mu_x3^2)*mu_x4 + 3*(38*B_x^3*b_x^4*b_y*mu_x2^5 -
                                                                                                    269*B_x^3*b_x^4*b_y*mu_x2^2*mu_x3^2 + 12*B_x^3*b_x^4*b_y*mu_x2*mu_x4^2 +
                                                                                                    22*B_x^3*b_x^4*b_y*mu_x2*mu_x3*mu_x5 + B_x^3*b_x^4*b_y*mu_x2^2*mu_x6 -
                                                                                                    51*(B_x^3*b_x^4*b_y*mu_x2^3 - B_x^3*b_x^4*b_y*mu_x3^2)*mu_x4)*mu_y2 -
               (33*B_x^2*B_y*b_x^4*b_y*mu_x2^2*mu_x5 - 3*(317*B_x^2*B_y*b_x^4*b_y +
                                                            11*B_x^3*b_x^3*b_y^2)*mu_x2^3*mu_x3 - 3*(11*B_y^3*b_x^2*b_y^3 +
                                                                                                       317*B_x*B_y^2*b_x*b_y^4)*mu_x3*mu_y2^3 + 40*(3*B_x^2*B_y*b_x^4*b_y +
                                                                                                                                                      B_x^3*b_x^3*b_y^2)*mu_x3^3 + 3*(91*B_x^2*B_y*b_x^4*b_y +
                                                                                                                                                                                        11*B_x^3*b_x^3*b_y^2)*mu_x2*mu_x3*mu_x4 +
                  3*(11*B_x^2*B_y*b_x^2*b_y^3*mu_x5 - (11*B_y^3*b_x^3*b_y^2 +
                                                         66*B_x^2*B_y*b_x^2*b_y^3 +
                                                         383*B_x*B_y^2*b_x^2*b_y^3)*mu_x2*mu_x3)*mu_y2^2 +
                  3*(22*B_x^2*B_y*b_x^3*b_y^2*mu_x2*mu_x5 - (383*B_x^2*B_y*b_x^3*b_y^2 +
                                                               66*B_x*B_y^2*b_x^3*b_y^2 + 11*B_x^3*b_x^2*b_y^3)*mu_x2^2*mu_x3 +
                       (91*B_x^2*B_y*b_x^3*b_y^2 +
                          11*B_x^3*b_x^2*b_y^3)*mu_x3*mu_x4)*mu_y2)*mu_y3 -
               3*(51*B_y^3*b_x*b_y^4*mu_x2*mu_y2^3 + 25*B_y^3*b_y^5*mu_y2^4 +
                    6*(3*B_x*B_y^2*b_x^4*b_y + B_x^2*B_y*b_x^3*b_y^2)*mu_x2^4 -
                    (11*B_x*B_y^2*b_x^4*b_y + 40*B_x^2*B_y*b_x^3*b_y^2)*mu_x2*mu_x3^2 -
                    6*(B_x*B_y^2*b_x^4*b_y + B_x^2*B_y*b_x^3*b_y^2)*mu_x2^2*mu_x4 +
                    3*((6*B_x*B_y^2*b_x^2*b_y^3 + 9*B_y^3*b_x^2*b_y^3 +
                          2*B_x^2*B_y*b_x*b_y^4)*mu_x2^2 - 2*(B_x*B_y^2*b_x^2*b_y^3 +
                                                                B_x^2*B_y*b_x*b_y^4)*mu_x4)*mu_y2^2 - 51*(B_y^3*b_x*b_y^4*mu_x2 +
                                                                                                            B_y^3*b_y^5*mu_y2)*mu_y3^2 + ((36*B_x*B_y^2*b_x^3*b_y^2 +
                                                                                                                                             B_y^3*b_x^3*b_y^2 + 12*B_x^2*B_y*b_x^2*b_y^3)*mu_x2^3 -
                                                                                                                                            (11*B_x*B_y^2*b_x^3*b_y^2 + 40*B_x^2*B_y*b_x^2*b_y^3)*mu_x3^2 -
                                                                                                                                            12*(B_x*B_y^2*b_x^3*b_y^2 + B_x^2*B_y*b_x^2*b_y^3)*mu_x2*mu_x4)*mu_y2 +
                    ((11*B_y^3*b_x^3*b_y^2 + 91*B_x*B_y^2*b_x^2*b_y^3)*mu_x2*mu_x3 +
                       (11*B_y^3*b_x^2*b_y^3 +
                          91*B_x*B_y^2*b_x*b_y^4)*mu_x3*mu_y2)*mu_y3)*mu_y4 -
               33*(B_x*B_y^2*b_x^3*b_y^2*mu_x2^2*mu_x3 +
                     2*B_x*B_y^2*b_x^2*b_y^3*mu_x2*mu_x3*mu_y2 +
                     B_x*B_y^2*b_x*b_y^4*mu_x3*mu_y2^2 - (B_y^3*b_x^2*b_y^3*mu_x2^2 +
                                                            2*B_y^3*b_x*b_y^4*mu_x2*mu_y2 + B_y^3*b_y^5*mu_y2^2)*mu_y3)*mu_y5 +
               (B_y^3*b_x^3*b_y^2*mu_x2^3 + 3*B_y^3*b_x^2*b_y^3*mu_x2^2*mu_y2 +
                  3*B_y^3*b_x*b_y^4*mu_x2*mu_y2^2 + B_y^3*b_y^5*mu_y2^3)*mu_y6)*A -
         45*(8*B_x^4*b_x^5*mu_x2^5 - 73*B_x^4*b_x^5*mu_x2^2*mu_x3^2)*mu_x4 +
         90*(6*B_x^4*b_x^4*b_y*mu_x2^6 - 73*B_x^4*b_x^4*b_y*mu_x2^3*mu_x3^2 +
               40*B_x^4*b_x^4*b_y*mu_x3^4 + 6*B_x^4*b_x^4*b_y*mu_x2^2*mu_x4^2 -
               (12*B_x^4*b_x^4*b_y*mu_x2^4 -
                  73*B_x^4*b_x^4*b_y*mu_x2*mu_x3^2)*mu_x4)*mu_y2 +
         90*(73*B_x^3*B_y*b_x^4*b_y*mu_x2^4*mu_x3 +
               146*B_x*B_y^3*b_x^2*b_y^3*mu_x2*mu_x3*mu_y2^3 +
               73*B_x*B_y^3*b_x*b_y^4*mu_x3*mu_y2^4 -
               160*B_x^3*B_y*b_x^4*b_y*mu_x2*mu_x3^3 -
               73*B_x^3*B_y*b_x^4*b_y*mu_x2^2*mu_x3*mu_x4 -
               73*(B_x^3*B_y*b_x^2*b_y^3*mu_x3*mu_x4 - (B_x*B_y^3*b_x^3*b_y^2 +
                                                          B_x^3*B_y*b_x^2*b_y^3)*mu_x2^2*mu_x3)*mu_y2^2 +
               2*(73*B_x^3*B_y*b_x^3*b_y^2*mu_x2^3*mu_x3 -
                    80*B_x^3*B_y*b_x^3*b_y^2*mu_x3^3 -
                    73*B_x^3*B_y*b_x^3*b_y^2*mu_x2*mu_x3*mu_x4)*mu_y2)*mu_y3 -
         45*(8*B_x^2*B_y^2*b_x^4*b_y*mu_x2^5 + 24*B_y^4*b_x*b_y^4*mu_x2*mu_y2^4 +
               8*B_y^4*b_y^5*mu_y2^5 - 73*B_x^2*B_y^2*b_x^4*b_y*mu_x2^2*mu_x3^2 -
               8*B_x^2*B_y^2*b_x^4*b_y*mu_x2^3*mu_x4 - 8*(B_x^2*B_y^2*b_x*b_y^4*mu_x4 -
                                                            (3*B_y^4*b_x^2*b_y^3 + B_x^2*B_y^2*b_x*b_y^4)*mu_x2^2)*mu_y2^3 -
               (73*B_x^2*B_y^2*b_x^2*b_y^3*mu_x3^2 +
                  24*B_x^2*B_y^2*b_x^2*b_y^3*mu_x2*mu_x4 - 8*(B_y^4*b_x^3*b_y^2 +
                                                                3*B_x^2*B_y^2*b_x^2*b_y^3)*mu_x2^3)*mu_y2^2 -
               73*(B_y^4*b_x^2*b_y^3*mu_x2^2 + 2*B_y^4*b_x*b_y^4*mu_x2*mu_y2 +
                     B_y^4*b_y^5*mu_y2^2)*mu_y3^2 + 2*(12*B_x^2*B_y^2*b_x^3*b_y^2*mu_x2^4 -
                                                         73*B_x^2*B_y^2*b_x^3*b_y^2*mu_x2*mu_x3^2 -
                                                         12*B_x^2*B_y^2*b_x^3*b_y^2*mu_x2^2*mu_x4)*mu_y2 +
               146*(B_x*B_y^3*b_x^3*b_y^2*mu_x2^2*mu_x3 +
                      2*B_x*B_y^3*b_x^2*b_y^3*mu_x2*mu_x3*mu_y2 +
                      B_x*B_y^3*b_x*b_y^4*mu_x3*mu_y2^2)*mu_y3)*mu_y4)/A^7
}	


