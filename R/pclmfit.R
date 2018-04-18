#' Penalised Composite Link Model by minimisation of AIC
#'
#' Fits a penalised composite link model (PCLM) to grouped count data. This
#' provides an estimate of the distribution of objects in a specified number
#' of bins within each group (typically equi-spaced).
#'
#' The PCLM is given by minimising the Akaike Information Criterion (to specify
#' the smoothing parameter) as outlined in the paper:
#'
#' Silvia Rizzi, Jutta Gampe, and Paul H. C. Eilers, 2015, Efficient estimation
#' of smooth distributions from coarsely grouped data, \emph{Am. J. Epidemiol.,
#' 182}(2), pp. 138-147, doi.10.1093/aje/kwv020
#'
#' Allows for a within-group weighting in each bin, so that additional
#' (relative) information about the rates within a group can be incorporated.
#' For example, given two bins within a group, it might be known that one bin
#' is twice as likely as the other (aside: this is not the same as knowing that
#' one bin is twice as 'wide' as the other).
#'
#' If there are a large number of bins, it can be useful to model the rates (in
#' each bin) using splines. This is possible by specifying \code{use_spline=T}
#' and selecting a number of knots using the \code{k} argument (defaults to
#' four). The knot placement is given by the default behaviour of
#' \code{\link[splines]{bs}}. For further discussion see:
#'
#' Paul H. C. Eilers and Brian D. Marx, 1996, Flexible smoothing with B-splines
#' and penalties, \emph{Statistical Science, 11}(2), pp. 89-102.
#'
#' @param n a numeric vector of counts in each group
#' @param x a numeric vector of the index into n for each bin
#' @param w a numeric vector of (within group) weights for each bin
#' @param use_spline a boolean indicating whether or not to use a spline for the
#'     basis of the expected count in each bin (default false)
#' @param k the number of knots to use in the spline describing the expected
#'     count in each bin (default 4)
#'
#' @return a list with the following named components
#' \itemize{
#'     \item{'mu'}{The (fitted) expected count in each group}
#'     \item{'gamma'}{The (fitted) expected count in each bin}
#'     \item{'trace'}{The effective dimension of the (fitted) model}
#'     \item{'dev'}{The deviance of gamma given mu}
#'     \item{'aic'}{The Aikake Information Criterion of the fitted model}
#'     \item{'converged'}{Convergence status of IRLS procedure}
#'     \item{'lambda'}{Value of smoothing parameter that minimised AIC}
#'     \item{'bin_group'}{The group to which each bin belongs}
#' }
#'
#' @importFrom stats optim
#' @importFrom splines ns
#' @export
pclmfit <- function(n, x, w=rep(1, length(x)), use_spline=F, k=4) {

    stopifnot(is.vector(n), is.numeric(n), all(n >= 0))
    stopifnot(is.vector(x), is.numeric(x), all(x %in% 1:length(n)))
    stopifnot(identical(x, sort(x)))
    stopifnot(all(1:length(n) %in% x))
    stopifnot(is.vector(w), is.numeric(w), all(w > 0))

    C <- matrix(0, nrow=length(n), ncol=length(x))
    C[cbind(x,1:length(x))] <- w

    if (use_spline)
        # Use natural cubic-spline basis with k knots
        X <- ns(1:length(x), knots=k, intercept=T)
    else
        # Identity matrix
        X <- diag(nrow=ncol(C))

    pclm_search <- function(lambda) {
        fit <- pclmIRLS(n, C, X, lambda=lambda)
        if (!fit$converged)
            Inf
        else
            fit$aic
    }

    lambda0 <- optim(par=1,
                     pclm_search,
                     method='L-BFGS-B',
                     lower=.Machine$double.eps^2)$par

    fit <- c(pclmIRLS(n, C, X, lambda=lambda0),
             list(lambda=lambda0,
                  bin_group=x))

}

