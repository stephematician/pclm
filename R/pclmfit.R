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
#' @param n numeric (vector); counts in each group
#' @param x numeric (vector); the index of the group for each bin (in sequence)
#' @param w numeric (vector); (within group) weights for each bin
#' @param lambda0 numeric (scalar); the (initial) value of penalty (default
#      \code{1})
#' @param do_optim boolean; whether to perform optimisation of AIC or not
#'     (default \code{T})
#' @param use_spline boolean; whether or not to use a spline from
#'      \code{MortSmooth_bbase} for the basis of the expected count in each 
#'      bin (default \code{F})
#' @param k integer (scalar); number of inner knots to use in the spline
#'     describing the expected count in each bin (default \code{4})
#' @param search_range (scalar): a factor for determining the range of the
#'     search (for optimal AIC) in log penalty space, default is `1e6` i.e. to
#'     search from `lambda0 / 1e3` to `1e3 * lambda0`
#' @param ... further arguments passed on to \code{optim} and \code{pclmIRLS}
#'
#' @return a list with the following named components
#' \describe{
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
#' @importFrom stats optimize
#' @importFrom utils modifyList
#' @importFrom splines ns
#' @export
pclmfit <- function(n, x, w=rep(1, length(x)), lambda0=1,
                    do_optim=T, use_spline=F, k=4, search_range=1e4, ...) {

    stopifnot(is.vector(n), is.numeric(n), all(n >= 0))
    stopifnot(is.vector(x), is.numeric(x), all(x %in% 1:length(n)))
    stopifnot(identical(x, sort(x)))
    stopifnot(all(1:length(n) %in% x))
    stopifnot(is.vector(w), is.numeric(w), all(w > 0))

    dx <- 1E-5 # TODO: let user modify

    C <- matrix(0, nrow=length(n), ncol=length(x))
    C[cbind(x,1:length(x))] <- w

    pclm_args <- list(n, C)

    if (use_spline) {
        # Use natural cubic-spline basis with k inner knots
        pclm_args$X <- ns(1:length(x), df=k + 2L, intercept=T)
        pclm_args$D <- (dx^-2) * with(pclm_args,
                                      predict(X, newx=2:(length(x)-1) + dx) +
                                          predict(X,
                                                  newx=2:(length(x)-1) - dx) -
                                          2 * X[-c(1, nrow(X)),])
    } else
        # Identity matrix
        pclm_args$X <- diag(nrow=ncol(C))

    pclm_search <- function(ln_lambda) {
        fit <- do.call(pclmIRLS,
                       modifyList(pclm_args,
                                  list(lambda=exp(ln_lambda))))
        if (!fit$converged)
            Inf
        else
            fit$aic
    }

    ln_lambda_fit <- log(lambda0)

    if (do_optim) {

        optim_args <- c(list(pclm_search),
                        modifyList(
                            list(lower=ln_lambda_fit - 0.5 * log(search_range),
                                 upper=ln_lambda_fit + 0.5 * log(search_range)),
                            list(...)
                        ))

        repeat {

            # search for optimal solution in current range
            ln_lambda_fit <- do.call(optimize, optim_args)$minimum

            # if we've searched as far as we realisticall can; break
            if (optim_args$lower < 0.5 * log(.Machine$double.eps)) break
            if (optim_args$upper > 0.5 * log(.Machine$double.xmax)) break

            # minimum distance between two points according to `optimize`
            min_d <- .Machine$double.eps^0.25 * abs(ln_lambda_fit) +
                         min(optim_args$tol,
                             .Machine$double.eps^0.25,
                             na.rm=T) / 3

            if (ln_lambda_fit - optim_args$lower < 1e1 * min_d) {
                # too close to LHS, expand search
                optim_args$lower <- ln_lambda_fit - log(search_range)
                optim_args$upper <- ln_lambda_fit + 1e1 * min_d
            } else if (optim_args$upper - ln_lambda_fit < 1e1 * min_d) {
                # too close to RHS, expand search 
                optim_args$lower <- ln_lambda_fit - 1e1 * min_d
                optim_args$upper <- ln_lambda_fit + log(search_range)
            } else
                # we are sufficiently far from boundary, end search
                break

        }

    }

    fit <- c(do.call(pclmIRLS,
                     modifyList(pclm_args,
                                list(lambda=exp(ln_lambda_fit)))),
             list(lambda=exp(ln_lambda_fit),
                  bin_group=x))

}

