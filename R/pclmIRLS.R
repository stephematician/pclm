#' Iteratively re-weighted least squares for Penalised Composite Link Models
#'
#' Applies IRLS to a find the maximum value of a penalised likelihood expression
#' representing the likelihood of a composite link model.
#'
#' Solves the system:
#'
#' \deqn{X^TWXb = X^TW[inv(W)(y - \mu) + X^Tb}
#' \deqn{\sqrt{lambda} D = 0}
#'
#' Where \eqn{X} is the transformation between the rate in each bin and the
#' subspace that the rates are being modelled on, by default the identity
#' matrix is used, but a spline basis can be used if there are too many
#' rates (bins) to estimate.
#'
#' \eqn{W} is constructed in each iteration so that the solution to the system
#' converges to the maximum likelihood estimate of \eqn{b}, using the usual
#' iteratively reweigthed least-squares for \eqn{L^1} norm minimisation.
#'
#' \eqn{D} is, by default, the second order difference matrix, a kind of
#' smoothing penalty applied to the system, and the value of lambda is a
#' smoothing parameter often determined by minimisation of the AIC.
#'
#' For further details on the use of splines and the penalty applied to the
#' coefficients, see
#'
#' Paul H. C. Eilers and Brian D. Marx, 1996, Flexible smoothing with B-splines
#' and penalties, \emph{Statistical Science, 11}(2), pp. 89-102,
#' doi:10.1214/ss/1038425655
#'
#' And for the original code upon which this was based, see (and erratum
#' thereof):
#'
#' Silvia Rizzi, Jutta Gampe, and Paul H. C. Eilers, 2015, Efficient estimation
#' of smooth distributions from coarsely grouped data, \emph{Am. J. Epidemiol.,
#' 182}(2), pp. 138-147, doi.10.1093/aje/kwv020
#'
#' @param n a numeric vector containing grouped count data
#' @param C a matrix for the mapping of bins to groups in the model
#' @param X a matrix to link a parameterisation of the probability
#'    distribution for the binned data
#' @param lambda a scalar smoothness penalty
#' @param deg the degree of the smoothness operator
#' @param D a matrix for penalising the likelihood (overrides \code{deg})
#'
#' @return a list with the following names components
#' \describe{
#'     \item{'mu'}{The (fitted) expected count in each group}
#'     \item{'gamma'}{The (fitted) expected count in each bin}
#'     \item{'trace'}{The effective dimension of the (fitted) model}
#'     \item{'dev'}{The deviance of gamma given mu}
#'     \item{'aic'}{The Aikake Information Criterion of the fitted model}
#'     \item{'converged'}{Convergence status of IRLS procedure}
#' }
#'
#' @importFrom stats lsfit
#' @export
pclmIRLS <- function(n, C, X=diag(nrow=ncol(C)), lambda=1, deg=2, D) {

    max_iter <- 100
    tolerance <- 1e-6

    n_bin <- ncol(X)
    if (missing(D))
        D_penalty <- diff(diag(n_bin), differences=deg)
    else
        D_penalty <- D

    # Initial guess - adapted to handle splines a bit better
    b <- matrix(lsfit(X,
                      t(C > 0) %*% pmax(log(n / rowSums(C)), -15), # magic
                      intercept=F)[['coefficients']],
                nrow=n_bin)

    if (lambda <= 0)
        warning(paste('pclmfit:pclmIRLS() lambda > 0 required, setting',
                      'lambda to sqrt(eps)'))

    lambda <- max(lambda, sqrt(.Machine[['double.eps']]))

    # Loop invariants
    b_penalty  <- matrix(0, nrow(D_penalty), 1)
    wt_penalty <- rep(lambda, nrow(D_penalty))

    if (nrow(X) == ncol(X))
        CX <- C %*% X

    converged <- T

    for (j in seq_len(max_iter)) {

        b_hat <- b
        # N * 1
        gamma_hat <- exp(X %*% b_hat)
        # M * 1
        mu_hat <- C %*% gamma_hat
        # R code in paper solves as the least squares fit of two systems:
        # CLM:     X'WXB = X'W[inv(W)(y - mu) + X'B]
        # penalty : sqrt(lambda) D = 0
        if (nrow(X) == ncol(X)) {
            Q <- CX %*% diag(as.vector(gamma_hat))
        } else
            Q <- C %*% diag(as.vector(gamma_hat)) %*% X
        z <- (n - mu_hat) + Q %*% b_hat
        wt_CLM <- 1 / mu_hat

        b <- lsfit(rbind(Q, D_penalty),
                   c(z, b_penalty),
                   wt=c(wt_CLM, wt_penalty),
                   intercept=F)[['coefficients']]

        delta_b <- max(abs(b - b_hat));

        if (delta_b < tolerance)
            break

        if (j == max_iter)
            converged <- F
    }

    gamma <- exp(X %*% b)
    mu <- C %*% gamma
    if (nrow(X) == ncol(X)) {
        Q <- CX %*% diag(as.vector(gamma)) 
    } else
        Q <- C %*% diag(as.vector(gamma)) %*% X
    R <- t(Q) %*% diag(as.vector(1 / mu)) %*% Q
    H <- solve(R + lambda * t(D_penalty) %*% D_penalty) %*% R
    tr <- sum(diag(H))
    nzv <- n > 0 & mu > 0
    dev <- 2 * sum(n[nzv] * log(n[nzv] / mu[nzv]))

    fit <- list(mu=mu,   gamma=gamma,      trace=tr,
                dev=dev, aic=dev + 2 * tr, converged=converged)

}

