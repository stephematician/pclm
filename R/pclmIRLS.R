#' Iteratively re-weighted least squares for Penalised Composite Link Models
#'
#' Applies IRLS to a find the maximum value of a penalised likelihood expression
#' representing the likelihood of a composite link model.
#'
#' Silvia Rizzi, Jutta Gampe, and Paul H. C. Eilers, 2015, Efficient Estimation
#' of Smooth Distributions From Coarsely Grouped Data Am. J. Epidemiol. first
#' published online June 16, 2015 doi.10.1093/aje/kwv020
#'
#' @param n a numeric vector containing grouped count data
#' @param C a matrix for the mapping of bins to groups in the model
#' @param X a matrix to link a parameterisation of the probability
#'    distribution for the binned data
#' @param lambda a scalar smoothness penalty
#' @param deg the degree of the smoothness operator
#'
#' @return a list with the following names components
#' \itemize{
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
pclmIRLS <- function(n, C, X=diag(nrow=ncol(C)), lambda=1, deg=2) {

    max_iter <- 100
    tolerance <- 1e-6

    n_bin <- ncol(X)
    D_penalty <- diff(diag(n_bin), diff=deg)

    # Initial guess
    b <- matrix(log(sum(n) / n_bin), n_bin, 1)

    if (lambda <= 0)
        warning(paste('pclmfit:pclmIRLS() lambda > 0 required, setting',
                      'sqrt(lambda) to eps'))

    lambda <- max(lambda, .Machine$double.eps^2)

    # Loop invariants
    b_penalty  <- matrix(0, nrow(D_penalty), 1)
    wt_penalty <- matrix(lambda, nrow(D_penalty), 1)

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
        Q <- CX %*% diag(as.vector(gamma_hat))
        z <- (n - mu_hat) + Q %*% b_hat
        wt_CLM <- 1 / mu_hat

        b <- lsfit(rbind(Q, D_penalty),
                   c(z, b_penalty),
                   wt=c(wt_CLM, wt_penalty),
                   intercept=F)$coef

        delta_b <- max(abs(b - b_hat));

        if (delta_b < tolerance)
            break

        if (j == max_iter)
            converged <- F
    }

    gamma <- exp(X %*% b)
    mu <- C %*% gamma
    plot(gamma)
    Q <- CX %*% diag(as.vector(gamma))
    R <- t(Q) %*% diag(as.vector(1 / mu)) %*% Q
    H <- solve(R + lambda * t(D_penalty) %*% D_penalty) %*% R
    tr <- sum(diag(H))
    nzv <- n > 0 & mu > 0
    dev <- 2 * sum(n[nzv] * log(n[nzv] / mu[nzv]))

    fit <- list(mu=mu,   gamma=gamma,      trace=tr,
                dev=dev, aic=dev + 2 * tr, converged=converged)

}

