#' Randomly draw samples from fitted PCLM
#'
#' Scales reasonably well (linear?) with
#' \itemize{
#'   \item number of observations
#'   \item number of unique (observed) intervals
#'   \item number of samples per interval
#' }
#'
#' Takes about 1 second on a 2016 desktop, Windows 7, for
#' \itemize{
#'   \item 100000 observations
#'   \item 2500 unique intervals
#'   \item 5 samples per interval
#' }
#'
#' Cannot handle PCLMs with more than \code{sqrt(.Machine$max.integer)} number
#' of bins due to a hashing mechanism used
#'
#' @param D A two-column data.frame of integers, each row represents an
#'          observation of an interval of bins, columns correspond to lower and
#'          upper bins of the interval (inclusive)
#' @param pclm_fit A list as returned by pclmfit
#' @param n_sample Number of samples to draw for each observation
#' @return Matrix of drawn bins
rpclmfit <- function(D, pclm_fit, n_sample=1L) {
    
    n_bin <- length(pclm_fit$gamma)

    # Check hashing scheme limitation
    stopifnot(n_bin < sqrt(.Machine$integer.max))

    hash_bin_range <- function(D_)
        (D_[,1] - 1L) + n_bin * (D_[,2] - 1L)

    unhash_bin_range <- function(x)
        1L + ((x %% n_bin)):(x %/% n_bin)

    # Check interval data are integer indices
    stopifnot(all(sapply(D, is.integer)))
    stopifnot(all(D >= 1L & D <= n_bin))
    stopifnot(all(D[,1L] <= D[,2L]))
    if (!(is.integer(n_sample)))
        n_sample <- as.integer(n_sample)

    h <- hash_bin_range(D)
    h_table <- table(h) * n_sample
    h_table_num <- as.integer(names(h_table))

    r <- matrix(NA_integer_, nrow=nrow(D), ncol=n_sample)

    for (j in seq_len(length(h_table))) {
        k <- unhash_bin_range(h_table_num[j])
        r[h_table_num[j] == h,] <- sample(k,
                                          size=h_table[j],
                                          replace=T,
                                          prob=pclm_fit$gamma[k])
    }

    r

}

