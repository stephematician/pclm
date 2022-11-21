#' Randomly draw samples from fitted PCLM
#'
#' Generate a sample for each observation (or otherwise) of a variable taken
#' from an interval of a PCLM.
#'
#' Define any interval of bins (not necessarily the same as the original
#' groupings) for each observation and then draw a bin from the PCLM using the
#' probabilities obtained from `pclmfit`.
#'
#' The flexibility to define any interval does come at a small cost in terms of
#' performance. This is most useful when incorporating some kind of censoring
#' into a prediction that provides additional information to narrow down the
#' range of possible bins. A simple example might be when 'age at survey' is
#' grouped, but 'age at event' is not, and therefore the participant must be at
#' least the latter age (in bins).
#'
#' @param D a two-column data.frame (or matrix) of integers, each row represents
#'   an observation of an interval of bins, columns correspond to lower and
#'   upper bins of the interval (inclusive)
#' @param pclm_fit a list as returned by pclmfit
#' @param n_sample number of samples to draw for each observation
#' @return matrix of drawn bins with one row per observation, and columns
#'   correspond to samples
#'
#' @export
#' @md
rpclmfit <- function(D, pclm_fit, n_sample=1L) {

    n_bin <- length(pclm_fit$gamma)

  # check interval data are integer indices
    stopifnot(all(apply(D, 2, is.integer)))
    stopifnot(min(D[,1:2]) >= 1L && max(D[,1:2]) <= n_bin)
    stopifnot(min(D[,2L] - D[,1L]) >= 0)
    if (!(is.integer(n_sample)))
        n_sample <- as.integer(n_sample)

    bins <- matrix(NA_integer_, nrow=nrow(D), ncol=n_sample)

    where_each <- split(1:nrow(D), D[,1L])

    for (start_j in seq_along(where_each)[sapply(where_each, length) > 0]) {

      # find which rows have the current interval 'start', and get the number
      # of draws required for each interval 'end'
        where_j <- where_each[[start_j]]
        max_end_j <- max(D[where_j,2L])
        n_draw <- tabulate(D[where_j,2L], nbins=max_end_j)[start_j:max_end_j] *
                      n_sample

        for (delta_k in seq_along(n_draw)[n_draw != 0]) {
          # for each [start, end] pair (i.e. index j,k), sample the values using
          # the probability given by the PCLM
            is_k <- D[where_j,2L] == (start_j - 1) + delta_k
            bins[where_j[is_k],] <- (start_j - 1) + sample.int(
                delta_k,
                size=n_draw[delta_k],
                replace=T,
                prob=pclm_fit$gamma[(start_j - 1) + (1:delta_k)]
            )
        }

    }

   bins 

}

