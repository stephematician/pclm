# Penalised Composite Link Model tools in R

_Stephen Wade_

Measurements of data in surveys are often given in groups. It can be useful to
think of each group as representing a collection of smaller bins which we want
to estimate the distribution for.

An example of where this might be useful is if two surveys have different
groupings of age, and they may not be a random sample of a known distribution.
We can then estimate the distribution of age for each survey, and for each
individual in those surveys, using, say, 1-year age groups if we wish. We could
then compare the results for 1-year age groups across the surveys.

Penalised composite link models are an efficient method to ungrouping binned
data where there is no additional information about the distribution.
The method implemented in this code is adapted from the paper:

Silvia Rizzi, Jutta Gampe, and Paul H. C. Eilers, 2015, Efficient estimation
of smooth distributions from coarsely grouped data, _Am. J. Epidemiol.,
182_(2), pp. 138-147,
[doi.10.1093/aje/kwv020](https://doi.org/10.1093/aje/kwv020)

Further discussion of the penalty applied to the rates (or spline coefficients)
can be found here:

Paul H. C. Eilers and Brian D. Marx, 1996, Flexible smoothing with B-splines
and penalties, _Statistical Science, 11_(2), pp. 89-102,
[doi:10.1214/ss/1038425655](https://doi.org/10.1214/ss/1038425655)


## Installation

Installation is easy using
[`devtools`](https://cran.r-project.org/package=devtools)

```r
library(devtools)
install_github('stephematician/pclm')
```

## Example

```r
test_data <- rpois(10000, 10)
pois_hist <- hist(test_data, breaks=10, plot=F)

n <- pois_hist$counts
x <- with(pois_hist,
          sapply(min(breaks):(max(breaks)-1), 
                 function(x) sum(x >= breaks)))
pois_pclm <- pclmfit(n, x)

pois_ungroup_hist <- hist(test_data,
                          breaks=with(pois_hist, min(breaks):max(breaks)),
                          plot=F)

plot(pois_ungroup_hist)
points(pois_ungroup_hist$mids, pois_pclm$gamma)
```

