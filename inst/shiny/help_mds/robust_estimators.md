## Robust Estimators

Check this option to use a cluster-robust sandwich estimate of the variance-covariance matrix of the regression coefficient estimates, allowing for corrections for heteroscedasticity-consistent variance estimators. 

In parallel, several small sample corrections are available:

- `CR0`: is the original form of the sandwich estimator ([Liang & Zeger, 1986](https://academic.oup.com/biomet/article/73/1/13/246001)), which does not make any small-sample correction.
- `CR1`: multiplies CR0 by m / (m - 1), where m is the number of clusters.
- `CR1p`: multiplies CR0 by m / (m - p), where m is the number of clusters and p is the number of covariates.
- `CR1S`: multiplies CR0 by (m (N-1)) / [(m - 1)(N - p)], where m is the number of clusters, N is the total number of observations, and p is the number of covariates. Some Stata commands use this correction by default.
- `CR2`: is the "bias-reduced linearization" adjustment proposed by [Bell and McCaffrey (2002)](https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X20020029058) and further developed in Pustejovsky and [Tipton (2017)](https://www.tandfonline.com/doi/full/10.1080/07350015.2016.1247004). The adjustment is chosen so that the variance-covariance estimator is exactly unbiased under a user-specified working model.
- `CR3`: approximates the leave-one-cluster-out jackknife variance estimator ([Bell & McCaffrey, 2002](https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X20020029058)).

`CR2` small sample correction is recommended and used by default.

More information can be found in [`vcovCR`](https://search.r-project.org/CRAN/refmans/clubSandwich/html/vcovCR.html) help page from the [`clubSandwich`](https://cran.r-project.org/web/packages/clubSandwich/index.html) package.

