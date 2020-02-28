# fastTopics

[![Travis Build Status](https://travis-ci.org/stephenslab/fastTopics.svg?branch=master)](https://travis-ci.org/stephenslab/fastTopics)
[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/224272mhk5fadgmt?svg=true)](https://ci.appveyor.com/project/pcarbo/fasttopics)
[![codecov](https://codecov.io/gh/stephenslab/fastTopics/branch/master/graph/badge.svg)](https://codecov.io/gh/stephenslab/fastTopics)

R package implementing fast optimization algorithms for fitting topic
models and non-negative matrix factorizations to count data.

If you find a bug, or you have a question or feedback on this software,
please post an [issue][issues].

## License

Copyright (c) 2019-2020, Peter Carbonetto and Matthew Stephens.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license].

## Quick Start

Install and load the package:

```R
devtools::install_github("stephenslab/fastTopics")
library(fastTopics)
```

Note that installing the package will require a C++ compiler setup
that is appropriate for the version of R installed on your
computer. For details, refer to the [CRAN documentation][cran].

Simulate a 80 x 100 counts matrix:

```R
set.seed(1)
X <- simulate_count_data(80,100,k = 3)$X
```

Fit the Poisson NMF model using the multiplicative ("mu") and
sequential co-ordinate ascent ("scd") updates:

```R
fit1 <- fit_poisson_nmf(X,k = 3,numiter = 500,method = "mu")
fit2 <- fit_poisson_nmf(X,k = 3,numiter = 200,method = "scd",
                        control = list(numiter = 4))
```

Compare the improvement in the solution over time:

```R
plot_progress_poisson_nmf(list(mu = fit1,scd = fit2))
```

For more, work through the `fit_poisson_nmf` example:

```R
example("fit_poisson_nmf")
```

## Credits

The fastTopics R package was developed by [Peter Carbonetto][peter] at
the [University of Chicago][uchicago], with guidance from
[Matthew Stephens][matthew].

[mit-license]: https://opensource.org/licenses/mit-license.html
[issues]: https://github.com/stephenslab/fastTopics/issues
[peter]: https://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
[uchicago]: https://www.uchicago.edu
[cran]: https://cran.r-project.org

