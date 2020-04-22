# fastTopics

[![Travis Build Status](https://travis-ci.org/stephenslab/fastTopics.svg?branch=master)](https://travis-ci.org/stephenslab/fastTopics)
[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/224272mhk5fadgmt?svg=true)](https://ci.appveyor.com/project/pcarbo/fasttopics)
[![CircleCI build status](https://circleci.com/gh/stephenslab/fastTopics.svg?style=svg)](https://circleci.com/gh/stephenslab/fastTopics)
[![codecov](https://codecov.io/gh/stephenslab/fastTopics/branch/master/graph/badge.svg)](https://codecov.io/gh/stephenslab/fastTopics)

fastTopics is an R package implementing fast optimization algorithms
for fitting topic models and non-negative matrix factorizations to
count data. The methods exploit the
[close relationship][vignette-close-relationship] between topic
modeling and Poisson non-negative matrix factorization.

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

Simulate a (sparse) 80 x 100 counts matrix:

```R
library(Matrix)
set.seed(1)
X <- simulate_count_data(80,100,k = 3,sparse = TRUE)$X
```

Remove columns (words) that do not appear in any row (document).

```R
X <- X[,colSums(X > 0) > 0]
```

Fit the Poisson NMF model by running 200 updates (you may want to
start with fewer iterations and see how it goes):

```R
fit <- fit_poisson_nmf(X,k = 3,numiter = 200)
```

Recover the topic model. After this step, the L matrix contains
the topic probabilities ("loadings"), and the F matrix contains the
word probabilities ("factors").

```R
fit.multinom <- poisson2multinom(fit)
```

For more, see the examples included with `help(fit_poisson_nmf)`:

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
[vignette-close-relationship]: https://stephenslab.github.io/fastTopics/articles/relationship.html
