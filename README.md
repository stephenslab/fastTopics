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

Install fastTopics using [devtools][devtools]:

```R
devtools::install_github("stephenslab/fastTopics")
```

Note that installing the package will require a C++ compiler setup
that is appropriate for the version of R installed on your
computer. For details, refer to the [CRAN documentation][cran].

Once you have installed the package, load the package in R:

```R
library(fastTopics)
```

Next, work through the `fit_poisson_nmf` example, which illustrates
the use of `fit_poisson_nmf` for fitting a non-negative matrix
factorization to a small (80 x 100) counts matrix:

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
[devtools]: https://github.com/r-lib/devtools
[cran]: https://cran.r-project.org

