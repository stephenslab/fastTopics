# fastTopics

[![Travis Build Status](https://travis-ci.org/stephenslab/fastTopics.svg?branch=master)](https://travis-ci.org/stephenslab/fastTopics)
[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/224272mhk5fadgmt?svg=true)](https://ci.appveyor.com/project/pcarbo/fasttopics)
[![CircleCI build status](https://circleci.com/gh/stephenslab/fastTopics.svg?style=svg)](https://circleci.com/gh/stephenslab/fastTopics)
[![codecov](https://codecov.io/gh/stephenslab/fastTopics/branch/master/graph/badge.svg)](https://codecov.io/gh/stephenslab/fastTopics)

fastTopics is an R package implementing fast, scalable optimization
algorithms for fitting topic models ("grade of membership" models) and
non-negative matrix factorizations to count data. The methods exploit
the [special relationship][vignette-close-relationship] between the
multinomial topic model (also "probabilistic latent semantic
indexing") and Poisson non-negative matrix factorization. The package
provides tools to compare, annotate and visualize model fits,
including functions to efficiently create "structure plots" and
identify key features in topics. The fastTopics package is a successor
to the [CountClust package][countclust].

If you find a bug, or you have a question or feedback on this software,
please post an [issue][issues].

## Citing this work

If you find the fastTopics package or any of the source code in this
repository useful for your work, please cite:

> Kushal K. Dey, Chiaowen Joyce Hsiao and Matthew Stephens (2017).
> [Visualizing the structure of RNA-seq expression data using grade of membership models.][countclust-paper]
> *PLoS Genetics* **13**, e1006599.
>
> Peter Carbonetto, Kevin Luo, Kushal Dey, Joyce Hsiao and Matthew
> Stephens (2020). fastTopics: fast algorithms for fitting topic models
> and non-negative matrix factorizations to count data. R package
> version 0.3-194. [https://github.com/stephenslab/fastTopics][fasttopics]

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

For guidance on using fastTopics to analyze gene expression data
(single-cell or bulk RNA-seq), see the [single-cell RNA-seq
vignette][vignette-scrnaseq].

Also, try running the small example that illustrates the fast model
fitting algorithms:

```R
example("fit_poisson_nmf")
```

See the [package documentation][pkgdown] for more information.

## Credits

The fastTopics R package was developed by [Peter Carbonetto][peter],
Kevin Luo, Kushal Dey, Joyce Hsiao and
[Matthew Stephens][matthew] at the [University of Chicago][uchicago].

[fasttopics]:  https://github.com/stephenslab/fastTopics
[mit-license]: https://opensource.org/licenses/mit-license.html
[issues]: https://github.com/stephenslab/fastTopics/issues
[peter]: https://pcarbo.github.io
[kevin]: https://github.com/kevinlkx
[matthew]: http://stephenslab.uchicago.edu
[uchicago]: https://www.uchicago.edu
[cran]: https://cran.r-project.org
[countclust]: https://github.com/kkdey/CountClust
[countclust-paper]: https://doi.org/10.1371/journal.pgen.1006599
[pkgdown]: https://stephenslab.github.io/fastTopics
[vignette-close-relationship]: https://stephenslab.github.io/fastTopics/articles/relationship.html
[vignette-scrnaseq]: https://stephenslab.github.io/fastTopics/articles/single_cell_rnaseq_basic.html
