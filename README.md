# fastTopics

[![R-CMD-check](https://github.com/stephenslab/fastTopics/workflows/R-CMD-check/badge.svg)](https://github.com/stephenslab/fastTopics/actions)
[![CircleCI](https://dl.circleci.com/status-badge/img/gh/stephenslab/fastTopics/tree/master.svg?style=svg)](https://app.circleci.com/pipelines/github/stephenslab/fastTopics?branch=master)
[![codecov](https://codecov.io/gh/stephenslab/fastTopics/branch/master/graph/badge.svg)](https://app.codecov.io/gh/stephenslab/fastTopics)

fastTopics is an R package implementing fast, scalable optimization
algorithms for fitting topic models and non-negative matrix 
factorizations to count data. The methods exploit the
[close relationship][vignette-close-relationship] between the topic
model and Poisson non-negative matrix factorization. The package also
provides tools to compare, annotate and visualize model fits,
including functions to create "structure plots" and functions to
identify distinctive features of topics. The fastTopics package is a
successor to the [CountClust package][countclust].

If you find a bug, or you have a question or feedback on this software,
please post an [issue][issues].

## Citing this work

If you find the fastTopics package or any of the source code in this
repository useful for your work, please cite:

> K. K. Dey, C. J. Hsiao and M. Stephens (2017). [Visualizing the
> structure of RNA-seq expression data using grade of membership 
> models.][countclust-paper] PLoS Genetics 13, e1006599.
>
> P. Carbonetto, A. Sarkar, Z. Wang and M. Stephens (2021).
> [Non-negative matrix factorization algorithms greatly improve topic
> model fits.][fasttopics-paper] arXiv 2105.13440.

If you used the `de_analysis` function in fastTopics, please cite:

> P. Carbonetto, K. Luo, A. Sarkar, A. Hung, K. Tayeb, S. Pott and
> M. Stephens (2023). [GoM DE: interpreting structure in sequence
> count data with differential expression analysis allowing for
> grades of membership.][singlecell-topics-paper]
> Genome Biology 24, 236.

## License

Copyright (c) 2019-2023, Peter Carbonetto and Matthew Stephens.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license].

## Quick Start

Install and load the package from CRAN:

```R
install.packages("fastTopics")
library(fastTopics)
```

Alternatively, install the latest version from GitHub:

```R
remotes::install_github("stephenslab/fastTopics")
library(fastTopics)
```

Note that installing the package will require a C++ compiler setup
that is appropriate for the version of R installed on your
computer. For details, refer to the documentation on the
[CRAN website][cran].

For guidance on using fastTopics to analyze gene expression data, see
the [single-cell RNA-seq vignette, part 1][vignette-scrnaseq-1] and
[part 2][vignette-scrnaseq-2].

Also, try running the small example that illustrates the fast model
fitting algorithms:

```R
example("fit_poisson_nmf")
```

See the [package documentation][pkgdown] for more information.

## Developer notes

To prepare the package for CRAN, remove both single-cell vignettes,
then run `R CMD build fastTopics` to build the source package.

This is the command used to check the package before submitting to
CRAN:

```r
library(rhub)
check_for_cran(".",show_status = TRUE,
  env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "false",
               `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = "true"))
```

## Credits

The fastTopics R package was developed by [Peter Carbonetto][peter],
[Matthew Stephens][matthew] and others.

[fasttopics]:  https://github.com/stephenslab/fastTopics
[mit-license]: https://opensource.org/license/mit
[issues]: https://github.com/stephenslab/fastTopics/issues
[peter]: https://pcarbo.github.io
[kevin]: https://github.com/kevinlkx
[matthew]: http://stephenslab.uchicago.edu
[uchicago]: https://www.uchicago.edu
[cran]: https://cran.r-project.org
[countclust]: https://github.com/kkdey/CountClust
[countclust-paper]: https://doi.org/10.1371/journal.pgen.1006599
[fasttopics-paper]: https://arxiv.org/abs/2105.13440
[singlecell-topics-paper]: https://doi.org/10.1186/s13059-023-03067-9
[pkgdown]: https://stephenslab.github.io/fastTopics/
[vignette-close-relationship]: https://stephenslab.github.io/fastTopics/articles/relationship.html
[vignette-scrnaseq-1]: https://stephenslab.github.io/fastTopics/articles/single_cell_rnaseq_basic.html
[vignette-scrnaseq-2]: https://stephenslab.github.io/fastTopics/articles/single_cell_rnaseq_practical.html
