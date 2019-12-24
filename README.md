# fastTopics

R package implementing fast optimization algorithms for fitting topic
models and non-negative matrix factorizations to count data.

If you find a bug, or you have a question or feedback on this software,
please post an [issue][issues].

## License

Copyright (c) 2019, Peter Carbonetto and Matthew Stephens.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license].

## Quick Start

The fastTopics package is currently undergoing major re-development in
the master branch, so it is recommended to install the version from
the v0.1 branch. To do this, clone the repository from GitHub, switch
to the v0.1 branch by running in the shell

```bash
git checkout v0.1
```

then use [devtools][devtools]:

```R
install.packages("devtools")
devtools::install_local("fastTopics")
```

This command should automatically install all required packages if
they are not installed already.

Alternatively, to install the latest (in-development) version of the
fastTopics package, download or clone the repository from GitHub, then
run:

```R
devtools::install_local("fastTopics")
```

Compiling this package from source will require a C++ compiler setup
that is appropriate for the version of R installed on your
computer. For details, refer to the [CRAN documentation][cran]. For
Mac computers, see [these notes][compiling-macos].

## Developer notes

### Testing and installing the package

To install and test the fastTopics R package, run the following
commands in the shell:

```bash
R CMD build fastTopics
R CMD INSTALL fastTopics_0.2-37.tar.gz
R CMD check --as-cran fastTopics_0.2-37.tar.gz
```

Note that these commands require that the dependencies have already
been installed. See the [DESCRIPTION](DESCRIPTION) file for details.

### Updating the pkgdown site

Run this line of R code to build the website (make sure you have an
Internet connection while running the code):

```R
pkgdown::build_site(mathjax = FALSE)
```

Version 1.1.0 of pkgdown was used.

### Building the Sweave PDF

To build the Sweave PDF, change your working directory to the
`inst/derivations` subdirectory, and run this R code to generate the
Latex from the Sweave document, then the PDF from the Latex:

```R
setwd("vignettes")
Sweave("algorithms.Rnw")
tools::texi2pdf("algorithms.tex")
```

Once the changes are committed, pushed to GitHub, and sync'ed to
Overleaf, the Latex on Overleaf will also be updated.

### Updating the C++ source and documentation

When any changes are made to [roxygen2][roxygen2] markup or to the C++
code in the `src` directory, simply run `devtools::document()` to 
update the [RcppExports.cpp](src/RcppExports.cpp), the NAMESPACE file,
and the package documentation files in the `man` directory.

## Credits

The fastTopics R package was developed by [Peter Carbonetto][peter] at
the [University of Chicago][uchicago], with guidance from
[Matthew Stephens][matthew].

[mit-license]: https://opensource.org/licenses/mit-license.html
[issues]: https://github.com/stephenslab/fastTopics/issues
[peter]: https://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
[roxygen2]: https://cran.r-project.org/package=roxygen2
[devtools]: https://github.com/r-lib/devtools
[cran]: https://cran.r-project.org
[compiling-macos]: https://pcarbo.github.io/pcarbo/r-macos.html
