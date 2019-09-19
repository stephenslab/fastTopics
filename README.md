# fastTopics

R package implementing fast optimization algorithms for topic modeling
and non-negative matrix factorization.

If you find a bug, or you have a question or feedback on this software,
please post an [issue][issues].

## License

Copyright (c) 2019, Peter Carbonetto and Matthew Stephens.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license].

## Quick Start

To install the latest version of the fastTopics package, first
download or clone the repository from GitHub, then use
[devtools][devtools]:

```R
install.packages("devtools")
devtools::install_local("fastTopics")
```

This command should automatically install all required packages if
they are not installed already.

Compiling this package from source will require a C++ compiler setup
that is appropriate for the the R installed on your computer. For
details, refer to the [CRAN documentation][cran]. For Mac computers,
see [these notes][compiling-macos].

## Developer notes

### Testing the package

To install and test the fastTopics R package, run the following
commands in the shell:

```bash
R CMD build fastTopics
R CMD INSTALL fastTopics_0.1-75.tar.gz
R CMD check --as-cran fastTopics_0.1-75.tar.gz
```

Note that these commands require that the dependencies have already
been installed. See the [DESCRIPTION](DESCRIPTION) file for details.

### Updating the C++ source and documentation

When any changes are made to [roxygen2][roxygen2] markup or to the C++
code in the `src` directory, simply run `devtools::document()` to 
update the [RcppExports.cpp](src/RcppExports.cpp), the NAMESPACE file,
and the package documentation files in the `man` directory.

## Credits

The fastTopics R package was developed by [Peter Carbonetto][peter] at
the [University of Chicago][uchicago], with guidance from
[Matthew Stephens][matthew].

[issues]: https://github.com/stephenslab/fastTopics/issues
[peter]: https://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
[roxygen2]: https://cran.r-project.org/package=roxygen2
[devtools]: https://github.com/r-lib/devtools
[cran]: https://cran.r-project.org
[compiling-macos]: https://pcarbo.github.io/pcarbo/r-macos.html
