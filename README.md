# fastTopics

R package implementing fast optimization algorithms for topic modeling
and non-negative matrix factorization.

## License

Copyright (c) 2019, Peter Carbonetto and Matthew Stephens.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license].

## Developer notes

### Testing the package

To install and test the mixsqp R package, run the following commands
in the shell:

```bash
R CMD build fastTopics
R CMD INSTALL fastTopics_0.1-6.tar.gz
R CMD check --as-cran fastTopics_0.1-6.tar.gz
```

Note that these commands require that the dependencies have already
been installed. See the [DESCRIPTION](DESCRIPTION) file for details.

### Updating the C++ source and documentation

When any changes are made to [roxygen2][roxygen2] markup or to the C++
code in the `src` directory, simply run `devtools::document()` to 
update the [RcppExports.cpp](src/RcppExports.cpp), the NAMESPACE file,
and the package documentation files in the `man` directory.

## Credits

The mixsqp R package was developed by [Peter Carbonetto][peter] at the
[University of Chicago][uchicago], with help from
[Matthew Stephens][matthew].

[peter]: https://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
[roxygen2]: https://cran.r-project.org/package=roxygen2
