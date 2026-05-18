# corrgram

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/corrgram)](https://cran.r-project.org/package=corrgram)
[![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/corrgram)](https://cranlogs.r-pkg.org/badges/corrgram)

Homepage: <https://kwstat.github.io/corrgram>

Repository: <https://github.com/kwstat/corrgram>

The `corrgram` package provides a simple way to create correlograms from
raw data or a correlation matrix.

## Package overview

The `corrgram` package provides functions for creating corrgrams using
three different graphics systems, base, grid, and lattice.

Base R graphics

- single function
  [`corrgram()`](http://kwstat.github.io/corrgram/reference/corrgram.md)
  for dataframes or matrices.
- Enables most features found in the paper by @friendly2002corrgrams.
- No automatic legend.
- Not easily combined with other graphics.

`lattice` graphics

- Separate panel functions for
  [`lattice::levelplot()`](https://rdrr.io/pkg/lattice/man/levelplot.html)
  for dataframes and
  [`lattice::splom()`](https://rdrr.io/pkg/lattice/man/splom.html) for
  correlation matrices.
- Enables automatic legend.
- Enables corrgrams conditioned on other variables.
- Can be combined with other lattice graphics for complex figures.
- Not feature complete compared to base R.

`grid` graphics

- single function
  [`corrgram2()`](http://kwstat.github.io/corrgram/reference/corrgram2.md)
  for either dataframes or correlation matrices.
- Enables automatic legend.
- Can be combined with other grid graphics for complex figures.
- Not feature complete compared to base R.
- Faster than base R when evaluated inside Positron.

## Installation

``` r

# Install the released version from CRAN:
install.packages("corrgram")

# Install the development version from GitHub:
install.packages("devtools")
devtools::install_github("kwstat/corrgram")
```

## Usage for base R

``` r

require(corrgram)
# Base R
corrgram(mtcars, order=TRUE, lower.panel=panel.shade, upper.panel=panel.pie,
         text.panel=panel.txt, main="mtcars")
```

![corrgram](reference/figures/corrgram_mtcars.png)

corrgram

## Usage for grid version

``` r

# grid
vars6 <- setdiff(colnames(auto), c("Model", "Origin"))
corrgram2(auto[, vars6], 
          lower.panel = grid_panel.shade, upper.panel=grid_panel.pie,
          title = "auto data", legend = TRUE)
```

![grid corrgram](reference/figures/corrgram_grid.png)

grid corrgram

## Usage for lattice version

``` r

library(lattice)
pengvars <- c("bill_len", "bill_dep", "flipper_len", "body_mass")
splom(~penguins[ , pengvars]|penguins$species, upper.panel=splom_panel.pie, 
      pscales=0, main="penguins data in lattice")
```

![lattice corrgram](reference/figures/corrgram_lattice.png)

lattice corrgram
