# corrgram

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/corrgram)](https://cran.r-project.org/package=corrgram)
[![Research software impact](http://depsy.org/api/package/cran/corrgram/badge.svg)](http://depsy.org/package/r/corrgram)

The `corrgram` package provides a simple way to create correlograms from raw data or a correlation matrix.

Key features:

* Stable, well-tested, widely used.

* Extensive examples show how to customize the display.

## Installation

```R
# Install the released version from CRAN:
install.packages("corrgram")

# Install the cutting edge development version from GitHub:
install.packages("devtools")
devtools::install_github("kwstat/corrgram")
```
## Example

Vignette:
[Lucid printing](https://rawgit.com/kwstat/corrgram/master/vignettes/lucid_printing.html)

```R
require(corrgram)
corrgram(mtcars, order=TRUE, lower.panel=panel.shade, upper.panel=panel.pie,
         text.panel=panel.txt, main="mtcars")
```
![corrgram](figure/corrgram_auto.png?raw=true)
