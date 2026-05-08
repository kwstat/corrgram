# Panel functions for lattice corrgrams via splom

These functions provide custom panel methods for
[`lattice::splom()`](https://rdrr.io/pkg/lattice/man/splom.html).

## Usage

``` r
splom.pie(x, y, z, subscripts, at = pretty(z), cor.method = "pearson", ...)

splom.shade(
  x,
  y,
  z,
  subscripts,
  at = pretty(z),
  col.regions = NULL,
  cor.method = "pearson",
  ...
)

splom.ellipse(x, y, col.regions, cor.method = "pearson", ...)
```

## Arguments

- x:

  Numeric coordinates from levelplot.

- y:

  Numeric coordinates from levelplot.

- z:

  Correlation values from levelplot.

- subscripts:

  Subscripts for lattice panel. (not used)

- at:

  Breaks for color levels.

- cor.method:

  Correlation method (default "pearson").

- ...:

  Additional arguments passed to panel functions.

- col.regions:

  Color palette for shading (default NULL uses internal red white blue
  palette).

## Details

`splom.pie` Draws pie glyphs representing correlation coefficients,
omitting the diagonal.

`splom.shade` Draws shaded rectangles with hash lines representing
correlation coefficients.

`splom.ellipse` Draws ellipses representing correlation coefficients.
The position of the ellipse is determined by the position of the data,
and the shape of the ellipse is determined by the correlation.

## Examples

``` r
library(lattice)
pengvars <- c("bill_len", "bill_dep", "flipper_len", "body_mass")
splom(~penguins[ , pengvars], upper.panel=splom.pie, pscales=0)

splom(~penguins[ , pengvars]|penguins$species, upper.panel=splom.pie, pscales=0)

splom(~penguins[ , pengvars], upper.panel=splom.shade, pscales=0)

splom(~penguins[ , pengvars]|penguins$species, upper.panel=splom.shade, pscales=0)

splom(~penguins[ , pengvars], upper.panel=splom.ellipse, pscales=0)

splom(~penguins[ , pengvars]|penguins$species, upper.panel=splom.ellipse, pscales=0)
```
