
# corrgram 1.12 - May 2017

New vignette about an invalid correlation matrix.

`corrgram()` gives a better warning if a symmetric matrix has values outside [-1,1].\

Added package logo on github.

# corrgram 1.11 - Apr 2017

New argument `outer.labels` for adding labels along outside edges of corrgram. (Request of Vanessa Bruat and others.)

Test coverage 95% using `covr` package.

# corrgram 1.10 - Nov 2016

Function `corrgram()` now returns the correlation matrix.

Fixed custom label ordering when `order=TRUE` for M.Bruneaux.

Began using `testthat` package.

# corrgram 1.9 - Jul 2016

New panel function `panel.cor` for colored correlation values.

Fixed minor bugs with no complete cases.

Added more cases to tests directory.

# corrgram 1.8 - Jul 2015

Namespace changes due to R devel changes.

# corrgram 1.7 - Feb 2015

Added more seriate options.

# corrgram 1.6 - Aug 2014

Moved packages from Depends to Imports.

Argument `label.pos` now defaults to c(.5, .5) for x,y positioning. (Request of Evan Williams).

New argument `cor.method`. (Request of Evan Williams).

# corrgram 1.5 - Aug 2013

Updated references links.  (Request of Michael Friendly).

Fixed small bug with test for correlation matrix. (Reported by F. Rosa)

# corrgram 1.4 - Nov 2012

New argument `col.regions` to specify panel colors.

Re-worked examples section for more variety.

# corrgram 1.3 - Aug 2012

New panel function `panel.bar`.  (Request of dadrivr)

Added example for unclipped labels (in the test suite).

# corrgram 1.2 - Mar 2012

Small bug.  Now accepts NA values in correlation matrices, and
NAs caused by missing combinations of data and cor(use="pair").

Non-numeric columns in the data will be ignored. (Request of JZ)

New panel function `panel.density`.

Test suite now includes corrgram of inverse correlation, partial
correlation matrices.

# corrgram 1.1 - Oct 2011

Added namespace.

# corrgram 1.0 - Jul 2011

New panel function `panel.conf`.

New data set `vote`.

New arguments `dir`, `type`, `label.srt`, `title`, `abs`.

Now works for either a data.frame or a correlation matrix using the
`type` argument.

New ordering method `OLO` using `seriation` package.

Ordering can now be done on absolute value of correlations.

# corrgram 0.1 - Dec 2006

First release to CRAN.

# corrgram 0.0  - Apr 2006

Package development begins.
