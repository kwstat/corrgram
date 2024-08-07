# corrgram 1.15 unpublished

* Switch to MIT license


# corrgram 1.14 (2021-04-01)

* The `seriation` package is moved from Imports to Suggests.  This makes `corrgram` lighter loading for most people. (Requested by P.Kiener and others)

* The `outer.labels` argument has improved default values.


# corrgram 1.13 (2018-07-09)

* Test coverage at 95% using `covr` package.

* New `panel.fill()` function, omits diagonal lines.

* `panel.conf` and `panel.cor` now auto-scale based on absolute value of correlations.

* `panel.conf` now only allows `cor.method="pearson"` (which is the default.) Fix issue #13.


# corrgram 1.12 (2017-05-07)

* New vignette about an invalid correlation matrix.

* `corrgram()` gives a better warning if a symmetric matrix has values outside [-1,1].

* Added package logo on github.


# corrgram 1.11 (2017-04-03)

* New argument `outer.labels` for adding labels along outside edges of corrgram. (Request of Vanessa Bruat and others.)


# corrgram 1.10 (2016-11-09)

* Function `corrgram()` now returns the correlation matrix.

* Fixed custom label ordering when `order=TRUE` for M.Bruneaux.

* Began using `testthat` package.


# corrgram 1.9 (2016-07-16)

* New panel function `panel.cor` for colored correlation values.

* Fixed minor bugs with no complete cases.

* Added more cases to tests directory.


# corrgram 1.8 (2015-07-03)

* Namespace changes due to R devel changes.


# corrgram 1.7 (2015-02-13)

* Added more seriate options.


# corrgram 1.6 (2014-08-29)

* Moved packages from Depends to Imports.

* Argument `label.pos` now defaults to c(.5, .5) for x,y positioning. (Request of Evan Williams).

* New argument `cor.method`. (Request of Evan Williams).


# corrgram 1.5 (2013-08-29)

* Updated references links.  (Request of Michael Friendly).

* Fixed small bug with test for correlation matrix. (Reported by F. Rosa)


# corrgram 1.4 (2012-11-07)

* New argument `col.regions` to specify panel colors.

* Re-worked examples section for more variety.


# corrgram 1.3 (2012-08-15)

* New panel function `panel.bar`.  (Request of dadrivr)

* Added example for unclipped labels (in the test suite).


# corrgram 1.2 (2012-03-28)

* Small bug.  Now accepts NA values in correlation matrices, and NAs caused by missing combinations of data and cor(use="pair").

* Non-numeric columns in the data will be ignored. (Request of JZ)

* New panel function `panel.density`.

* Test suite now includes corrgram of inverse correlation, partial correlation matrices.


# corrgram 1.1 (2011-10-20)

* Added namespace.


# corrgram 1.0 (2011-07-02)

* New panel function `panel.conf`.

* New data set `vote`.

* New arguments `dir`, `type`, `label.srt`, `title`, `abs`.

* Now works for either a data.frame or a correlation matrix using the `type` argument.

* New ordering method `OLO` using `seriation` package.

* Ordering can now be done on absolute value of correlations.


# corrgram 0.1 (2006-12-01)

* First release to CRAN.


# corrgram 0.0  (2006-04-01)

* Package development begins.
