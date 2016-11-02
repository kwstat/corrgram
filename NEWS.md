
# corrgram 1.10 - unpublished

`corrgram()` now returns the correlation matrix.

# corrgram 1.9 - Jul 2016

Change vignette from Rnw to Rmd.

Switch to Roxygen for documentation.

Add `panel.cor` for colored correlation values.

Fix minor bugs with no complete cases

Added more cases to tests directory.

# corrgram 1.8 - Jul 2015

Namespace changes due to R devel changes.

# corrgram 1.7 - Feb 2015

Support more seriate options.

# corrgram 1.6 - Aug 2014

Move packages from Depends to Imports.

Argument label.pos now defaults to c(.5, .5) for x,y positioning. (Request of Evan Williams).

New argument 'cor.method'. (Request of Evan Williams).

# corrgram 1.5 - Aug 2013

Update references links.  (Request of Michael Friendly).

Small bug with test for correlation matrix (Reported by F. Rosa)

# corrgram 1.4 - Nov 2012

Panel colors are now specified via a 'col.regions' argument.  The old
method of using 'col.corrgram' is ignored. Could cause slight
incompatability.

Re-worked examples section for more variety.

# corrgram 1.3 - Aug 2012

Added `panel.bar`.  (Request of dadrivr)

Added example for unclipped labels (in the test suite).

# corrgram 1.2

Small bug.  Now accepts NA values in correlation matrices, and
NAs caused by missing combinations of data and cor( use="pair").

Non-numeric columns in the data will be ignored. (Request of JZ)

New panel function panel.density

Test suite now includes corrgram of inVerse correlation, partial
correlation matrices.

# corrgram 1.1

Added namespace.

# corrgram 1.0

New panel function 'panel.conf'

New data set 'vote'

New function arguments 'dir' (alternative to 'row1attop'), 'type',
'label.srt', 'title', 'abs'.

Now works for either a data.frame or a correlation matrix using the
'type' argument.  Panel functions need a 'corr' argument.

New ordering method 'OLO' using 'seriation' package.

Ordering can now be done on absolute value of correlations.

Fixed bug: Strongly negative correlations used no color.

# corrgram 0.1 - Nov 2006

First release to CRAN.

# corrgram 0.0  - Apr 2006

Package development begins.
