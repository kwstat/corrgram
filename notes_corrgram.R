# notes_corrgram.R

# ----------------------------------------------------------------------------

# GGally

lib(GGally)

ggpairs(iris)

ggpairs(auto[,c(3:14)]) # slow

# ----------------------------------------------------------------------------

gpairs is generalized pairs plot...very nice.  Based on grid graphics

libs("gpairs")
data(Leaves)
gpairs(Leaves) # Impressive speed for large splom

gpairs(Leaves[1:10], upper.pars = list(scatter = 'stats'),
         lower.pars = list(scatter = 'corrgram'),
         stat.pars = list(verbose = FALSE), gap = 0)
gpairs::corrgram(Leaves[,-33])

gpairs(iris) # Nice

head(auto)
gpairs(auto)  # has problems...

# ----------------------------------------------------------------------------

Portfolio Optimization With R/Rmetrics uses my corrgram as the basis for theirs

See the following page (through page 161)

http://books.google.com/books?id=keDLpq86DU4C&pg=PA407&dq=corrgram+%22r-project%22&num=8&client=internal-uds&cd=1&source=uds#v=snippet&q=wright&f=false

# ----------------------------------------------------------------------------

The book "R in Action" (Quick-R site) uses corrgram.  Page 284-287.

http://www.amazon.com/gp/product/1935182390/ref=as_li_ss_tl?ie=UTF8&tag=luisapiolaswe-20&linkCode=as2&camp=1789&creative=390957&creativeASIN=1935182390

# ----------------------------------------------------------------------------

See also: corrplot,

http://cos.name/wp-content/uploads/2009/12/An-Introduction-to-Matrix-Visualization-and-corrplot-Package.pdf

