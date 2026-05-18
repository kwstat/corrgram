# corrgram_grid.R

# For a correlation matrix, grid_panel.ellipse, grid_panel.pie do nothing.

# Currently there are 14 color bins, but this could be made an argument.
# The legend also has 14 bins, which is not ideal.

#' Draw a correlogram using `grid` graphics
#'
#' `corrgram2()` produces a correlogram using grid graphics. The off-diagonal cells
#' can be shaded or filled with custom panel functions to show the correlation
#' structure of a data matrix or correlation matrix.
#'
#' This function is a grid-graphics variant of [corrgram()]. It accepts either a
#' data matrix/data frame with one observation per row or a correlation matrix.
#' When raw data are supplied, correlations are computed with
#' `use = "pairwise.complete.obs"`.
#'
#' Variable reordering can be used to improve the display by placing related
#' variables near each other. `order = TRUE` and `order = "PC"` use the
#' PCA-based angular ordering described by Friendly (2002). `order =
#' "seriation"` uses `cba::seriation(..., method = "Optimal")` and requires the
#' `cba` package.
#'
#' @aliases corrgram2 grid_panel.conf grid_panel.ellipse grid_panel.fill
#' grid_panel.pie grid_panel.pts grid_panel.shade grid_text.panel
#' 
#' @param x A data frame or matrix with one observation per row, or a
#'   correlation matrix.
#' 
#' @param type Use 'data' or 'cor'/'corr' to explicitly specify whether
#'   `x` is raw data or a correlation matrix. Usually this is inferred.
#' 
#' @param order Should variables be reordered? Use `FALSE` for no reordering,
#'   `TRUE` or `"PC"` for PCA-based angular ordering, or `"seriation"` for
#'   optimal seriation via the `cba` package.
#' 
#' @param labels Labels to use on the diagonal instead of column names.
#' 
#' @param panel Default panel function used for both `lower.panel` and
#'   `upper.panel`.
#' 
#' @param ... Additional arguments passed to the panel functions.
#' 
#' @param lower.panel,upper.panel Separate panel functions used below and above
#'   the diagonal.
#' 
#' @param diag.panel Optional panel function used on the diagonal before drawing
#'   diagonal labels.
#' 
#' @param text.panel Included for API compatibility with [corrgram()]. Diagonal
#'   labels are currently drawn with [grid_text.panel()].
#' 
#' @param label.pos Horizontal and vertical placement of the diagonal label.
#' 
#' @param label.srt Rotation for diagonal labels.
#' 
#' @param cex.labels Label size for diagonal labels. Use `"fit"` to size labels
#'   to the panel width.
#' 
#' @param dir Direction of the main diagonal. Use `"left"` or `"\\"` for a
#'   descending diagonal, and `"right"` or `"/"` for an ascending diagonal.
#' 
#' @param legend If `TRUE`, draw a legend for the color scale.
#' 
#' @param col.regions A function returning a vector of colors.
#' 
#' @param cor.method Correlation method passed to panel functions. Default is
#'   `"pearson"`.
#' 
#' @param title Optional title drawn above the correlogram.
#' 
#' @param abs Logical; if `TRUE`, use absolute correlations for variable
#'   reordering.
#'
#' @return Invisibly returns `NULL`.
#'
#' @references
#' Friendly, Michael. 2002. Corrgrams: Exploratory Displays for Correlation
#' Matrices. *The American Statistician*, 56, 316--324.
#' \url{http://datavis.ca/papers/corrgram.pdf}
#'
#' D. J. Murdoch and E. D. Chow. 1996. A Graphical Display of Large Correlation
#' Matrices. The American Statistician, 50, 178-180.
#'
#' @examples
#' # Draw a grid-based correlogram from data
#' vars6 <- setdiff(colnames(auto), c("Model", "Origin"))
#' corrgram2(auto[vars6], order = TRUE,
#'     lower.panel = grid_panel.shade,
#'     upper.panel = grid_panel.pie)
#'
#' # 'vote' is a correlation matrix
#' corrgram2(vote, order = TRUE, upper.panel = grid_panel.conf)
#'
#' @keywords hplot
#' @import grid
#' @export
corrgram2 <- function(x,
                type=NULL,
                order=FALSE, 
                labels,
                panel=grid_panel.shade,
                ...,
                lower.panel = panel,
                upper.panel = panel,
                diag.panel = NULL,
                text.panel = grid_text.panel,
                label.pos = c(0.5, 0.5),
                label.srt=0,
                cex.labels="fit",
                dir="left", 
                legend=FALSE,
                col.regions=colorRampPalette(c("red","salmon","white","royalblue","navy")),
                cor.method="pearson",
                title=NULL,
                abs=FALSE) {
    
  # Direction
  if (dir=="\\") dir <-  "left"
  if (dir=="/") dir <-  "right"
  
  if (ncol(x) < 2) stop("Only one column in the argument to 'corrgram'")
  
  # if a matrix x has only colnames, isSymmetric(x) reports FALSE
  # use unname(x) so it will report TRUE
  if(is.matrix(x) && isSymmetric(unname(x))) {
    maybeCorr <- TRUE
    if( min(x, na.rm=TRUE) < -1-.Machine$double.eps ) {
      warning("Matrix symmetric, but some values < -1. Treated as data.frame.")
      maybeCorr <- FALSE
    }
    if( max(x, na.rm=TRUE) >  1+.Machine$double.eps ) {
      warning("Matrix symmetric, but some values > +1. Treated as data.frame.")
      maybeCorr <- FALSE
    }
  } else maybeCorr <- FALSE

  if(is.null(type)){
    type <- if(maybeCorr) "corr" else "data"
  } else if(type=="data"){
    if(maybeCorr)
      warning("This looks like a correlation matrix.")
  } else if(type=="cor" || type=="corr") {
    type <- "corr"
    if(!maybeCorr)
      warning("This is NOT a correlation matrix.")
  } else {
    stop("unknown data type in 'corrgram'")
  }

  # Remove non-numeric columns from data frames
  if(type=="data" && !is.matrix(x)) x <- x[ , vapply(x, is.numeric, logical(1))]

   # If a data matrix, then calculate the correlation matrix
  if(type=="data") {
    cmat <- cor(x, use="pairwise.complete.obs", method=cor.method)
  } else {
    cmat <- x
  }

  # Save the correlation matrix for returning to user, re-ordered below
  cmat.return <- cmat

  # should we use absolute correlations for determining ordering?
  cmat <- if(abs) abs(cmat) else cmat


 
  # If a data matrix, then calculate the correlation matrix
  x.cor <- if(type=="data") cor(x, use="pairwise.complete.obs") else x
  
  has.labs <- TRUE
  if (missing(labels)) {
    labels <- colnames(x)
    if (is.null(labels)) labels <- paste("var", 1:ncol(x))
  }
  else if(is.null(labels)) has.labs <- FALSE
  

  # default order
  if(order==FALSE){
    ord <- 1:nrow(cmat)
    # Re-order the data to group highly correlated variables
  } else if(order==TRUE || order=="PC" || order=="PCA"){
    # calculate the size of the angle between the horizontal axis and the
    # PC vectors
    x.eigen <- eigen(cmat)$vectors[,1:2]
    e1 <- x.eigen[,1]
    e2 <- x.eigen[,2]
    alpha <- ifelse(e1>0, atan(e2/e1), atan(e2/e1)+pi)
    ord <- order(alpha)
    x <- if(type=="data") x[,ord] else x[ord, ord]
    cmat.return <- cmat.return[ord, ord]
  } else {
    # Any other order method is delegated to seriation
    # The 'seriation' package is very heavy (many dependencies), so we
    # do not import it, but only check to see if it is installed.
    if (!("seriation" %in% rownames(installed.packages())))
      stop("Please use install.packages('seriation') for this 'order' option.")
    # "OLO" is used in this book
    # R Visualizations: Derive Meaning from Data By David Gerbing · 2020 p. 125
    distx <- as.dist(sqrt(1 - cmat))
    ss <- seriation::seriate(distx, method=order) # from seriation package
    ord <- seriation::get_order(ss)
    x <- if(type=="data") x[,ord] else x[ord,ord]
    cmat.return <- cmat.return[ord,ord]
  }
 
  # Don't pass some arguments on to the panel functions
  localLowerPanel <- function(..., main, oma, font.main, cex.main)
    lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main)
    upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main)
    diag.panel(...)
  #  localTextPanel <- function(...)
  
  #dots <- list(...)
  #nmdots <- names(dots)
   
  nc <- ncol(x)
   
  # New page
  grid::grid.newpage()
  grid::grid.rect(gp=grid::gpar(fill="light blue"))

  if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
    lower.panel <- match.fun(lower.panel)
  if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
    upper.panel <- match.fun(upper.panel)
  
  has.diag  <- !is.null(diag.panel)
  if(has.diag && !missing( diag.panel))
    diag.panel <- match.fun( diag.panel)
  
  if(dir=="left") {
    tmp <- lower.panel; lower.panel <- upper.panel; upper.panel <- tmp
    tmp <- has.lower; has.lower <- has.upper; has.upper <- tmp
  }
  
  # overall
  pushViewport(plotViewport(margins=c(1.1,1.1,3.1,1.1)))
  # Title
  if(!is.null(title))
    grid.text(title, x=unit(0.5, "npc"), 
                     y = unit(1, "npc") + unit(1.5, "lines"),
                     gp = gpar(fontsize = 16))
  
  # Splom. respect=TRUE keeps the panels square
  scatplot <- grid.layout(nc, nc, respect=TRUE)
  
  # Legend
  if(legend==TRUE){
    # Keep space for legend
    pushViewport(viewport(layout=grid.layout(1,2,widths=unit(c(.8,.2),
                                                             c("npc","npc")))))
    
    pushViewport(plotViewport(layout=scatplot,
                              layout.pos.col=1, layout.pos.row=1))
    grid.rect(gp=gpar(fill="blue"))
  } else {
    # No legend
    pushViewport(plotViewport(margins=c(1.1,1.1,4.1,2.1), layout=scatplot))
  }
  
  grid::grid.rect(gp=grid::gpar(fill="gray95"))
  #  grid.show.layout(scatplot)
  for(i in 1:nc) {
    for(j in 1:nc){
      
      xvals <- if(type=="data") as.vector(x[,j]) else NULL
      yvals <- if(type=="data") as.vector(x[,i]) else NULL
      corr <- if(type=="corr") x[j,i] else NULL
      
      colpos <- ifelse(dir=="left", j, nc-j+1)
      if(type=="data")
        pushViewport(dataViewport(xvals, yvals,
                                  layout.pos.row=i, layout.pos.col=colpos))
      else
        pushViewport(viewport(layout.pos.row=i, layout.pos.col=colpos))
      
      if(i==j) {
        # Diagonal panel (non-text)
        if (has.diag) {
          localDiagPanel(xvals, ...)
        }
        
        # Diagonal text
        if (has.labs) {
          if(cex.labels=="fit") {
            tw <- convertWidth(max(stringWidth(labels)), "inches",
                               valueOnly=TRUE)
            vw <- convertWidth(unit(1, "npc"), "inches", valueOnly=TRUE)
            cex.labels <- vw/tw
          }
          grid_text.panel(label.pos[1], label.pos[2], txt=labels[i], cex=cex.labels, rot=label.srt)
        }
        
      } else if (i<j) {
        # Lower panel
        if(type=="data")
          localLowerPanel(xvals, yvals, NULL, col.regions, cor.method, ...)
        else
          localLowerPanel(NULL, NULL, corr, col.regions, cor.method, ...)
      } else {
        # Upper panel
        if(type=="data")
          localUpperPanel(xvals, yvals, NULL, col.regions, cor.method, ...)
        else
          localUpperPanel(NULL, NULL, corr, col.regions, cor.method, ...)            
      }
      popViewport() # Panel
    }
  }
  popViewport() # Scatter plot
  
  # Draw legend at right
  if(legend) {
    n.colors <- 14
    legplot <- grid.layout(n.colors, 2)
    pushViewport(plotViewport(margins=c(1,1,1,3),
                              layout=legplot, layout.pos.row=1, layout.pos.col=2))
    #grid.rect(gp=gpar(fill="yellow"))
    pal <- col.regions(n.colors)
    midpoints <- seq(-1, 1, length.out=n.colors+1)[-1] - diff(seq(-1, 1, length.out=n.colors+1))[1]/2
    for(i in 1:n.colors) {
      # colored cells between break points
      pushViewport(viewport(layout.pos.row=i, layout.pos.col=1))
      grid.rect(gp=gpar(fill=pal[i]))
      popViewport()
      # labels for break points
      pushViewport(viewport(layout.pos.row=i, layout.pos.col=2))
      #grid.text(midpoints[i], 2, just="left", gp=gpar(cex=0.5))
      #grid.text(midpoints[i])
      grid.text(formatC(midpoints[i], format = "f", digits = 2))
      popViewport()
    }

    popViewport() # For legend
  }
  
  popViewport() # Entire graph

  invisible(cmat.return)
}


# ----------------------------------------------------------------------------

#' @export
grid_panel.fill <- function(x, y, corr=NULL, col.regions, cor.method="pearson") {
  # Same as grid_panel.shade, but without diagonal lines

  # if corr not given, try to calculate it
  if(is.null(corr)) {
    if(sum(complete.cases(x,y)) < 2) {
      warning("Need at least 2 complete cases for cor()")
      return()
    } else {
      corr <- cor(x, y, use='pair', method=cor.method)
    }
  }
  
  n_bins <- 14
  pal <- col.regions(n_bins)
  col.ind <- as.numeric(cut(corr, breaks=seq(from=-1, to=1, length.out=n_bins+1),
                            include.lowest=TRUE))
#browser()
  # Solid fill
  grid.rect(gp=gpar(col=pal[col.ind], fill=pal[col.ind]))
}
#corrgram2(auto[, vars6], abs=TRUE, panel=grid_panel.fill)

#' @export
grid_panel.shade <- function(x, y, corr=NULL, col.regions, cor.method="pearson") {
    # If corr not given, try to calculate it
  if(is.null(corr)) {
    if(sum(complete.cases(x,y)) < 2) {
      warning("Need at least 2 complete cases for cor()")
      return()
    } else {
      corr <- cor(x, y, use='pair', method=cor.method)
    }
  }

  n_bins <- 14
  pal <- col.regions(n_bins)
  col.ind <- as.numeric(cut(corr, breaks=seq(from=-1, to=1, length.out=n_bins+1),
                            include.lowest=TRUE))
  
  # Solid fill
  grid.rect(gp=gpar(col="lightgray", fill=pal[col.ind]))
  # Add diagonal lines in white
  if(corr>0){
    grid.polyline(c(0,1,.5,1,0,.5), c(0,1,0,.5,.5,1), id=c(1,1,2,2,3,3),
                  gp=gpar(col="white"))
  } else if(corr<0) {
    grid.polyline(c(1,0,.5,0,1,.5), c(0,1,0,.5,.5,1), id=c(1,1,2,2,3,3),
                  gp=gpar(col="white"))
  }
  
  # Bounding box needs to plot on top of the shading, so do it last.
  grid.rect(gp=gpar(col="lightgray", fill=NA))
}
#corrgram2(auto[, vars6], abs=TRUE, panel=grid_panel.shade)

#' @export
grid_panel.pts <- function(x, y, corr=NULL, col.regions, cor.method="pearson", ...){
  if(!is.null(corr)) return()
  # use hollow circles for the points
  grid.points(x, y, pch=1, gp = gpar(cex=0.5, fill=NA))
  grid.rect(gp=gpar(col="lightgray",fill=NA))
}
#corrgram2(auto[, vars6], abs=TRUE, panel=grid_panel.pts) # broken 

  
#' @export
grid_panel.pie <- function(x, y, corr=NULL, col.regions, cor.method="pearson"){

  # if corr not given, try to calculate it
  if(is.null(corr)) {
    if(sum(complete.cases(x,y)) < 2) {
      warning("Need at least 2 complete cases for cor()")
      return()
    } else {
      corr <- cor(x, y, use='pair', method=cor.method)
    }
  }
  
  # Start with the 'pie' piece
  # Multiply the radius by .97 so the circles do not overlap
  rad <- 0.5 * .97

  n_bins <- 14
  pal <- col.regions(n_bins)
  col.ind <- as.numeric(cut(corr, breaks=seq(from=-1, to=1, length.out=n_bins+1),
                            include.lowest=TRUE))
  col.pie <- pal[col.ind]
  
  numseg <- round(60*abs(corr),0)
  if(numseg>0){ # Watch out for the case with 0 numseg
    angles <- seq(pi/2, pi/2+(2*pi* -corr), length=numseg)
    circx <- c(0.5 + cos(angles)*rad, 0.5)
    circy <- c(0.5 + sin(angles)*rad, 0.5)
    grid.polygon(unit(circx, "npc"), unit(circy, "npc"),
                 gp=gpar(col=col.pie, fill=col.pie))
  }   
  
  # overlay the circle
  numseg <- 60
  angles <- seq(0, 2*pi, length=numseg)
  circx <- 0.5 + cos(angles)*rad
  circy <- 0.5 + sin(angles)*rad
  grid.polygon(unit(circx, "npc"), unit(circy, "npc"),
               gp=gpar(col="gray30", fill=NA))
  
}
#corrgram2(auto[, vars6], abs=TRUE, panel=grid_panel.pie) # NOT right fixme

#' @export
grid_panel.conf <- function(x, y, corr=NULL, col.regions, cor.method="pearson", ...){

  # if corr not given, try to calculate it
  if(is.null(corr)) {
    if(sum(complete.cases(x,y)) < 2) {
      warning("Need at least 2 complete cases for cor()")
      return()
    } else {
      corr <- cor(x, y, use='pair', method=cor.method)
    }
  }
  
  n_bins <- 14
  pal <- col.regions(n_bins)
  col.ind <- as.numeric(cut(corr, breaks=seq(from=-1, to=1, length.out=n_bins+1),
                            include.lowest=TRUE))
  col.pie <- pal[col.ind]
  
  if(is.null(x)) { # No data x, cannot calculate conf int
    grid.text(formatC(corr, 2, format='f'), gp=gpar(col=pal[col.ind]))      
  } else {
    if (corr < 1 - 1e-09) {
      conf.level <- .95
      z <- atanh(corr)
      b <- qnorm((1 - conf.level)/2)/sqrt(length(x) - 3)
      ci.z <- c(z + b, z - b)
      confint <- tanh(ci.z)
      confCol <- "black"
      if(confint[1]>0 | confint[2]<0) confCol <- pal[col.ind]
      confint <- formatC(confint, 2, format='f')
      confint <- paste("(",confint[1], ",", confint[2], ")",sep="")
    } else {
      confint <- "(NA,1)"
    }
    grid.text(formatC(corr, 2, format='f'), gp=gpar(col=pal[col.ind]))      
    grid.text(confint, y=.25, gp=gpar(col=confCol, cex=.5))
  }
}
#corrgram2(auto[, vars6], abs=TRUE, panel=grid_panel.conf)

#' @export
grid_panel.ellipse <- function(x, y, corr=NULL, col.regions, cor.method="pearson", ...){

  # For correlation matrix, do nothing
  if(!is.null(corr)) return()

  # If too few points, do nothing
  if(sum(complete.cases(x,y)) < 2) {
    warning("Need at least 2 complete cases to draw ellipse.")
    return()
  }

  # Min/Max of bounding box
  xrng <- extendrange(x, f=.05)
  yrng <- extendrange(y, f=.05)
  # Remove missing values
  keep <- (!is.na(x) & !is.na(y))
  x <- x[keep]
  y <- y[keep]
  # Draw an ellipse
  dfn <- 2
  dfd <- length(x)-1
  shape <- var(cbind(x,y))
  center <- c(mean(x),mean(y))
  radius <- sqrt(dfn*qf(.68,dfn,dfd))
  segments <- 75
  angles <- seq(0,2*pi,length=segments)
  unit.circle <- cbind(cos(angles),sin(angles))
  ellipse.pts <- t(center+radius*t(unit.circle%*%chol(shape)))
  ellx <- ellipse.pts[,1]
  elly <- ellipse.pts[,2]
  # Truncate ellipse at bounding box
  ellx <- ifelse(ellx < xrng[1], xrng[1], ellx) 
  ellx <- ifelse(ellx > xrng[2], xrng[2], ellx)
  elly <- ifelse(elly < yrng[1], yrng[1], elly)
  elly <- ifelse(elly > yrng[2], yrng[2], elly)
  grid.lines(unit(ellx, "native"), unit(elly, "native"),
             gp=gpar(col="gray30"))
  
  # Add a lowess line through the ellipse (if we have data)
  if(!is.null(x)){
    low <- stats::lowess(x, y, f=2/3, iter=3)
    grid.lines(unit(low$x,"native"), unit(low$y,"native"),
               gp=gpar(col="red"))
  }
  
}
#corrgram2(auto[, vars6], abs=TRUE, panel=grid_panel.ellipse)

#' @export
grid_text.panel <- function(x=0.5, y=0.5, txt, cex=1, rot=0){
  # diagonal panel labels
  grid.text(txt, x=x, y=y, just=c(.5,.5), rot=rot, gp=gpar(cex=cex))
}

# ----------------------------------------------------------------------------


if(FALSE){
  #data(auto, package="corrgram")
  #data(vote, package="corrgram")
  library(corrgram)
  vars6 <- setdiff(colnames(auto), c("Model", "Origin"))
  
  # compare speed of corrgram and corrgram2. 
  # grid version is much faster in Positron, but not in Rgui. 
  corrgram::corrgram(auto[,vars6], panel=corrgram::panel.pts)
  corrgram2(auto[, vars6], abs=TRUE, panel=grid_panel.pts)
  
  corrgram2(auto[, vars6], abs=TRUE, panel=grid_panel.fill, title="auto data")
  corrgram2(auto[, vars6], abs=TRUE, panel=grid_panel.shade)
  corrgram2(auto[, vars6], abs=TRUE, panel=grid_panel.pie)
  corrgram2(auto[, vars6], abs=TRUE, panel=grid_panel.conf)
  corrgram2(auto[, vars6], abs=TRUE, panel=grid_panel.ellipse)
  
  corrgram2(vote,abs=TRUE)
  corrgram2(vote,abs=TRUE, labels=1:12)
  corrgram2(vote,abs=TRUE, label.srt=0)
  corrgram2(vote,abs=TRUE, label.srt=60, cex.labels=1.5)
  corrgram2(vote,abs=TRUE,order=TRUE, legend=TRUE)
  corrgram2(vote, order="OLO", lower.panel=grid_panel.pie)
  corrgram2(vote, panel=grid_panel.conf)
  corrgram2(vote, panel=grid_panel.ellipse)
  corrgram2(vote, lower.panel=grid_panel.fill)
  corrgram2(vote, panel=grid_panel.pie, order=TRUE)
  corrgram2(vote, panel=grid_panel.pts) # fixme should err
  corrgram2(vote, upper.panel=grid_panel.shade)
  corrgram2(vote, lower.panel=grid_panel.fill, upper.panel=grid_panel.shade)
    
}