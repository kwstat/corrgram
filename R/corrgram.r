# corrgram.r
# Time-stamp: <05 May 2017 16:04:24 c:/x/rpack/corrgram/R/corrgram.R>

# color key
# http://stackoverflow.com/questions/9852343/how-to-add-a-color-key-to-a-pairs-plot

# To do: Add a legend/ribbon.  See 'corrplot' package for a different
# way to calculate the layout.

# hexbin
# https://procomun.wordpress.com/2011/03/18/splomr/

# ----------------------------------------------------------------------------

#' Draw a correlogram
#' 
#' The corrgram function produces a graphical display of a correlation matrix,
#' called a correlogram.  The cells of the matrix can be shaded or colored to
#' show the correlation value.
#' 
#' Note: Use the 'col.regions' argument to specify colors.
#' 
#' Non-numeric columns in the data will be ignored.
#' 
#' The off-diagonal panels are specified with \code{panel.pts},
#' \code{panel.pie}, \code{panel.shade}, \code{panel.bar},
#' \code{panel.ellipse}, \code{panel.conf}. \code{panel.cor}.
#' 
#' Diagonal panels are specified with \code{panel.txt}, \code{panel.minmax},
#' \code{panel.density}.
#' 
#' Use a NULL panel to omit drawing the panel.
#' 
#' This function is basically a modification of the \code{pairs.default}
#' function with the use of customized panel functions.
#' 
#' The panel.conf function uses \code{cor.test} and calculates pearson
#' correlations.  Confidence intervals are not available in \code{cor.test} for
#' other methods (kendall, spearman).
#' 
#' You can create your own panel functions by starting with one of the included
#' panel functions and making suitable modifications.  Note that because of the
#' way the panel functions are called inside the main function, your custom
#' panel function must include the arguments shown in the \code{panel.pts}
#' function, even if the custom panel function does not use those arguments!
#' 
#' TODO: legend, grid graphics version.
#' 
#' @aliases corrgram panel.bar panel.conf panel.cor panel.density panel.ellipse
#' panel.minmax panel.pie panel.pts panel.shade panel.txt
#' 
#' @param x A \emph{tall} data frame with one observation per row, or a
#' correlation matrix.
#' 
#' @param type Use 'data' or 'cor'/'corr' to explicitly specify that 'x' is
#' data or a correlation matrix.  Rarely needed.
#' 
#' @param order Should variables be re-ordered?  Use TRUE/"PCA" for PCA-based
#' re-ordering.  Options from the 'seriate' package include "OLO" for optimal
#' leaf ordering, "GW", and "HC".
#' 
#' @param labels Labels to use (instead of data frame variable names) for
#' diagonal panels. If 'order' option is used, this vector of labels will be
#' also be appropriately reordered by the function.
#' 
#' @param panel Function used to plot the contents of each panel.
#' 
#' @param lower.panel,upper.panel Separate panel functions used below/above the
#' diagonal.
#' 
#' @param diag.panel,text.panel Panel function used on the diagonal.
#' 
#' @param label.pos Horizontal and vertical placement of label in diagonal
#' panels.
#' 
#' @param label.srt String rotation for diagonal labels.
#' 
#' @param cex.labels,font.labels Graphics parameter for diagonal panels.
#' 
#' @param row1attop TRUE for diagonal like " \ ", FALSE for diagonal like " / ".
#' 
#' @param dir Use \code{dir="left"} instead of 'row1attop'.
#' 
#' @param gap Distance between panels.
#' 
#' @param abs Use absolute value of correlations for clustering?  Default FALSE.
#' 
#' @param col.regions A \emph{function} returning a vector of colors.
#' 
#' @param cor.method Correlation method to use in panel functions.  Default is
#' 'pearson'.  Alternatives: 'spearman', 'kendall'.
#'
#' @param outer.labels A list of the form 'list(bottom,left,top,right)',
#' each component of which is a list of the form 'list(labels,cex,srt)'.
#' This is used to add labels along the outside edges of the corrgram.
#' Defaults: 'cex=1', 'srt=90' (bottom/top), 'srt=0' (left/right).
#' 
#' @param ... Additional arguments passed to plotting methods.
#' 
#' @return The correlation matrix used for plotting is returned. The 'order' and 'abs'
#' arguments affect the returned value.
#' 
#' @author Kevin Wright
#' 
#' @references
#' Friendly, Michael.  2002.  Corrgrams: Exploratory Displays for
#' Correlation Matrices.  \emph{The American Statistician}, 56, 316--324.
#' \url{http://datavis.ca/papers/corrgram.pdf}
#' 
#' D. J. Murdoch and E. D. Chow. 1996. A Graphical Display of Large Correlation
#' Matrices.  The American Statistician, 50, 178-180.
#' 
#' @keywords hplot
#' 
#' @examples
#' 
#' # To reproduce the figures in Michael Friendly's paper, see the
#' # vignette, or see the file 'friendly.r' in this package's
#' # test directory.
#' 
#' # Demonstrate density panel, correlation confidence panel
#' corrgram(iris, lower.panel=panel.pts, upper.panel=panel.conf,
#'          diag.panel=panel.density)
#' 
#' # Demonstrate panel.shade, panel.pie, principal component ordering
#' vars2 <- c("Assists","Atbat","Errors","Hits","Homer","logSal",
#'            "Putouts","RBI","Runs","Walks","Years")
#' corrgram(baseball[vars2], order=TRUE, main="Baseball data PC2/PC1 order",
#'          lower.panel=panel.shade, upper.panel=panel.pie)
#' 
#' # CAUTION: The latticeExtra package also has a 'panel.ellipse' function
#' # that clashes with the same-named function in corrgram. In order to use
#' # the right one, the example below uses 'lower.panel=corrgram::panel.ellipse'.
#' # If you do not have latticeExtra loaded, you can just use
#' # 'lower.panel=panel.ellipse'.
#' 
#' # Demonstrate panel.bar, panel.ellipse, panel.minmax, col.regions
#' corrgram(auto, order=TRUE, main="Auto data (PC order)",
#'          lower.panel=corrgram::panel.ellipse,
#'          upper.panel=panel.bar, diag.panel=panel.minmax,
#'          col.regions=colorRampPalette(c("darkgoldenrod4", "burlywood1",
#'                                         "darkkhaki", "darkgreen")))
#' 
#' # 'vote' is a correlation matrix, not a data frame
#' corrgram(vote, order=TRUE, upper.panel=panel.cor)
#'
#'
#' # outer labels, all options, larger margins, xlab, ylab
#' labs=colnames(state.x77)
#' corrgram(state.x77, oma=c(7, 7, 2, 2),
#'          outer.labels=list(bottom=list(labels=labs,cex=1.5,srt=60),
#'                            left=list(labels=labs,cex=1.5,srt=30)))
#' mtext("Bottom", side=1, cex=2, line = -1.5, outer=TRUE, xpd=NA)
#' mtext("Left", side=2, cex=2, line = -1.5, outer=TRUE, xpd=NA)
#' 
#' @import graphics
#' @import grDevices
#' @import seriation
#' @import stats
#' @export corrgram
corrgram <- function (x, type=NULL,
            order=FALSE, labels, panel = panel.shade,
            lower.panel = panel, upper.panel = panel,
            diag.panel = NULL, text.panel = textPanel,
            label.pos = c(0.5, 0.5), label.srt=0,
            cex.labels = NULL, font.labels = 1,
            row1attop = TRUE, dir="", gap = 0,
            abs=FALSE,
            col.regions = colorRampPalette(c("red","salmon","white","royalblue","navy")),
            cor.method="pearson",
            outer.labels=NULL,
            ...) {
  # Need graphics
  # Need grDevices
  # Need seriation for seriate
  # Need stats for cor, qf
    
  if(is.null(order)) order <- FALSE

  # Former versions used label.pos=0.5 for vertical positioning.
  if(length(label.pos) < 2) stop("label.pos needs a vector of length 2")

  # Direction
  if(dir=="") {
    dir <- if(row1attop) "left" else "right"
  }
  if (dir=="\\") dir <-  "left"
  if (dir=="/") dir <-  "right"

  if (ncol(x) < 2) stop("Only one column in the argument to 'corrgram'")

  # Do we have a data.frame or correlation matrix?
  # Note: Important to use "<=" instead of "<" (for example).
  ## if(is.matrix(x) && isSymmetric(x) &&
  ##    min(x, na.rm=TRUE) >= -1 - .Machine$double.eps &&
  ##    max(x, na.rm=TRUE) <= 1 + .Machine$double.eps)
  ##   maybeCorr <- TRUE
  ## else
  ##   maybeCorr <- FALSE
  
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
  cmat <- if(type=="data")
    cor(x, use="pairwise.complete.obs", method=cor.method)
  else
    x

  # Save the correlation matrix for returning to user, re-ordered below
  cmat.return <- cmat

  # should we use absolute correlations for determining ordering?
  cmat <- if(abs) abs(cmat) else cmat

  # Default order
  if(order==FALSE) ord <- seq_len(nrow(cmat))
  # Re-order the data to group highly correlated variables
  if(order==TRUE || order=="PC" || order=="PCA"){
    # Order by angle size between PCAs (first two) of correlation matrix
    x.eigen <- eigen(cmat)$vectors[,1:2]
    e1 <- x.eigen[,1]
    e2 <- x.eigen[,2]
    alpha <- ifelse(e1>0, atan(e2/e1), atan(e2/e1)+pi)
    ord <- order(alpha)
    x <- if(type=="data") x[,ord] else x[ord, ord]
    cmat.return <- cmat.return[ord,ord]
  } else if (order=="OLO") {
    distx <- dist(cmat)
    ss <- seriate(distx, method="OLO") # from seriation package
    ord <- get_order(ss)
    x <- if(type=="data") x[,ord] else x[ord,ord]
    cmat.return <- cmat.return[ord,ord]
  } else if (order=="GW"){ # GW order
    distx <- dist(cmat)
    ss <- seriate(distx, method="GW")
    ord <- get_order(ss)
    x <- if(type=="data") x[,ord] else x[ord,ord]
    cmat.return <- cmat.return[ord,ord]
  } else if (order=="HC"){ # HC ... just for comparision really
    distx <- dist(cmat)
    ss <- seriate(distx, method="HC")
    ord <- get_order(ss)
    x <- if(type=="data") x[,ord] else x[ord,ord]
    cmat.return <- cmat.return[ord,ord]
  } else if(order!=FALSE){
    stop("Unknown order argument in 'corrgram'.")
  }

  textPanel <- function(x = 0.5, y = 0.5, txt, cex, font, srt) {
    text(x, y, txt, cex=cex, font=font, srt=srt)
  }

  localAxis <- function(side, x, y, xpd, bg, col=NULL, main, oma, ...) {
    ## Explicitly ignore any color argument passed in as
    ## it was most likely meant for the data points and
    ## not for the axis.
    if(side %%2 == 1) Axis(x, side=side, xpd=NA, ...)
    else Axis(y, side=side, xpd=NA, ...)
  }

  # Don't pass some arguments on to the panel functions via the '...'
  localPlot <- function(..., main, oma, font.main, cex.main)
    plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main)
    lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main)
    upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main)
    diag.panel(...)

  dots <- list(...)
  nmdots <- names(dots)

  # Check for non-numeric data
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for(i in seq(along=names(x))) {
      if(is.factor(x[[i]]) || is.logical(x[[i]]))
        x[[i]] <- as.numeric(x[[i]])
      if(!is.numeric(unclass(x[[i]])))
        stop("non-numeric argument to 'corrgram'")
    }
  } else if (!is.numeric(x)) stop("non-numeric argument to 'corrgram'")

  # Get panel functions
  panel <- match.fun(panel)
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

  nc <- ncol(x)
  has.labs <- TRUE
  if (missing(labels)) {
    labels <- colnames(x)
    if (is.null(labels)) labels <- paste("var", seq_len(nc))
  } else if (order!=FALSE) {
      labels = labels[ord]
  } else if(is.null(labels)) has.labs <- FALSE
  if(is.null(text.panel)) has.labs <- FALSE

  oma <- if("oma" %in% nmdots) dots$oma else NULL
  main <- if("main" %in% nmdots) dots$main else NULL
  
  # Plot layout

  if (is.null(oma)) {
    oma <- c(4, 4, 4, 4)
    if (!is.null(main)) oma[3] <- 6 # Space for the title
  }
  # corrgram uses mfrow, plotting down the columns
  opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
  on.exit(par(opar))

  # Main loop to draw each panel
  for (i in if(dir=="left") seq_len(nc) else seq.int(from = nc, to = 1))
    for (j in seq_len(nc)) {
      # Set up plotting area
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, type = "n", ...)
      if(i == j || (i < j && has.lower) || (i > j && has.upper) ) {

        if(i == j) {
          # Diagonal panel
          if (has.diag) {
            if(type=="data")
              localDiagPanel(as.vector(x[, i]), NULL, ...)
            else
              localDiagPanel(NULL, x[i,i], ...)
          }

          # Diagonal text
          if (has.labs) {
            par(usr = c(0, 1, 0, 1))
            if(is.null(cex.labels)) {
              l.wid <- strwidth(labels, "user")
              cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
            }
            text.panel(label.pos[1], label.pos[2], labels[i],
                       cex = cex.labels, font = font.labels, srt=label.srt)
          }
        } else if(i < j) { # Lower panel
          if(type=="data")
            localLowerPanel(as.vector(x[, j]), as.vector(x[, i]), NULL, col.regions, cor.method, ...)
          else
            localLowerPanel(NULL, NULL, x[j,i], col.regions, cor.method, ...)
        } else { # Upper panel
          if(type=="data")
            localUpperPanel(as.vector(x[, j]), as.vector(x[, i]), NULL, col.regions, cor.method, ...)
          else
            localUpperPanel(NULL, NULL, x[j,i], col.regions, cor.method, ...)
        }

      } else { # No panel drawn
        par(new = FALSE)
      }

    }

  if (!is.null(main)) {
    font.main <- if("font.main" %in% nmdots) dots$font.main else par("font.main")
    cex.main <- if("cex.main" %in% nmdots) dots$cex.main else par("cex.main")
    mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
  }

  # ----------------------------------------------------------------------------
  
  corrgram.outer.labels(1, nc, ord, outer.labels$bottom)
  corrgram.outer.labels(2, nc, ord, outer.labels$left)
  corrgram.outer.labels(3, nc, ord, outer.labels$top)
  corrgram.outer.labels(4, nc, ord, outer.labels$right)
  
  ## # Add overall labels
  ## mtext(xlab, side = 1, line = -1.5, outer=TRUE, xpd=NA)
  ## mtext(ylab, side = 2, line = -1.5, outer=TRUE, xpd=NA)

  invisible(cmat.return)
}

corrgram.outer.labels <- function(side,nc,ord,ll){
  # ll=list(labels,cex,las,srt)
  # inspired by Leo Leopold, with modifications to rotate text from
  # http://menugget.blogspot.com/2014/08/rotated-axis-labels-in-r-plots.html
  
  if(is.null(ll)) return()
  
  if(length(ll$labels) != nc)
    stop("The length of labels of side ", side, " does not match the number of columns of the corrgram.")
  
  # default cex, srt
  ll$labels = ll$labels[ord]
  if(is.null(ll$cex)) ll$cex=1
  if(is.element(side, c(1, 3)) && is.null(ll$srt)) ll$srt=90 # vert
  if(is.element(side, c(2, 4)) && is.null(ll$srt)) ll$srt=0 # horiz
  
  for(i in seq_len(nc)){
    # row/column grid position down/right 
    if(side==1) { # bottom
      par(mfg=c(nc, i))
      # without 'clip', only the first label is added
      clip(0, -2, 0, 1)
      text(x=0.5, y= 0 - 0.05*(1-0), labels=ll$labels[i],
           cex=ll$cex, srt=ll$srt, adj=1, xpd=NA)
    } else if (side==2){ # left
      par(mfg=c(i, 1))
      clip(0, -2, 0, 1)
      text(x=0 - 0.05*(1-0), y=0.5, labels=ll$labels[i],
           cex=ll$cex, srt=ll$srt, adj=1, xpd=NA)
    } else if (side==3) { # top
      par(mfg=c(1, i))
      clip(0, -2, 0, 1)
      text(x=0.5, y= 1 + 0.05*(1-0), labels=ll$labels[i],
           cex=ll$cex, srt=ll$srt, adj=0, xpd=NA)
    } else if (side==4) { # right
      par(mfg=c(i, nc))
      clip(0, -2, 0, 1)
      text(x=1 + .05*(1-0), y=0.5, labels=ll$labels[i],
           cex=ll$cex, srt=ll$srt, adj=0, xpd=NA)
    } else {
      stop("'side' must be 1, 2, 3, or 4")
    }
    
  }
  
  return()
}

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Panel functions

#' @export
panel.pts <- function(x, y, corr=NULL, col.regions, cor.method, ...){

  # For correlation matrix, do nothing
  if(!is.null(corr)) return()

  plot.xy(xy.coords(x, y), type="p", ...)
  box(col="lightgray")
}

#' @export
panel.pie <- function(x, y, corr=NULL, col.regions, cor.method, ...){

  # coordinates of box
  usr <- par()$usr  # par is in graphics package
  minx <- usr[1]; maxx <- usr[2]
  miny <- usr[3]; maxy <- usr[4]
  # multiply the radius by .96 so the circles do not overlap
  rx <- (maxx-minx)/2 * .96
  ry <- (maxy-miny)/2 * .96
  centerx <- (minx+maxx)/2
  centery <- (miny+maxy)/2

  # if corr not given, try to calculate it
  if(is.null(corr)) {
    if(sum(complete.cases(x,y)) < 2) {
      warning("Need at least 2 complete cases for cor()")
      return()
    } else {
      corr <- cor(x, y, use='pair', method=cor.method)
    }
  }

  # draw circle
  segments <- 180
  angles <- seq.int(0, 2*pi, length.out=segments)
  circ <- cbind(centerx + cos(angles)*rx, centery + sin(angles)*ry)
  lines(circ[,1], circ[,2], col='gray40',...)

  ncol <- 14
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks=seq.int(from=-1, to=1, length.out=ncol+1),
                            include.lowest=TRUE))
  col.pie <- pal[col.ind]

  # overlay a colored polygon
  segments <- round(180*abs(corr),0)
  if(segments>0){ # Watch out for the case with 0 segments
    angles <- seq.int(pi/2, pi/2+(2*pi* -corr), length.out=segments)
    circ <- cbind(centerx + cos(angles)*rx, centery + sin(angles)*ry)
    circ <- rbind(circ, c(centerx, centery), circ[1,])
    polygon(circ[,1], circ[,2], col=col.pie,border="gray40")
  }

}

#' @export
panel.shade <- function(x, y, corr=NULL, col.regions, cor.method, ...){

  # If corr not given, try to calculate it
  if(is.null(corr)) {
    if(sum(complete.cases(x,y)) < 2) {
      warning("Need at least 2 complete cases for cor()")
      return()
    } else {
      corr <- cor(x, y, use='pair', method=cor.method)
    }
  }

  ncol <- 14
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks=seq.int(from=-1, to=1, length.out=ncol+1),
                            include.lowest=TRUE))
  usr <- par("usr")
  # Solid fill
  rect(usr[1], usr[3], usr[2], usr[4], col=pal[col.ind], border=NA)
  # Add diagonal lines

  if(!is.na(corr)) {
    rect(usr[1], usr[3], usr[2], usr[4], density=5,
         angle=ifelse(corr>0, 45, 135), col="white")
  }
  # Boounding box needs to plot on top of the shading, so do it last.
  box(col='lightgray')
}

#' @export
panel.ellipse <- function(x,y, corr=NULL, col.regions, cor.method, ...){

  # For correlation matrix, do nothing
  if(!is.null(corr)) return()

  # If too few points, do nothing
  if(sum(complete.cases(x,y)) < 2) {
    warning("Need at least 2 complete cases to draw ellipse.")
    return()
  }

  # Draw an ellipse
  dfn <- 2
  dfd <- length(x)-1
  shape <- var(cbind(x,y),na.rm=TRUE)
  keep <- (!is.na(x) & !is.na(y))
  center <- c(mean(x[keep]),mean(y[keep]))
  radius <- sqrt(dfn*qf(.68,dfn,dfd))
  segments <- 180
  angles <- seq.int(0,2*pi,length.out=segments)
  unit.circle <- cbind(cos(angles),sin(angles))
  ellipse.pts <- t(center+radius*t(unit.circle%*%chol(shape)))
  ellx <- ellipse.pts[,1]
  elly <- ellipse.pts[,2]
  # Truncate ellipse at min/max or at bounding box
  usr <- par()$usr
  minx <- usr[1] #min(x, na.rm=TRUE)
  maxx <- usr[2] #max(x, na.rm=TRUE)
  miny <- usr[3] #min(y, na.rm=TRUE)
  maxy <- usr[4] #max(y, na.rm=TRUE)
  ellx <- ifelse(ellx < minx, minx, ellx)
  ellx <- ifelse(ellx > maxx, maxx, ellx)
  elly <- ifelse(elly < miny, miny, elly)
  elly <- ifelse(elly > maxy, maxy, elly)
  lines(ellx, elly, col='gray40',...)

  # Filled ellipse
  # polygon(ellx, elly, col="blue", ...)

  # Add a lowess line through the ellipse.  Use 'ok' to remove NAs
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    lines(stats::lowess(x[ok], y[ok], f = 2/3, iter = 3),
          col = "red", ...)
}

#' @export
panel.bar <- function(x, y, corr=NULL, col.regions, cor.method, ...){
  # Use 'bars' as in Friendly, figure 1

  usr <- par()$usr
  minx <- usr[1]; maxx <- usr[2]
  miny <- usr[3];  maxy <- usr[4]

  if (is.null(corr)) {
    if(sum(complete.cases(x,y)) < 2) {
      warning("Need at least 2 complete cases for cor()")
      return()
    } else {
      corr <- cor(x, y, use = "pair", method=cor.method)
    }
  }
  
  ncol <- 14
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks = seq.int(from = -1, to = 1,
                                    length.out = ncol + 1), include.lowest = TRUE))
  col.bar <- pal[col.ind]
  if(corr < 0) {
    # Draw up from bottom
    maxy <- miny + (maxy-miny) *  abs(corr)
    rect(minx, miny, maxx, maxy, col = pal[col.ind],
         border = "lightgray")
  } else if (corr > 0){
    # Draw down from top
    miny <- maxy - (maxy-miny)*corr
    rect(minx, miny, maxx, maxy, col = pal[col.ind],
         border = "lightgray")
  }
  
}

#' @export
panel.cor <- function(x, y, corr=NULL, col.regions, cor.method, digits=2, cex.cor, ...){
  # Correlation values only, colored

  # If corr not given, try to calculate it
  if(is.null(corr)) {
    if(sum(complete.cases(x,y)) < 2) {
      warning("Need at least 2 complete cases for cor()")
      return()
    } else {
      corr <- cor(x, y, use='pair', method=cor.method)
    }
  }

  auto <- missing(cex.cor)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  ncol <- 14
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks=seq.int(from=-1, to=1, length.out=ncol+1),
                            include.lowest=TRUE))
  
  corr <- formatC(corr, digits=digits, format='f')
  if(auto) cex.cor <- 0.7/strwidth(corr)
  text(0.5, 0.5, corr, cex=cex.cor, col=pal[col.ind])

}

#' @export
panel.conf <- function(x, y, corr=NULL, col.regions, cor.method, digits=2, cex.cor, ...){

  auto <- missing(cex.cor)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))

  # ncol <- 14
  # pal <- col.regions(ncol)
  
  # For correlation matrix, only show the correlation
  if(!is.null(corr)) {
    est <- corr
    #col.ind <- as.numeric(cut(est, breaks=seq(from=-1, to=1, length=ncol+1),
    #                          include.lowest=TRUE))
    est <- formatC(est, digits=digits, format='f')
    if(auto) cex.cor <- 0.7/strwidth(est)
    text(0.5, 0.6, est, cex=cex.cor) #, col=pal[col.ind])

  } else { # Calculate correlation and confidence interval
    if(sum(complete.cases(x,y)) < 4) {
      warning("Need at least 4 complete cases for cor.test()")
    } else {
      results <- cor.test(x, y, alternative = "two.sided")
      
      # First, the estimate
      est <- results$estimate
      #col.ind <- as.numeric(cut(est, breaks=seq(from=-1, to=1, length=ncol+1),
      #                          include.lowest=TRUE))
      est <- formatC(est, digits=digits, format='f')
      if(auto) cex.cor <- 0.7/strwidth(est)
      text(0.5, 0.6, est, cex=cex.cor) #, col=pal[col.ind])
      
      ci <- results$conf.int
      ci <- formatC(ci, digits=2, format='f')
      ci <- paste0("(",ci[1],",",ci[2],")")
      if(auto) cex.cor <- 0.8/strwidth(ci)
      text(0.5, 0.3, ci, cex=cex.cor) # , col=pal[col.ind])
    }
    
  }
}

#' @export
panel.txt <- function(x=0.5, y=0.5, txt, cex, font, srt){
  text(x, y, txt, cex=cex, font=font, srt=srt)
}

#' @export
panel.density <- function(x, corr=NULL, ...){
  # For correlation matrix, do nothing
  if(!is.null(corr)) return()

  dd = density(x, na.rm=TRUE)
  xr=range(dd$x)
  yr=range(dd$y)
  par(usr = c(min(xr), max(xr), min(yr), max(yr)*1.1))
  plot.xy(xy.coords(dd$x, dd$y), type="l", col="black", ...)
  box(col="lightgray")
}

#' @export
panel.minmax <- function(x, corr=NULL, ...){
  # For correlation matrix, do nothing
  if(!is.null(corr)) return()
  # Put the minimum in the lower-left corner and the
  # maximum in the upper-right corner
  minx <- round(min(x, na.rm=TRUE),2)
  maxx <- round(max(x, na.rm=TRUE),2)
  text(minx, minx, minx, cex=1, adj=c(0,0))
  text(maxx, maxx, maxx, cex=1, adj=c(1,1))
}
