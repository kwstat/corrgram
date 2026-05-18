# corrgram_splom.R

#' Panel functions for lattice corrgrams via splom
#'
#' These functions provide custom panel methods for \code{lattice::splom()}.
#'
#' @param x Numeric coordinates from levelplot.
#' @param y Numeric coordinates from levelplot.
#' @param z Correlation values from levelplot.
#' @param subscripts Subscripts for lattice panel. (not used)
#' @param at Breaks for color levels.
#' @param cor.method Correlation method (default "pearson").
#' @param ... Additional arguments passed to panel functions.
#'
#' @details
#' 
#' \code{splom_panel.pie} Draws pie glyphs representing correlation coefficients,
#' omitting the diagonal.
#' 
#' \code{splom_panel.shade} Draws shaded rectangles with hash lines representing correlation
#' coefficients.
#' 
#' \code{splom_panel.ellipse} Draws ellipses representing correlation coefficients.
#' The position of the ellipse is determined by the position of the data,
#' and the shape of the ellipse is determined by the correlation.
#' 
#' @rdname splom_panel.pie
#' @importFrom ellipse ellipse
#' @importFrom lattice level.colors panel.polygon panel.text
#' @examples
#' library(lattice)
#' pengvars <- c("bill_len", "bill_dep", "flipper_len", "body_mass")
#' splom(~penguins[ , pengvars], upper.panel=splom_panel.pie, pscales=0)
#' splom(~penguins[ , pengvars]|penguins$species, upper.panel=splom_panel.pie, pscales=0)
#' splom(~penguins[ , pengvars], upper.panel=splom_panel.shade, pscales=0)
#' splom(~penguins[ , pengvars]|penguins$species, upper.panel=splom_panel.shade, pscales=0)
#' splom(~penguins[ , pengvars], upper.panel=splom_panel.ellipse, pscales=0)
#' splom(~penguins[ , pengvars]|penguins$species, upper.panel=splom_panel.ellipse, pscales=0)
#' @importFrom lattice current.panel.limits panel.xyplot panel.polygon panel.rect splom
#' @export
splom_panel.pie <- function(x, y, z, subscripts, at = pretty(z), 
                      cor.method="pearson", ...) {
  # Each glyph is in a single panel
  # x and y are for only this panel...do not subset by subscripts

  # Corner coordinates of panel box, which is NOT square
  cpl <- lattice::current.panel.limits()
  xmin <- cpl$xlim[1]
  xmax <- cpl$xlim[2]
  ymin <- cpl$ylim[1]
  ymax <- cpl$ylim[2]

  # Multiply the radius by .97 so the circles do not overlap
  rx <- (xmax-xmin)/2 * .97
  ry <- (ymax-ymin)/2 * .97
  centerx <- (xmin+xmax)/2
  centery <- (ymin+ymax)/2

  # Don't think z is passed in as an argument
  if(missing(z)) {
    z <- cor(x,y, use="pairwise.complete.obs", method=cor.method)
  } else {
    z <- as.numeric(z)
  }
  zcol <- level.colors(z, at = at, ...)
  
  # Draw circle
  segments <- 60
  angles <- seq(0,2*pi,length=segments)
  circ <- cbind(centerx + cos(angles)*rx, centery + sin(angles)*ry)
  panel.xyplot(circ[,1], circ[,2], type = "l",
               lty = 1, col="gray30", ...)

  # Overlay a colored pacman polygon
  n_bins <- 14
  pal <- colorRampPalette(c("red","salmon","white","royalblue","navy"))(n_bins)
  col.ind <- as.numeric(cut(z, breaks=seq(from=-1, to=1, length=n_bins+1),
                            include.lowest=TRUE))
  col.pie <- pal[col.ind]
  segments <- round(60*abs(z),0)
  if(segments>0){ # Watch out for the case with 0 segments
    angles <- seq(pi/2, pi/2+(2*pi* -z), length=segments)
    circ <- cbind(centerx + cos(angles)*rx, centery + sin(angles)*ry)
    circ <- rbind(circ, c(centerx, centery), circ[1,])
    panel.polygon(circ[,1], circ[,2], col=col.pie)
  }
}
#splom(~penguins[ , pengvars], upper.panel=splom_panel.pie, pscales=0)
#splom(~penguins[ , pengvars]|penguins$species, upper.panel=splom_panel.pie, pscales=0)

# -----------------------------------------------------------------------------

#' @rdname splom_panel.pie
#' @param col.regions Color palette for shading (default NULL uses internal
#' red white blue palette).
#' @export
splom_panel.shade <- function(x, y, z, subscripts, at = pretty(z),
  col.regions = NULL,
  cor.method = "pearson", ...) {
  # x and y are already subsetted for this panel
  if (length(x) < 2 || length(y) < 2 || 
    all(is.na(x)) || all(is.na(y))) return()

  # Compute correlation
  # Don't think z is passed in as an argument
  if(missing(z)) {
    z <- cor(x,y, use="pairwise.complete.obs", method=cor.method)
  } else {
    z <- as.numeric(z)
  }
  zcol <- level.colors(z, at = at, ...)


  # Color palette
  if (is.null(col.regions)) {
    pal <- colorRampPalette(c("red","salmon","white","royalblue","navy"))(14)
  } else {
    pal <- col.regions(14)
  }
  col.ind <- as.numeric(cut(z, breaks = seq(-1, 1, length.out = 15), include.lowest = TRUE))
  col.shade <- pal[col.ind]

  # Panel limits
  cpl <- lattice::current.panel.limits()
  xmin <- cpl$xlim[1]
  xmax <- cpl$xlim[2]
  ymin <- cpl$ylim[1]
  ymax <- cpl$ylim[2]

  # Solid fill
  panel.rect(xmin, ymin, xmax, ymax, col = col.shade)

  # Hash lines using [0,1] coordinates for the panel
  if(z >0) {
    # 5 lines with slope 1
    x0 = c(0, 0, 0, .33, .67)
    x1 = c(.33, .67, 1, 1, 1)
    y0 = c(.67, .33, 0, 0, 0)
    y1 = c(1, 1, 1, .67, .33)
  } else {
    # slope -1
    x0 = c(0, 0, 0, .33, .67)
    x1 = c(.33, .67, 1, 1, 1)
    y0 = c(.33, .67, 1, 1, 1)
    y1 = c(0, 0, 0, .33, .67)
  }
  
  grid::grid.segments(
      x0 = x0, y0 = y0,
      x1 = x1, y1 = y1,
      default.units = "npc",
      gp = grid::gpar(col = "white", lwd = 1, lty = 1, ...)
    )

}

#splom(~penguins[ , pengvars], upper.panel=splom_panel.shade, pscales=0)
#splom(~penguins[ , pengvars]|penguins$species, upper.panel=splom_panel.shade, pscales=0)

# ----------------------------------------------------------------------------

#' @rdname splom_panel.pie
#' @importFrom lattice llines panel.xyplot panel.rect
#' @export
splom_panel.ellipse <- function(x,y, col.regions, cor.method="pearson", ...){

  # For correlation matrix, do nothing
  #if(!is.null(corr)) return()

  # If too few points in the panel, do nothing
  if(sum(complete.cases(x,y)) < 2) {
    warning("Need at least 2 complete cases to draw ellipse.")
    return()
  }

  # Draw an unfilled ellipse
  dfn <- 2
  dfd <- length(x)-1
  shape <- var(cbind(x,y),na.rm=TRUE)
  keep <- (!is.na(x) & !is.na(y))
  center <- c(mean(x[keep]),mean(y[keep]))
  radius <- sqrt(dfn*qf(.68,dfn,dfd))
  segments <- 180
  angles <- seq(from=0, to=2*pi, length.out=segments)
  unit.circle <- cbind(cos(angles),sin(angles))
  ellipse.pts <- t(center+radius*t(unit.circle%*%chol(shape)))
  ellx <- ellipse.pts[,1]
  elly <- ellipse.pts[,2]

  # Truncate ellipse at min/max or at bounding box
  cpl <- lattice::current.panel.limits()
  xmin <- cpl$xlim[1]
  xmax <- cpl$xlim[2]
  ymin <- cpl$ylim[1]
  ymax <- cpl$ylim[2]
  ellx <- ifelse(ellx < xmin, xmin, ellx)
  ellx <- ifelse(ellx > xmax, xmax, ellx)
  elly <- ifelse(elly < ymin, ymin, elly)
  elly <- ifelse(elly > ymax, ymax, elly)

  llines(ellx, elly, col='gray40')

}
#splom(~penguins[ , pengvars], upper.panel=splom_panel.ellipse, pscales=0)
#splom(~penguins[ , pengvars]|penguins$species, upper.panel=splom_panel.ellipse, pscales=0)

