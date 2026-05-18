# corrgram_levelplot.R

#' Panel functions for lattice corrgrams via levelplot
#'
#' These functions provide custom panel methods for \code{lattice::levelplot()}.
#'
#' @param x Numeric coordinates from levelplot.
#' @param y Numeric coordinates from levelplot.
#' @param z Correlation values from levelplot.
#' @param subscripts Subscripts for lattice panel. (not used)
#' @param at Breaks for color levels.
#' @param level Confidence level for ellipse (default 0.9).
#' @param label Logical; if TRUE, show correlation values as text.
#' @param ... Additional arguments passed to panel functions.
#'
#' @details
#' 
#' \code{levelplot_panel.ellipse} Draws ellipses representing correlation coefficients
#' in the upper triangle of a matrix. Optionally adds numeric labels in the 
#' lower triangle.
#' 
#' \code{levelplot_panel.pie} Draws pie glyphs representing correlation coefficients,
#' omitting the diagonal.
#' 
#' @rdname levelplot_panel.ellipse
#' @importFrom ellipse ellipse
#' @importFrom lattice levelplot level.colors panel.polygon panel.text
#' @examples
#' library(lattice)
#' levelplot(vote, at = do.breaks(c(-1.01, 1.01), 20),
#'             xlab = NULL, ylab = NULL, colorkey = list(space = "top"),
#'             scales = list(x = list(rot = 90)), 
#'             panel = levelplot_panel.ellipse, label = TRUE)
#' levelplot(vote, at = do.breaks(c(-1.01, 1.01), 20),
#'             xlab = NULL, ylab = NULL, colorkey = list(space = "top"),
#'             scales = list(x = list(rot = 90)), 
#'             panel = levelplot_panel.pie, label = TRUE)
#' @export
levelplot_panel.ellipse <- function(
  x, y, z, subscripts, 
  at, level = 0.9, label = FALSE, ...) {
  x <- as.numeric(x)[subscripts]
  y <- as.numeric(y)[subscripts]
  z <- as.numeric(z)[subscripts]
  zcol <- level.colors(z, at = at, ...)

  # Draw the glyph for each correlation value
  for (i in seq(along = z)) {
    # Only dry glyphs for the upper triangle
    if(x[i] < y[i]) {
      xy.ell <- ellipse::ellipse(z[i], level = level, npoints = 50,
                   scale = c(.2, .2), centre = c(x[i], y[i]))
      panel.polygon(xy.ell, col = zcol[i], border = zcol[i], ...)
    }
  }

  if (label) {
    ixskip <- x <= y # only label the lower triangle
    panel.text(x = x[!ixskip], y = y[!ixskip],
               lab = 100 * round(z[!ixskip], 2), cex = 0.8,
               col = ifelse(z[!ixskip] < 0, "red", "navy"))
  }
}

# --------------------------------------------------------------------

#' @param scale Numeric; scaling factor for pie size (default 0.8).
#' 
#' @rdname levelplot_panel.ellipse
#' 
#' @export
levelplot_panel.pie <- function(x, y, z, subscripts, at = pretty(z), scale = 0.9, ...) {
  # x,y are coordinates of the matrix, z is the correlation value
  x <- as.numeric(x)[subscripts]
  y <- as.numeric(y)[subscripts]
  z <- as.numeric(z)[subscripts]

  zcol <- level.colors(z, at = at, ...)

  for (i in seq(along = z)) {
  if(x[i] == y[i]) next # skip the diagonal
    # The full circle
  grid::grid.circle(x = x[i], y = y[i], r = .5 * scale, default.units = "native")

  # pie-shaped piece for each correlation value
    numseg <- round(60*abs(z[i]),0) # number of segments for each pie piece
    if(numseg>0){ # Watch out for the case with 0 segments
      angles <- seq(pi/2, pi/2+(2*pi* -z[i]), length=numseg)
      xy.pie <- cbind(x[i] + cos(angles)*0.5*scale,
                      y[i] + sin(angles)*0.5*scale)
      xy.pie <- rbind(xy.pie, c(x[i], y[i]), xy.pie[1,])
      panel.polygon( xy.pie, col=zcol[i])
    }
  }
}

