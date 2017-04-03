# misc.r
# Time-stamp: <03 Apr 2017 11:51:31 c:/x/rpack/corrgram/tests/testthat/test_raw.R>

require(corrgram)

# type
corrgram(vote)
corrgram(vote, type='corr')
corrgram(vote, type='data') # warn user
corrgram(vote, lower.panel=panel.conf)

# cor.method
corrgram(auto)
corrgram(auto, cor.method="pearson")
corrgram(auto, cor.method="spearman") # Slight change in colors

# ignore non-numeric columns
corrgram(iris)

# labels
corrgram(mtcars[2:6],
         labels=c('Axle ratio','Weight','Displacement','Cylinders','Horsepower'))

# label.srt, diagonal labels rotated 45 degrees
corrgram(auto, label.srt=-45)
  
# label.cex, label.pos
corrgram(auto, label.srt=45, label.pos=c(.75,.75), cex.labels=2.5, upper=NULL)


# order argument
corrgram(mtcars)
corrgram(mtcars, order=NULL)
corrgram(mtcars, order=FALSE)
corrgram(mtcars, order=TRUE)
corrgram(mtcars, order="GW")
corrgram(mtcars, order="HC")
corrgram(mtcars, order="PC")
corrgram(mtcars, order="OLO")
corrgram(mtcars, order="PC", abs=TRUE)
corrgram(mtcars, order="OLO", abs=TRUE)

# make sure 'labels' works correctly with 'order'
myLabels = names(mtcars)
myLabels[myLabels == "hp"] = "horse\npower"
corrgram(mtcars, lower.panel = panel.conf, labels = myLabels)
cmat1 <- corrgram(mtcars, lower.panel = panel.conf, labels = myLabels, order = TRUE)


# diagonal direction
corrgram(auto, order=TRUE, dir="right")
corrgram(auto, order=TRUE, dir="/")


# off-diagonal panels
corrgram(auto, panel=panel.bar)
corrgram(auto, panel=panel.conf)
corrgram(auto, panel=panel.cor)
corrgram(auto, panel=panel.ellipse) # note: latticeExtra also has panel.ellipse
corrgram(auto, panel=panel.pie)
corrgram(auto, panel=panel.pts)
corrgram(auto, panel=panel.shade)
# text/diag panels
corrgram(auto, text.panel=NULL, diag.panel=panel.density)
corrgram(auto, text.panel=panel.txt, diag.panel=panel.minmax)
  

# col.regions with all panels
col.earth <- colorRampPalette(c("darkgoldenrod4", "burlywood1", "darkkhaki", "darkgreen"))
# off-diagonal panels
corrgram(auto, panel=panel.bar,col.regions=col.earth)
corrgram(auto, panel=panel.conf,col.regions=col.earth)
corrgram(auto, panel=panel.cor,col.regions=col.earth)
corrgram(auto, panel=panel.ellipse,col.regions=col.earth)
corrgram(auto, panel=panel.pie,col.regions=col.earth)
corrgram(auto, panel=panel.pts,col.regions=col.earth)
corrgram(auto, panel=panel.shade,col.regions=col.earth)
# text/diag panels
corrgram(auto, text.panel=NULL, diag.panel=panel.density,col.regions=col.earth)
corrgram(auto, text.panel=panel.txt, diag.panel=panel.minmax,col.regions=col.earth)

corrgram(mtcars, order=TRUE, lower.panel=panel.shade, upper.panel=panel.pie,
         main="A Corrgram of a Different Color",
         col.regions=col.earth)

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if(FALSE){ # No need to test automatically

  # Split long variable names on two lines
  corrgram(mtcars[2:6], order=TRUE, upper.panel=NULL,
           lower.panel=panel.pie,
           text.panel=panel.txt,
           labels=rep('A very long \n variable name',4))

  # Bug with negative correlation
  
  set.seed(123)
  a = seq(1,100)
  b = jitter(seq(1,100), 80)
  cor(a,b) # r about .95
  ab=as.data.frame(cbind(a,b))
  ab$c = -1 * ab$b # flip direction of correlation
  cor(ab$a, ab$c) # r now about -.95
  corrgram(ab, order=NULL, lower.panel=panel.pie, upper.panel=NULL,
           text.panel=panel.txt)
  corrgram(ab)

  # missing value in a correlation matrix.
  vote2 <- vote
  vote2[2:6,2:6] <- NA
  corrgram(vote2)

  # missing combinations cause cor( , use="pair") to be NAs
  dat <- data.frame(E1=c(NA,NA,NA,NA,NA,6,7,8,9,10),
                    E2=c(1,2,3,4,5,NA,NA,NA,NA,NA),
                    E3=c(1,2,3,4,5,6,7,8,9,10)+.1,
                    E4=c(2,1,5,6,8,7,9,4,5,3))
  cor(dat, use="pair")
  corrgram(dat)

  # diagonal labels unclipped.
  # This has a slight quirk...the red box is only drawn the first time.  Calling
  # corrgram a 2nd time doesn't draw the red box.
  require('grid')
  require('gridBase')
  unclipped.txt <- function(x=0.1, y=0.5, txt, cex, font, srt){
    vps <- baseViewports()
    vps$figure$clip <- NA # Hack. Do NOT clip text that falls outside the ploting region
    pushViewport(vps$inner) # Figure region
    #grid.rect(gp=gpar(lwd=3, col="red"))
    pushViewport(vps$figure) # The diagonal box region
    #grid.rect(gp=gpar(lwd=3, col="blue"))
    grid.text(txt, x=0.1, y=y, just='left', gp=gpar(cex=cex))
    popViewport(2)
  }
  ## corrgram(mtcars[2:6], order=FALSE,
  ##          lower.panel=panel.conf)
  # Order of labels need to match original data
  corrgram(mtcars[2:6], order=TRUE,
           labels=c("Cylinders","Displacement","Horsepower","Axle ratio","Weight"),
           cex.labels=2, adj=0,
           upper.panel=NULL, lower.panel=panel.conf,
           diag.panel=NULL, text.panel=unclipped.txt)
  
  # Manually add a legend for coloring points
  panel.colpts <- function(x, y, corr=NULL, col.regions, ...){
    # For correlation matrix, do nothing
    if(!is.null(corr)) return()
    plot.xy(xy.coords(x, y), type="p", ..., col=1:2)
    box(col="lightgray")
  }
  corrgram(auto, lower.panel=panel.conf, upper.panel=panel.colpts)
  require(grid)
  grid.clip()
  pushViewport(viewport(.5, .95, width=stringWidth("Group1"),
                        height=unit(2,"lines"),
                        name="pagenum", gp=gpar(fontsize=8)))
  grid.legend(pch=1:2, labels=c("Group1","Group2"), gp=gpar(col=c('red')))
  popViewport()

  # one pair of variables had no complete observations
  dati <- iris[1:50,]
  dati[seq(from=2, to=50, by=2),1] <- NA
  dati[seq(from=1, to=49, by=2),2] <- NA
  # off-diagonal panels
  corrgram(dati, panel=panel.bar)
  corrgram(dati, panel=panel.conf)
  corrgram(dati, panel=panel.cor)
  corrgram(dati, panel=panel.ellipse)
  corrgram(dati, panel=panel.pie)
  corrgram(dati, panel=panel.pts)
  corrgram(dati, panel=panel.shade)
  # text/diag panels
  corrgram(dati, text.panel=NULL, diag.panel=panel.density)
  corrgram(dati, text.panel=NULL, diag.panel=panel.minmax)
  
} # end if

