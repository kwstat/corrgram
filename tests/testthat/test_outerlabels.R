# test_outerlabels.R
# Time-stamp: <05 May 2017 16:29:11 c:/x/rpack/corrgram/tests/testthat/test_outerlabels.R>

context("test_outerlabels.R")

require(corrgram)

# short syntax for outer labels
corrgram(state.x77, outer.labels=list(bottom=TRUE, right=TRUE))

# use default labels in outer margin
corrgram(state.x77, outer.labels=list(bottom=TRUE, right=list(srt=25)))

labs=c("Population", "Income", "Illiteracy", "Life Exp", "Murder", "HS Grad", "Frost", "Area")

# outer.labels not given
corrgram(state.x77)

# outer labels, one side at a time
corrgram(state.x77, outer.labels=list(bottom=list(labels=labs)))
corrgram(state.x77, outer.labels=list(left=list(labels=labs)))
corrgram(state.x77, outer.labels=list(top=list(labels=labs)))
corrgram(state.x77, outer.labels=list(right=list(labels=labs)))

# outer labels with no diagonal labels
corrgram(state.x77, text.panel=NULL, 
         outer.labels=list(bottom=list(labels=labs)))

# outer.labels, all 4 sides at once
corrgram(state.x77,
         outer.labels=list(bottom=list(labels=labs),
                           left=list(labels=labs),
                           top=list(labels=labs),
                           right=list(labels=labs)))

# outer.labels, all 4 sides at once, re-ordered
corrgram(state.x77, order=TRUE,
         outer.labels=list(bottom=list(labels=labs),
                           left=list(labels=labs),
                           top=list(labels=labs),
                           right=list(labels=labs)))

# outer labels, srt
corrgram(state.x77,
         outer.labels=list(bottom=list(labels=labs,srt=60),
                           left=list(labels=labs,srt=30),
                           top=list(labels=labs,srt=90),
                           right=list(labels=labs,srt=0)))

# outer labels, cex
corrgram(state.x77, outer.labels=list(bottom=list(labels=labs,cex=0.5)))
corrgram(state.x77, outer.labels=list(left=list(labels=labs,cex=1)))
corrgram(state.x77, outer.labels=list(top=list(labels=labs,cex=1.5)))
corrgram(state.x77, outer.labels=list(right=list(labels=labs,cex=2)))

# outer labels, all options, larger margins, xlab, ylab
corrgram(state.x77, oma=c(7, 7, 2, 2), main="state.x77",
         outer.labels=list(bottom=list(labels=labs,cex=1.5,srt=60),
                           left=list(labels=labs,cex=1.5,srt=30)))
mtext("Bottom", side=1, cex=2, line = -1.5, outer=TRUE, xpd=NA)
mtext("Left", side=2, cex=2, line = -1.5, outer=TRUE, xpd=NA)

test_that("outer labels are wrong length", {
  expect_error(corrgram(state.x77, outer.labels=list(bottom=list(labels=labs[-1]))))
  expect_error(corrgram(state.x77, outer.labels=list(left=list(labels=labs[-1]))))
  expect_error(corrgram(state.x77, outer.labels=list(top=list(labels=labs[-1]))))
  expect_error(corrgram(state.x77, outer.labels=list(right=list(labels=labs[-1]))))
})
