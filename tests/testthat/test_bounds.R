# test_bounds.R
# Time-stamp: <23 Apr 2019 14:49:16 c:/x/rpack/corrgram/tests/testthat/test_bounds.R>

dat <- data.frame(
  a=c(4.427371, 4.652635, 5.117894, 4.648367),
  b=c(3.065796, 2.652594, 2.702682, 2.664237),
  c=c(1.820023, 1.248203, 1.583516, 0.943176),
  d=c(1.0981531, 1.3303935, 1.4942022, 0.9183384),
  e=c(2.445383, 2.230344, 2.796065, 2.108716),
  f=c(-3.058478, -3.287118, -3.063749, -2.809525),
  g=c(3.512448, 3.196442, 3.450073, 3.190861),
  h=c(2.16616, 2.601121, 2.175354, 2.87997),
  i=c(3.223686, 2.84551, 2.959626, 2.845965),
  j=c(2.298442, 1.562345, 1.897684, 1.541668),
  k=c(-1.481442, -2.59257, -1.751966, -1.667766),
  l=c(-3.321619, -3.422294, -3.501768, -3.175093),
  m=c(1.898439, 1.331196, 2.141605, 1.928749),
  n=c(2.0098, 1.713265, 2.110812, 2.506027),
  o=c(3.180492, 3.101148, 3.488329, 3.329727),
  p=c(2.397239, 2.091732, 2.521272, 2.510106),
  q=c(1.3390112, 0.8297595, 1.1517592, 1.0869653)
)

cmat <- cor(dat, method='spearman')

# In R versions before 2015-12-23, the cmat matrix had a correlation < 1
# exactly -1 - .Machine$double.eps
# Update: cor() now returns values in [-1,1].  See
# https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=16638
# print(min(cmat), digits=17) # -1.0000000000000002
# is.matrix(cmat) # TRUE
# isSymmetric(cmat) # TRUE
# min(cmat) > -1 - .Machine$double.eps # FALSE
# min(cmat) >= -1 - .Machine$double.eps # TRUE
# max(cmat) < 1 + .Machine$double.eps # TRUE

require(corrgram)
corrgram(cmat, lower.panel=panel.shade, upper=NULL,
         main="Correlation matrix")

test_that("Calculated spearman correlations are correct", {
  cc <- corrgram(dat, cor.method='spearman', upper=NULL)
  expect_equal(cmat, cc)
})


# This correlation matrix was produced by cov2cor, but has values > 1
# Not sure if cov2cor could still do this after bug fixed above
cmat2 <-
  structure(c(1, 0.952878953343098, 0.0373396513876667,
              0.668048524607397, 0.952888464873119, 0.952876206931305,
              0.556626045335446, 0.952878953343098, 1,
              0.0391818179842369, 0.701082452069751, 1.00000600150022,
              0.999994904466531, 0.584159785791795, 0.0373396513876667,
              0.0391818179842369, 1, 0.0274687144977842,
              0.0391956602070947, 0.0391788599536022,
              0.0228895131694637, 0.668048524607397, 0.701082452069751,
              0.0274687144977842, 1, 0.701101131756589,
              0.701081039882255, 0.409550411443535, 0.952888464873119,
              1.00000600150022, 0.0391956602070947, 0.701101131756589,
              1, 0.999980461387252, 0.584144381934955,
              0.952876206931305, 0.999994904466531, 0.0391788599536022,
              0.701081039882255, 0.999980461387252, 1,
              0.584163309675686, 0.556626045335446, 0.584159785791795,
              0.0228895131694637, 0.409550411443535, 0.584144381934955,
              0.584163309675686, 1), .Dim = c(7L, 7L),
            .Dimnames = list(NULL, NULL))
test_that("Correlation values outside [-1,1] are flagged.", {
  expect_warning(corrgram(cmat2))
  expect_warning(corrgram(-cmat2))
})
