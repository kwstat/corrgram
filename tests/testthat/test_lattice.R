
require(corrgram)
require(lattice)

test_that("lattice levelplot panel functions do not err", {

  # vote data is in corrgram
  levelplot(vote, at = do.breaks(c(-1.01, 1.01), 20),
            xlab = NULL, ylab = NULL, colorkey = list(space = "top"),
            scales = list(x = list(rot = 90)), 
            panel = levelplot_panel.ellipse, label = TRUE)
  levelplot(vote, at = do.breaks(c(-1.01, 1.01), 20),
            xlab = NULL, ylab = NULL, colorkey = list(space = "top"),
            scales = list(x = list(rot = 90)), 
            panel = levelplot_panel.pie, label = TRUE)

})

test_that("lattice splom panel functions do not err", {

  pengvars <- c("bill_len", "bill_dep", "flipper_len", "body_mass")
  splom(~penguins[ , pengvars], upper.panel=splom_panel.pie, pscales=0)
  splom(~penguins[ , pengvars]|penguins$species, upper.panel=splom_panel.pie, pscales=0)
  splom(~penguins[ , pengvars], upper.panel=splom_panel.shade, pscales=0)
  splom(~penguins[ , pengvars]|penguins$species, upper.panel=splom_panel.shade, pscales=0)
  splom(~penguins[ , pengvars], upper.panel=splom_panel.ellipse, pscales=0)
  splom(~penguins[ , pengvars]|penguins$species, upper.panel=splom_panel.ellipse, pscales=0)

})