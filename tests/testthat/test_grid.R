require(corrgram)

vars6 <- setdiff(colnames(auto), c("Model", "Origin"))

test_that("corrgram2 dataframe panel functions do not err", {

  corrgram2(auto[, vars6], panel=grid_panel.fill, title="auto data")
  corrgram2(auto[, vars6], panel=grid_panel.shade)
  corrgram2(auto[, vars6], panel=grid_panel.pie)
  corrgram2(auto[, vars6], panel=grid_panel.conf)
  corrgram2(auto[, vars6], panel=grid_panel.ellipse)

})

test_that("corrgram2 matrix panel functions do not err", {

  corrgram2(vote)
  corrgram2(vote, abs=TRUE, labels=1:12, label.srt=30, cex.labels=1.5)
  corrgram2(vote, abs=TRUE, order=TRUE, legend=TRUE)
  corrgram2(vote, panel=grid_panel.conf)
  corrgram2(vote, panel=grid_panel.ellipse) # does nothing
  corrgram2(vote, lower.panel=grid_panel.fill)
  corrgram2(vote, panel=grid_panel.pie, order=TRUE)
  corrgram2(vote, upper.panel=grid_panel.shade)
  corrgram2(vote, lower.panel=grid_panel.fill, upper.panel=grid_panel.shade)

})

test_that("corrgram2 vote with OLO order does not err", {
  skip_if_not_installed("seriation")
  corrgram2(vote, order="OLO", lower.panel=grid_panel.pie, 
  title="vote data with OLO order", cex.labels="fit", label.srt=45)
})
