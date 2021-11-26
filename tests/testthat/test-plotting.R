test_that("plots execute error-free", {
  data("exrates", package = "stochvol")
  m <- 6
  n <- 50
  y <- 100 * logret(tail(exrates[, seq_len(m)], n + 1))
  res <- fsvsample(y, factors = 2, draws = 100, zeromean = FALSE,
                   thin = 1, quiet = TRUE)
  expect_invisible(comtimeplot(res))
  expect_invisible(voltimeplot(res))
  expect_invisible(corimageplot(res, plottype = "corrplot"))
  expect_invisible(corimageplot(res, plottype = "imageplot"))
  expect_invisible(corimageplot(res, plotCI = "circle"))
  expect_invisible(cortimeplot(res, series = 2))
  expect_invisible(covtimeplot(res, series = 2))
  expect_invisible(facloadpairplot(res))
  expect_invisible(facloadcredplot(res))
  expect_invisible(facloadpointplot(res))
  expect_invisible(logvartimeplot(res))
  expect_invisible(paratraceplot(res))
  expect_invisible(facloadtraceplot(res))
  expect_invisible(facloaddensplot(res))
  expect_invisible(corplot(res))
  expect_invisible(evdiag(res))
})
