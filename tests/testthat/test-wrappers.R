context("Samplers")

test_that("fsvsample works", {
  expect_warning(fsvsample(y, draws = draws, burnin = burnin, quiet = TRUE), NA) %>%
    expect_is("fsvdraws")
  for (th in thin_values) {
    for (pflt in priorfacloadtype_values) {
      for (fs in factors_values) {
        expect_warning(fsvsample(y, draws = draws, burnin = burnin, quiet = TRUE,
                                     factors = fs, thin = th, priorfacloadtype = pflt,
                                     restrict = "none",
                                     runningstore = if (fs == 0) 1 else 6), NA) %>%
          expect_is("fsvdraws")
        if (fs > 1) {
          for (rst in restrict_values) {
            expect_warning(fsvsample(y, draws = draws, burnin = burnin, quiet = TRUE,
                                         factors = fs, thin = th, priorfacloadtype = pflt,
                                         restrict = rst), NA) %>%
              expect_is("fsvdraws")
          }
        }
      }
    }
  }
})
