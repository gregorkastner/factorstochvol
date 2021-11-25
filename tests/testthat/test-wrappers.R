context("Samplers")

test_that("fsvsample works", {
  expect_warning(fsvsample(y, draws = draws, burnin = burnin, quiet = TRUE), NA) %>%
    expect_is("fsvdraws")
  for (th in thin_values) {
    for (pflt in priorfacloadtype_values) {
      for (hsk in heteroskedastic_values) {
        for (fs in factors_values) {
          expect_warning(fsvsample(y, draws = draws, burnin = burnin, quiet = TRUE,
                                   factors = fs, thin = th, priorfacloadtype = pflt,
                                   restrict = "none",
                                   heteroskedastic = hsk,
                                   priorhomoskedastic = if (!isTRUE(hsk)) priorhomoskedastic[seq_len(NCOL(y)), ] else NA,
                                   interweaving = if (!isTRUE(hsk)) 0 else 4,
                                   runningstore = if (fs == 0) 1 else 6), NA) %>%
            expect_is("fsvdraws")
          if (fs > 1) {
            for (rst in restrict_values) {
              expect_warning(fsvsample(y, draws = draws, burnin = burnin, quiet = TRUE,
                                       factors = fs, thin = th, priorfacloadtype = pflt,
                                       heteroskedastic = hsk,
                                       priorhomoskedastic = if (!isTRUE(hsk)) priorhomoskedastic[seq_len(NCOL(y)), ] else NA,
                                       interweaving = if (!isTRUE(hsk)) 0 else 4,
                                       restrict = rst), NA) %>%
                expect_is("fsvdraws")
            }
          }
        }
      }
    }
  }

  # bug report from Kevin Kurt
  sim_res <- fsvsim(n = 10, series = 3, factors = 2)
  expect_warning(fsvsample(y = sim_res$y, factors = 2, draws = 10, quiet = TRUE,
                           priorsigmafac = rep(1, 2)), NA)
})

