data("exrates", package = "stochvol")

y <- logret(exrates[1:21, 1:7], demean = TRUE)
draws <- 20
burnin <- 10

thin_values <- c(1, 3)
factors_values <- c(0, 1, 3)
restrict_mat <- matrix(FALSE, nrow = NCOL(y), ncol = max(factors_values))
restrict_mat[1, 1] <- TRUE  # restrict the upper left element to zero
restrict_values <- list("upper", "auto", restrict_mat)
priorfacloadtype_values <- c("normal", "rowwiseng", "colwiseng") #, "dl" is not supported anymore
priorhomoskedastic <- matrix(c(1.1, 1.1), nrow = NCOL(y),
                             ncol = 2, byrow = TRUE)
heteroskedastic_values <- list(TRUE, c(FALSE, FALSE))

res <- list(fsvsample(y, draws = draws, burnin = burnin, quiet = TRUE))
for (th in thin_values) {
  for (pflt in priorfacloadtype_values) {
    for (hsk in heteroskedastic_values) {
      for (fs in factors_values) {
        res <- c(res, list(fsvsample(y, draws = draws, burnin = burnin, quiet = TRUE,
                                     factors = fs, thin = th, priorfacloadtype = pflt,
                                     restrict = "none",
                                     heteroskedastic = hsk,
                                     priorhomoskedastic = if (!isTRUE(hsk)) priorhomoskedastic else NA,
                                     interweaving = if (!isTRUE(hsk)) 0 else 4,
                                     runningstore = if (fs == 0) 1 else 6)))
        if (fs > 1) {
          for (rst in restrict_values) {
            res <- c(res, list(fsvsample(y, draws = draws, burnin = burnin, quiet = TRUE,
                                         factors = fs, thin = th, priorfacloadtype = pflt,
                                         heteroskedastic = hsk,
                                         priorhomoskedastic = if (!isTRUE(hsk)) priorhomoskedastic else NA,
                                         interweaving = if (!isTRUE(hsk)) 0 else 4,
                                         restrict = rst)))
          }
        }
      }
    }
  }
}

