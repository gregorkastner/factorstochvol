context("Iterative-use")

test_that("Drawing from the posterior in one as well as in multiple calls to fsvsample()", {

factors <- 2
iterations <- 2

sim <- fsvsim(factors = factors, series = 3)

# Do all draws in one call to fsvsample
cat("Drawing in one call to fsvsample...\n")
set.seed(1)
tim <- system.time({
  res <- fsvsample(sim$y, factors = factors, draws = iterations, burnin = 0,
		   keeptime = "all", quiet = TRUE, signident = FALSE)
})
print(tim)


# Do the draws in 'iterations' calls to fsvsample

cat("Drawing in many calls to fsvsample...\n")
set.seed(1)
tim2 <- system.time({
  res2 <- fsvsample(sim$y, factors = factors, draws = 1, burnin = 0,
		    keeptime = "all", quiet = TRUE, signident = FALSE)

  for (i in seq_len(iterations - 1)) {
    res2 <- fsvsample(sim$y, factors = factors, draws = 1, burnin = 0,
		      keeptime = "all", quiet = TRUE, signident = FALSE,
		      startfacloadvar = res2$latestauxiliary$facloadvar,
		      startfacload = res2$facload[,,1],
		      startfac = res2$fac[,,1],
		      startlogvar = res2$logvar[,,1],
		      startlogvar0 = res2$logvar0[,1],
		      startpara = res2$para[,,1])
  }
})
print(tim2)


# In principle, res and res2 should contain exactly the same (last) draws.
# Numerical differences may arise due to, e.g., rounding errors and their
# propagation when iterations is larger than a handful or so.
# Stochastically the results should be equivalent.

cat("Are they the same?\n")

cat("Latest auxiliaries:", all.equal(res$latestauxiliary, res2$latestauxiliary), "\n")
expect_equal(res$latestauxiliary, res2$latestauxiliary)

cat("Latest (absolute) loadings:", all.equal(abs(res$facload[,,iterations]), abs(res2$facload[,,1])), "\n")
expect_equal(abs(res$facload[,,iterations]), abs(res2$facload[,,1]))

cat("Latest paras:", all.equal(res$para[,,iterations], res2$para[,,1]), "\n")
expect_equal(res$para[,,iterations], res2$para[,,1])

# If you are interested in speed, set iterations to a few hundred at least
cat("Absolute time difference:", tim2 - tim, "\n")
cat("Relative time difference:", (tim2 - tim) / tim, "\n")
})
