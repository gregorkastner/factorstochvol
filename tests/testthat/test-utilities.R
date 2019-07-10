context("Utilities")

test_that("covmat works", {
  lapply(res, function (x) {
           expect_warning(covmat(x), NA) %>%
             expect_is("array")
  })
})

test_that("cormat works", {
  lapply(res, function (x) {
           expect_warning(cormat(x), NA) %>%
             expect_is("array")
  })
})

test_that("print works", {
  lapply(res, function (x) {
           expect_output(print(x), NULL)
  })
})

test_that("predh works", {
  lapply(res, function (x) {
           expect_warning(predh(x), NA)
           expect_warning(predh(x, ahead = 3), NA)
           expect_warning(predh(x, ahead = c(1, 4), each = 4), NA)
  })
})

test_that("predcov works", {
  lapply(res, function (x) {
           expect_warning(predcov(x), NA)
           expect_warning(predcov(x, ahead = 3), NA)
           expect_warning(predcov(x, ahead = c(1, 4), each = 4), NA)
  })
})

test_that("predcor works", {
  lapply(res, function (x) {
           expect_warning(predcor(x), NA)
           expect_warning(predcor(x, ahead = 3), NA)
           expect_warning(predcor(x, ahead = c(1, 4), each = 4), NA)
  })
})

test_that("predprecWB works", {
  lapply(res, function (x) {
           expect_warning(predprecWB(x), NA)
           expect_warning(predprecWB(x, ahead = 3), NA)
           expect_warning(predprecWB(x, ahead = c(1, 4), each = 4), NA)
  })
})

test_that("predloglik works", {
  lapply(res, function (x) {
           expect_warning(predloglik(x, matrix(0,   1, NCOL(y))), NA)
           expect_warning(predloglik(x, matrix(0.1, 1, NCOL(y)), ahead = 3), NA)
           expect_warning(predloglik(x, matrix(-1,  2, NCOL(y)), ahead = c(1, 4), each = 4), NA)
  })
})

test_that("predloglikWB works", {
  lapply(res, function (x) {
           expect_warning(predloglikWB(x, matrix(0,   1, NCOL(y))), NA)
           expect_warning(predloglikWB(x, matrix(0.1, 1, NCOL(y)), ahead = 3), NA)
           expect_warning(predloglikWB(x, matrix(-1,  2, NCOL(y)), ahead = c(1, 4), each = 4), NA)
  })
})

test_that("signident works", {
  possiblemethods <- c("diagonal", "maximin")
  implementation <- 1:3
  lapply(res, function (x, pms, ims) {
           for (method in pms) {
             for (im in ims) {
               expect_warning(signident(x, method = method, implementation = im), NA)
             }
           }
  }, possiblemethods, implementation)
})

test_that("orderident works", {
  possiblemethods <- c("summean", "summeaninv", "summeanabs", "summedabs", "summed",
                       "summedinv", "maxmed", "maxmedinv", "maxmedrel", "maxmedabsrel")
  lapply(res, function (x, pms) {
           for (method in pms) {
             expect_warning(orderident(x, method = method), NA)
           }
  }, possiblemethods)
})
